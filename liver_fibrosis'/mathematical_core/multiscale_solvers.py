"""
Multi-Scale Numerical Solvers Module

Implements advanced numerical methods for solving coupled multi-scale systems
including reaction-diffusion equations, stiff ODEs, and adaptive time stepping.

Key Features:
- ADI (Alternating Direction Implicit) solver for 3D diffusion
- Runge-Kutta-Fehlberg adaptive method for stiff systems
- Operator splitting for coupled reaction-diffusion systems
- Multiple time stepping for multi-scale problems
- Embedded error estimation and adaptive control
"""

import numpy as np
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve
from typing import Callable, Tuple, Optional, Dict, Any
import warnings


class ADISolver:
    """
    Alternating Direction Implicit solver for 3D diffusion equations.
    
    Solves: ∂C/∂t = D(∇²C) + R(C) using operator splitting:
    - Step 1: ∂C/∂t = D(∂²C/∂x²) + D(∂²C/∂y²) + D(∂²C/∂z²)
    - Step 2: Apply reaction terms
    """
    
    def __init__(self, grid_shape: Tuple[int, int, int], 
                 grid_spacing: Tuple[float, float, float],
                 boundary_conditions: str = "neumann"):
        """
        Initialize ADI solver.
        
        Args:
            grid_shape: (nx, ny, nz) grid dimensions
            grid_spacing: (dx, dy, dz) spatial resolution
            boundary_conditions: "neumann", "dirichlet", or "periodic"
        """
        self.nx, self.ny, self.nz = grid_shape
        self.dx, self.dy, self.dz = grid_spacing
        self.bc_type = boundary_conditions
        
        # Pre-compute finite difference matrices
        self._setup_difference_matrices()
        
    def _setup_difference_matrices(self):
        """Setup finite difference matrices for each direction."""
        # X-direction matrix
        diag_x = np.ones(self.nx)
        off_diag_x = np.ones(self.nx - 1)
        
        if self.bc_type == "neumann":
            # Neumann boundary conditions: zero flux
            self.Ax = diags([-off_diag_x, 2*diag_x, -off_diag_x], 
                           [-1, 0, 1], shape=(self.nx, self.nx))
            self.Ax[0, 1] = -2  # Adjust for boundary
            self.Ax[-1, -2] = -2
        elif self.bc_type == "dirichlet":
            # Dirichlet boundary conditions: fixed values
            self.Ax = diags([-off_diag_x, 2*diag_x, -off_diag_x],
                           [-1, 0, 1], shape=(self.nx, self.nx))
        else:
            raise ValueError(f"Unsupported boundary condition: {self.bc_type}")
            
        # Y and Z direction matrices (similar structure)
        self.Ay = diags([-np.ones(self.ny-1), 2*np.ones(self.ny), -np.ones(self.ny-1)],
                       [-1, 0, 1], shape=(self.ny, self.ny))
        self.Az = diags([-np.ones(self.nz-1), 2*np.ones(self.nz), -np.ones(self.nz-1)], 
                       [-1, 0, 1], shape=(self.nz, self.nz))
        
        # Convert to CSC format for efficient solving
        self.Ax = csc_matrix(self.Ax)
        self.Ay = csc_matrix(self.Ay)
        self.Az = csc_matrix(self.Az)
        
    def solve_step(self, concentration: np.ndarray, diffusion_coeff: float, 
                   dt: float) -> np.ndarray:
        """
        Solve one ADI time step.
        
        Args:
            concentration: Current concentration field (nx, ny, nz)
            diffusion_coeff: Diffusion coefficient
            dt: Time step size
            
        Returns:
            Updated concentration field
        """
        C = concentration.copy()
        
        # ADI parameters
        r_x = diffusion_coeff * dt / (2 * self.dx**2)
        r_y = diffusion_coeff * dt / (2 * self.dy**2) 
        r_z = diffusion_coeff * dt / (2 * self.dz**2)
        
        # Step 1: Implicit in x, explicit in y,z
        I_x = np.eye(self.nx)
        A_x_impl = I_x + r_x * self.Ax
        
        for j in range(self.ny):
            for k in range(self.nz):
                rhs = C[:, j, k] - r_y * (self.Ay @ C[:, j, :]).sum(axis=1) - \
                      r_z * (self.Az @ C[:, :, k]).sum(axis=0)
                C[:, j, k] = spsolve(A_x_impl, rhs)
        
        # Step 2: Implicit in y, explicit in x,z  
        I_y = np.eye(self.ny)
        A_y_impl = I_y + r_y * self.Ay
        
        for i in range(self.nx):
            for k in range(self.nz):
                rhs = C[i, :, k] - r_x * (self.Ax @ C[:, :, k]).sum(axis=0) - \
                      r_z * (self.Az @ C[i, :, :]).sum(axis=1)
                C[i, :, k] = spsolve(A_y_impl, rhs)
                
        # Step 3: Implicit in z, explicit in x,y
        I_z = np.eye(self.nz)
        A_z_impl = I_z + r_z * self.Az
        
        for i in range(self.nx):
            for j in range(self.ny):
                rhs = C[i, j, :] - r_x * (self.Ax @ C[:, j, :]).sum(axis=0) - \
                      r_y * (self.Ay @ C[i, :, :]).sum(axis=1)
                C[i, j, :] = spsolve(A_z_impl, rhs)
                
        return C


class RungeKuttaFehlbergSolver:
    """
    Runge-Kutta-Fehlberg adaptive solver for stiff ODE systems.
    
    Implements RKF45 method with embedded error estimation for
    automatic step size control.
    """
    
    def __init__(self, rtol: float = 1e-6, atol: float = 1e-9,
                 max_step: float = np.inf, min_step: float = 1e-12):
        """
        Initialize RKF solver.
        
        Args:
            rtol: Relative tolerance
            atol: Absolute tolerance  
            max_step: Maximum step size
            min_step: Minimum step size
        """
        self.rtol = rtol
        self.atol = atol
        self.max_step = max_step
        self.min_step = min_step
        
        # RKF45 coefficients
        self.a = np.array([
            [0, 0, 0, 0, 0, 0],
            [1/4, 0, 0, 0, 0, 0],
            [3/32, 9/32, 0, 0, 0, 0],
            [1932/2197, -7200/2197, 7296/2197, 0, 0, 0],
            [439/216, -8, 3680/513, -845/4104, 0, 0],
            [-8/27, 2, -3544/2565, 1859/4104, -11/40, 0]
        ])
        
        self.b4 = np.array([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])
        self.b5 = np.array([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])
        
    def solve_step(self, f: Callable, t: float, y: np.ndarray, 
                   h: float) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Solve one RKF step with error estimation.
        
        Args:
            f: Right-hand side function dy/dt = f(t, y)
            t: Current time
            y: Current state vector
            h: Step size
            
        Returns:
            y_new: Updated state (4th order)
            y_err: Error estimate
            h_new: Suggested new step size
        """
        # Compute k values
        k = np.zeros((6, len(y)))
        k[0] = h * f(t, y)
        
        for i in range(1, 6):
            y_temp = y + np.sum(self.a[i, :i, np.newaxis] * k[:i], axis=0)
            k[i] = h * f(t + np.sum(self.a[i, :i]) * h, y_temp)
        
        # 4th and 5th order solutions
        y4 = y + np.sum(self.b4[:, np.newaxis] * k, axis=0)
        y5 = y + np.sum(self.b5[:, np.newaxis] * k, axis=0)
        
        # Error estimate
        error = np.abs(y5 - y4)
        
        # Compute error tolerance
        tol = self.atol + self.rtol * np.maximum(np.abs(y), np.abs(y4))
        
        # Error ratio
        error_ratio = np.max(error / tol)
        
        # Adaptive step size
        if error_ratio <= 1.0:
            # Accept step
            h_new = min(self.max_step, 0.9 * h * (1.0 / error_ratio)**0.2)
        else:
            # Reject step
            h_new = max(self.min_step, 0.9 * h * (1.0 / error_ratio)**0.25)
            
        return y4, error, h_new
    
    def integrate(self, f: Callable, t_span: Tuple[float, float], 
                  y0: np.ndarray, h0: float = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Integrate ODE system over time interval.
        
        Args:
            f: Right-hand side function
            t_span: (t_start, t_end) integration interval
            y0: Initial conditions
            h0: Initial step size (auto if None)
            
        Returns:
            t_values: Time points
            y_values: Solution values
        """
        t_start, t_end = t_span
        
        if h0 is None:
            h0 = (t_end - t_start) / 100
            
        t = t_start
        y = y0.copy()
        h = h0
        
        t_values = [t]
        y_values = [y.copy()]
        
        while t < t_end:
            if t + h > t_end:
                h = t_end - t
                
            y_new, error, h_new = self.solve_step(f, t, y, h)
            
            # Check if step is acceptable
            tol = self.atol + self.rtol * np.maximum(np.abs(y), np.abs(y_new))
            error_ratio = np.max(error / tol)
            
            if error_ratio <= 1.0:
                # Accept step
                t += h
                y = y_new
                t_values.append(t)
                y_values.append(y.copy())
                
            h = h_new
            
        return np.array(t_values), np.array(y_values)


class OperatorSplittingSolver:
    """
    Operator splitting solver for coupled reaction-diffusion systems.
    
    Splits the equation ∂C/∂t = D∇²C + R(C) into:
    1. Diffusion step: ∂C/∂t = D∇²C  
    2. Reaction step: ∂C/∂t = R(C)
    """
    
    def __init__(self, diffusion_solver: ADISolver, 
                 reaction_solver: RungeKuttaFehlbergSolver):
        """
        Initialize operator splitting solver.
        
        Args:
            diffusion_solver: Solver for diffusion step
            reaction_solver: Solver for reaction step
        """
        self.diffusion_solver = diffusion_solver
        self.reaction_solver = reaction_solver
        
    def solve_step(self, concentration: np.ndarray, 
                   diffusion_coeffs: Dict[str, float],
                   reaction_func: Callable,
                   dt: float, 
                   species_names: list) -> np.ndarray:
        """
        Solve one operator splitting step.
        
        Args:
            concentration: Current concentration fields
            diffusion_coeffs: Diffusion coefficients for each species
            reaction_func: Reaction function R(C)
            dt: Time step
            species_names: Names of chemical species
            
        Returns:
            Updated concentration fields
        """
        C = concentration.copy()
        
        # Step 1: Diffusion (half time step)
        for i, species in enumerate(species_names):
            if species in diffusion_coeffs:
                D = diffusion_coeffs[species]
                C[..., i] = self.diffusion_solver.solve_step(
                    C[..., i], D, dt/2)
        
        # Step 2: Reaction (full time step)
        shape = C.shape
        C_flat = C.reshape(-1, len(species_names))
        
        for idx in range(C_flat.shape[0]):
            def reaction_ode(t, y):
                return reaction_func(y)
            
            t_vals, y_vals = self.reaction_solver.integrate(
                reaction_ode, (0, dt), C_flat[idx], dt/10)
            C_flat[idx] = y_vals[-1]
            
        C = C_flat.reshape(shape)
        
        # Step 3: Diffusion (half time step)
        for i, species in enumerate(species_names):
            if species in diffusion_coeffs:
                D = diffusion_coeffs[species]
                C[..., i] = self.diffusion_solver.solve_step(
                    C[..., i], D, dt/2)
                
        return C


class MultipleTimeSteppingSolver:
    """
    Multiple time stepping solver for multi-scale problems.
    
    Uses different time steps for processes with different characteristic
    time scales (diffusion, mechanics, biochemistry).
    """
    
    def __init__(self, dt_diffusion: float, dt_mechanics: float, 
                 dt_biochemistry: float):
        """
        Initialize multiple time stepping solver.
        
        Args:
            dt_diffusion: Time step for diffusion processes
            dt_mechanics: Time step for mechanical processes  
            dt_biochemistry: Time step for biochemical processes
        """
        self.dt_diff = dt_diffusion
        self.dt_mech = dt_mechanics
        self.dt_bio = dt_biochemistry
        
        # Find least common multiple for synchronization
        self.dt_sync = np.lcm.reduce([
            int(dt_diffusion * 1000),
            int(dt_mechanics * 1000), 
            int(dt_biochemistry * 1000)
        ]) / 1000.0
        
        self.n_diff = int(self.dt_sync / dt_diffusion)
        self.n_mech = int(self.dt_sync / dt_mechanics)
        self.n_bio = int(self.dt_sync / dt_biochemistry)
        
    def solve_step(self, state: Dict[str, np.ndarray],
                   diffusion_func: Callable,
                   mechanics_func: Callable, 
                   biochemistry_func: Callable) -> Dict[str, np.ndarray]:
        """
        Solve one synchronized time step.
        
        Args:
            state: Current system state
            diffusion_func: Function for diffusion updates
            mechanics_func: Function for mechanics updates
            biochemistry_func: Function for biochemistry updates
            
        Returns:
            Updated system state
        """
        new_state = {key: val.copy() for key, val in state.items()}
        
        # Sub-stepping for each process
        for i in range(max(self.n_diff, self.n_mech, self.n_bio)):
            
            # Diffusion updates
            if i < self.n_diff:
                new_state = diffusion_func(new_state, self.dt_diff)
                
            # Mechanics updates  
            if i < self.n_mech:
                new_state = mechanics_func(new_state, self.dt_mech)
                
            # Biochemistry updates
            if i < self.n_bio:
                new_state = biochemistry_func(new_state, self.dt_bio)
                
        return new_state


class AdaptiveTimeStepController:
    """
    Adaptive time step controller with embedded error estimation.
    
    Implements automatic step size adjustment based on error estimates
    and stability criteria.
    """
    
    def __init__(self, rtol: float = 1e-6, atol: float = 1e-9,
                 safety_factor: float = 0.9, max_factor: float = 2.0,
                 min_factor: float = 0.1):
        """
        Initialize adaptive controller.
        
        Args:
            rtol: Relative tolerance
            atol: Absolute tolerance
            safety_factor: Safety factor for step size adjustment
            max_factor: Maximum step size increase factor
            min_factor: Minimum step size decrease factor
        """
        self.rtol = rtol
        self.atol = atol
        self.safety = safety_factor
        self.max_factor = max_factor
        self.min_factor = min_factor
        
    def estimate_error(self, y_high: np.ndarray, y_low: np.ndarray,
                      y_current: np.ndarray) -> float:
        """
        Estimate truncation error using embedded methods.
        
        Args:
            y_high: Higher order solution
            y_low: Lower order solution  
            y_current: Current solution
            
        Returns:
            Error estimate
        """
        error = np.abs(y_high - y_low)
        tolerance = self.atol + self.rtol * np.maximum(
            np.abs(y_current), np.abs(y_high))
        
        error_ratio = np.max(error / tolerance)
        return error_ratio
        
    def adjust_step_size(self, error_ratio: float, current_dt: float,
                        order: int = 4) -> Tuple[float, bool]:
        """
        Adjust step size based on error estimate.
        
        Args:
            error_ratio: Ratio of error to tolerance
            current_dt: Current time step
            order: Order of the method
            
        Returns:
            new_dt: New time step
            accept: Whether to accept current step
        """
        if error_ratio <= 1.0:
            # Accept step
            accept = True
            if error_ratio == 0:
                factor = self.max_factor
            else:
                factor = min(self.max_factor, 
                           self.safety * (1.0 / error_ratio)**(1.0/(order+1)))
        else:
            # Reject step
            accept = False
            factor = max(self.min_factor,
                        self.safety * (1.0 / error_ratio)**(1.0/order))
            
        new_dt = factor * current_dt
        return new_dt, accept
        
    def cfl_condition(self, diffusion_coeff: float, grid_spacing: float,
                     dt: float) -> bool:
        """
        Check CFL stability condition for diffusion.
        
        Args:
            diffusion_coeff: Diffusion coefficient
            grid_spacing: Spatial grid spacing
            dt: Time step
            
        Returns:
            True if CFL condition is satisfied
        """
        cfl_number = diffusion_coeff * dt / (grid_spacing**2)
        return cfl_number <= 0.5  # Stability limit for explicit schemes
