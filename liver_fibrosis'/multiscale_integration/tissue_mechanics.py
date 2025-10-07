"""
Tissue Mechanics Module

Implements continuum mechanics models for tissue-level mechanical behavior
including cell density fields, stress calculations, and mechanical signal
transduction in biological systems.

Key Features:
- Continuum mechanics modeling with elastic-viscous behavior
- Cell density field evolution and dynamics
- Stress field calculations and mechanical equilibrium
- Mechanical signal transduction pathways
- Cell migration mechanics and chemomechanical coupling
- Tissue remodeling and adaptation mechanisms
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable, Union
from dataclasses import dataclass
from enum import Enum
import warnings
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.ndimage import gaussian_filter


class MaterialModel(Enum):
    """Types of material models."""
    LINEAR_ELASTIC = "linear_elastic"
    HYPERELASTIC = "hyperelastic"
    VISCOELASTIC = "viscoelastic"
    POROELASTIC = "poroelastic"


class BoundaryCondition(Enum):
    """Types of mechanical boundary conditions."""
    FIXED_DISPLACEMENT = "fixed_displacement"
    FIXED_FORCE = "fixed_force"
    FREE_SURFACE = "free_surface"
    SYMMETRY = "symmetry"


@dataclass
class MaterialProperties:
    """Material properties for tissue mechanics."""
    youngs_modulus: float = 1000.0      # Pa
    poissons_ratio: float = 0.45        # dimensionless
    bulk_modulus: float = None          # Pa
    shear_modulus: float = None         # Pa
    viscosity: float = 100.0            # Pa·s
    density: float = 1000.0             # kg/m³
    
    def __post_init__(self):
        """Calculate derived properties."""
        if self.bulk_modulus is None:
            self.bulk_modulus = self.youngs_modulus / (3 * (1 - 2 * self.poissons_ratio))
        if self.shear_modulus is None:
            self.shear_modulus = self.youngs_modulus / (2 * (1 + self.poissons_ratio))


@dataclass
class MechanicalState:
    """Container for mechanical state variables."""
    displacement: np.ndarray = None
    velocity: np.ndarray = None
    stress: np.ndarray = None
    strain: np.ndarray = None
    pressure: float = 0.0
    
    def __post_init__(self):
        """Initialize arrays if not provided."""
        if self.displacement is None:
            self.displacement = np.zeros(2)  # Default 2D
        if self.velocity is None:
            self.velocity = np.zeros_like(self.displacement)


class ContinuumMechanicsModel:
    """
    Continuum mechanics model for tissue-level mechanical behavior.
    """
    
    def __init__(self, grid_shape: Tuple[int, ...],
                 grid_spacing: Union[float, Tuple[float, ...]],
                 material_properties: MaterialProperties):
        """
        Initialize continuum mechanics model.
        
        Args:
            grid_shape: Shape of computational grid
            grid_spacing: Spatial grid spacing
            material_properties: Material properties
        """
        self.grid_shape = grid_shape
        self.ndim = len(grid_shape)
        
        if isinstance(grid_spacing, (int, float)):
            self.grid_spacing = np.full(self.ndim, grid_spacing)
        else:
            self.grid_spacing = np.array(grid_spacing)
        
        self.material_props = material_properties
        
        # Initialize mechanical fields
        self.displacement_field = np.zeros(grid_shape + (self.ndim,))
        self.velocity_field = np.zeros(grid_shape + (self.ndim,))
        self.stress_field = np.zeros(grid_shape + (self.ndim, self.ndim))
        self.strain_field = np.zeros(grid_shape + (self.ndim, self.ndim))
        
        # Initialize finite element matrices
        self._setup_finite_element_matrices()
    
    def _setup_finite_element_matrices(self):
        """Setup finite element stiffness and mass matrices."""
        # This is a simplified implementation
        # In practice, proper finite element assembly would be used
        
        total_nodes = np.prod(self.grid_shape)
        total_dofs = total_nodes * self.ndim
        
        # Simplified stiffness matrix (identity for now)
        self.stiffness_matrix = csr_matrix((total_dofs, total_dofs))
        self.mass_matrix = csr_matrix((total_dofs, total_dofs))
        
        # Build simplified matrices
        E = self.material_props.youngs_modulus
        nu = self.material_props.poissons_ratio
        
        # Plane stress/strain material matrix
        if self.ndim == 2:
            D_factor = E / (1 - nu**2)
            D_matrix = D_factor * np.array([
                [1, nu, 0],
                [nu, 1, 0],
                [0, 0, (1-nu)/2]
            ])
        elif self.ndim == 3:
            D_factor = E / ((1 + nu) * (1 - 2*nu))
            D_matrix = D_factor * np.array([
                [1-nu, nu, nu, 0, 0, 0],
                [nu, 1-nu, nu, 0, 0, 0],
                [nu, nu, 1-nu, 0, 0, 0],
                [0, 0, 0, (1-2*nu)/2, 0, 0],
                [0, 0, 0, 0, (1-2*nu)/2, 0],
                [0, 0, 0, 0, 0, (1-2*nu)/2]
            ])
        
        self.material_matrix = D_matrix
    
    def calculate_strain_from_displacement(self, displacement_field: np.ndarray) -> np.ndarray:
        """
        Calculate strain tensor from displacement field.
        
        Args:
            displacement_field: Displacement field
            
        Returns:
            Strain tensor field
        """
        strain_field = np.zeros(self.grid_shape + (self.ndim, self.ndim))
        
        # Calculate displacement gradients
        for i in range(self.ndim):
            displacement_component = displacement_field[..., i]
            
            for j in range(self.ndim):
                # Calculate partial derivative ∂u_i/∂x_j
                if j == 0:  # x-direction
                    grad = np.gradient(displacement_component, self.grid_spacing[0], axis=0)
                elif j == 1:  # y-direction
                    grad = np.gradient(displacement_component, self.grid_spacing[1], axis=1)
                elif j == 2 and self.ndim == 3:  # z-direction
                    grad = np.gradient(displacement_component, self.grid_spacing[2], axis=2)
                else:
                    continue
                
                # Strain tensor: ε_ij = 1/2 * (∂u_i/∂x_j + ∂u_j/∂x_i)
                strain_field[..., i, j] += 0.5 * grad
                strain_field[..., j, i] += 0.5 * grad
        
        return strain_field
    
    def calculate_stress_from_strain(self, strain_field: np.ndarray) -> np.ndarray:
        """
        Calculate stress tensor from strain tensor using constitutive law.
        
        Args:
            strain_field: Strain tensor field
            
        Returns:
            Stress tensor field
        """
        stress_field = np.zeros_like(strain_field)
        
        # Apply constitutive law at each grid point
        for idx in np.ndindex(self.grid_shape):
            strain_tensor = strain_field[idx]
            
            if self.ndim == 2:
                # Convert strain tensor to Voigt notation
                strain_voigt = np.array([
                    strain_tensor[0, 0],  # ε_xx
                    strain_tensor[1, 1],  # ε_yy
                    2 * strain_tensor[0, 1]  # γ_xy = 2*ε_xy
                ])
                
                # Apply material law: σ = D * ε
                stress_voigt = self.material_matrix @ strain_voigt
                
                # Convert back to tensor form
                stress_tensor = np.array([
                    [stress_voigt[0], stress_voigt[2]/2],
                    [stress_voigt[2]/2, stress_voigt[1]]
                ])
                
            elif self.ndim == 3:
                # Convert strain tensor to Voigt notation
                strain_voigt = np.array([
                    strain_tensor[0, 0],      # ε_xx
                    strain_tensor[1, 1],      # ε_yy
                    strain_tensor[2, 2],      # ε_zz
                    2 * strain_tensor[1, 2],  # γ_yz
                    2 * strain_tensor[0, 2],  # γ_xz
                    2 * strain_tensor[0, 1]   # γ_xy
                ])
                
                # Apply material law
                stress_voigt = self.material_matrix @ strain_voigt
                
                # Convert back to tensor form
                stress_tensor = np.array([
                    [stress_voigt[0], stress_voigt[5]/2, stress_voigt[4]/2],
                    [stress_voigt[5]/2, stress_voigt[1], stress_voigt[3]/2],
                    [stress_voigt[4]/2, stress_voigt[3]/2, stress_voigt[2]]
                ])
            
            stress_field[idx] = stress_tensor
        
        return stress_field
    
    def solve_equilibrium(self, body_forces: np.ndarray,
                         boundary_conditions: Dict[str, Tuple[BoundaryCondition, np.ndarray]],
                         tolerance: float = 1e-6,
                         max_iterations: int = 100) -> np.ndarray:
        """
        Solve mechanical equilibrium: ∇·σ + f = 0
        
        Args:
            body_forces: Body force field
            boundary_conditions: Boundary conditions
            tolerance: Convergence tolerance
            max_iterations: Maximum iterations
            
        Returns:
            Equilibrium displacement field
        """
        # This is a simplified iterative solver
        # In practice, proper finite element methods would be used
        
        displacement = self.displacement_field.copy()
        
        for iteration in range(max_iterations):
            # Calculate strain and stress
            strain = self.calculate_strain_from_displacement(displacement)
            stress = self.calculate_stress_from_strain(strain)
            
            # Calculate stress divergence (force balance)
            stress_divergence = self._calculate_stress_divergence(stress)
            
            # Calculate residual: R = ∇·σ + f
            residual = stress_divergence + body_forces
            
            # Check convergence
            residual_norm = np.linalg.norm(residual)
            if residual_norm < tolerance:
                break
            
            # Update displacement (simplified Newton-Raphson)
            displacement_correction = -0.01 * residual  # Simplified update
            displacement += displacement_correction
            
            # Apply boundary conditions
            displacement = self._apply_mechanical_boundary_conditions(
                displacement, boundary_conditions)
        
        # Update fields
        self.displacement_field = displacement
        self.strain_field = self.calculate_strain_from_displacement(displacement)
        self.stress_field = self.calculate_stress_from_strain(self.strain_field)
        
        return displacement
    
    def _calculate_stress_divergence(self, stress_field: np.ndarray) -> np.ndarray:
        """Calculate divergence of stress tensor."""
        divergence = np.zeros(self.grid_shape + (self.ndim,))
        
        for i in range(self.ndim):
            for j in range(self.ndim):
                stress_component = stress_field[..., i, j]
                
                if j == 0:  # ∂σ_ij/∂x
                    grad = np.gradient(stress_component, self.grid_spacing[0], axis=0)
                elif j == 1:  # ∂σ_ij/∂y
                    grad = np.gradient(stress_component, self.grid_spacing[1], axis=1)
                elif j == 2 and self.ndim == 3:  # ∂σ_ij/∂z
                    grad = np.gradient(stress_component, self.grid_spacing[2], axis=2)
                else:
                    continue
                
                divergence[..., i] += grad
        
        return divergence
    
    def _apply_mechanical_boundary_conditions(self, displacement: np.ndarray,
                                            boundary_conditions: Dict) -> np.ndarray:
        """Apply mechanical boundary conditions."""
        updated_displacement = displacement.copy()
        
        # This is a simplified implementation
        # Proper boundary condition application would depend on the specific
        # finite element discretization
        
        for bc_name, (bc_type, bc_value) in boundary_conditions.items():
            if bc_type == BoundaryCondition.FIXED_DISPLACEMENT:
                # Apply fixed displacement at boundaries
                if bc_name == "left" and self.ndim >= 1:
                    updated_displacement[0, :] = bc_value
                elif bc_name == "right" and self.ndim >= 1:
                    updated_displacement[-1, :] = bc_value
                elif bc_name == "bottom" and self.ndim >= 2:
                    updated_displacement[:, 0] = bc_value
                elif bc_name == "top" and self.ndim >= 2:
                    updated_displacement[:, -1] = bc_value
        
        return updated_displacement


class CellDensityFieldModel:
    """
    Model for cell density field evolution and dynamics.
    """
    
    def __init__(self, grid_shape: Tuple[int, ...],
                 grid_spacing: Union[float, Tuple[float, ...]]):
        """Initialize cell density field model."""
        self.grid_shape = grid_shape
        self.ndim = len(grid_shape)
        
        if isinstance(grid_spacing, (int, float)):
            self.grid_spacing = np.full(self.ndim, grid_spacing)
        else:
            self.grid_spacing = np.array(grid_spacing)
        
        # Initialize cell density field
        self.cell_density = np.zeros(grid_shape)
        self.cell_velocity = np.zeros(grid_shape + (self.ndim,))
        
        # Model parameters
        self.diffusion_coefficient = 1e-12  # m²/s
        self.proliferation_rate = 1e-5      # 1/s
        self.death_rate = 1e-6              # 1/s
        self.carrying_capacity = 1000.0     # cells/m³
    
    def update_cell_density(self, dt: float,
                           chemical_gradients: Optional[Dict[str, np.ndarray]] = None,
                           mechanical_stress: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Update cell density field.
        
        Args:
            dt: Time step
            chemical_gradients: Chemical gradient fields for chemotaxis
            mechanical_stress: Mechanical stress field
            
        Returns:
            Updated cell density field
        """
        # Calculate cell migration velocity
        migration_velocity = self._calculate_cell_migration_velocity(
            chemical_gradients, mechanical_stress)
        
        # Calculate diffusion term
        diffusion_term = self._calculate_cell_diffusion()
        
        # Calculate proliferation and death terms
        proliferation_term = self._calculate_proliferation()
        death_term = self._calculate_death()
        
        # Calculate advection term (cell transport)
        advection_term = self._calculate_advection(migration_velocity)
        
        # Update cell density: ∂ρ/∂t = D∇²ρ - ∇·(ρv) + P - D
        density_change = (diffusion_term - advection_term + 
                         proliferation_term - death_term) * dt
        
        self.cell_density = np.maximum(0.0, self.cell_density + density_change)
        self.cell_velocity = migration_velocity
        
        return self.cell_density
    
    def _calculate_cell_migration_velocity(self, chemical_gradients: Optional[Dict[str, np.ndarray]],
                                         mechanical_stress: Optional[np.ndarray]) -> np.ndarray:
        """Calculate cell migration velocity from various stimuli."""
        velocity = np.zeros(self.grid_shape + (self.ndim,))
        
        # Chemotaxis component
        if chemical_gradients is not None:
            for species_name, gradient in chemical_gradients.items():
                # Chemotactic sensitivity (species-dependent)
                if species_name == "TGF_beta":
                    chi = 1e-10  # m²/(s·M)
                elif species_name == "PDGF":
                    chi = 2e-10
                else:
                    chi = 1e-11
                
                # Add chemotactic velocity: v = χ∇C
                velocity += chi * gradient
        
        # Mechanotaxis component (durotaxis)
        if mechanical_stress is not None:
            # Calculate stress gradient for durotaxis
            stress_magnitude = np.sqrt(np.sum(mechanical_stress**2, axis=(-2, -1)))
            
            if self.ndim == 2:
                stress_gradient = np.stack([
                    np.gradient(stress_magnitude, self.grid_spacing[0], axis=0),
                    np.gradient(stress_magnitude, self.grid_spacing[1], axis=1)
                ], axis=-1)
            elif self.ndim == 3:
                stress_gradient = np.stack([
                    np.gradient(stress_magnitude, self.grid_spacing[0], axis=0),
                    np.gradient(stress_magnitude, self.grid_spacing[1], axis=1),
                    np.gradient(stress_magnitude, self.grid_spacing[2], axis=2)
                ], axis=-1)
            
            # Durotactic sensitivity
            xi = 1e-8  # m³/(Pa·s)
            velocity += xi * stress_gradient
        
        # Random motility component
        random_velocity = np.random.normal(0, 1e-8, velocity.shape)
        velocity += random_velocity
        
        return velocity
    
    def _calculate_cell_diffusion(self) -> np.ndarray:
        """Calculate cell diffusion term."""
        # Calculate Laplacian of cell density
        if self.ndim == 2:
            d2rho_dx2 = np.zeros_like(self.cell_density)
            d2rho_dy2 = np.zeros_like(self.cell_density)
            
            d2rho_dx2[1:-1, :] = ((self.cell_density[2:, :] - 2*self.cell_density[1:-1, :] + 
                                  self.cell_density[:-2, :]) / self.grid_spacing[0]**2)
            d2rho_dy2[:, 1:-1] = ((self.cell_density[:, 2:] - 2*self.cell_density[:, 1:-1] + 
                                  self.cell_density[:, :-2]) / self.grid_spacing[1]**2)
            
            laplacian = d2rho_dx2 + d2rho_dy2
            
        elif self.ndim == 3:
            d2rho_dx2 = np.zeros_like(self.cell_density)
            d2rho_dy2 = np.zeros_like(self.cell_density)
            d2rho_dz2 = np.zeros_like(self.cell_density)
            
            d2rho_dx2[1:-1, :, :] = ((self.cell_density[2:, :, :] - 2*self.cell_density[1:-1, :, :] + 
                                     self.cell_density[:-2, :, :]) / self.grid_spacing[0]**2)
            d2rho_dy2[:, 1:-1, :] = ((self.cell_density[:, 2:, :] - 2*self.cell_density[:, 1:-1, :] + 
                                     self.cell_density[:, :-2, :]) / self.grid_spacing[1]**2)
            d2rho_dz2[:, :, 1:-1] = ((self.cell_density[:, :, 2:] - 2*self.cell_density[:, :, 1:-1] + 
                                     self.cell_density[:, :, :-2]) / self.grid_spacing[2]**2)
            
            laplacian = d2rho_dx2 + d2rho_dy2 + d2rho_dz2
        
        return self.diffusion_coefficient * laplacian
    
    def _calculate_proliferation(self) -> np.ndarray:
        """Calculate cell proliferation term."""
        # Logistic growth with carrying capacity
        growth_factor = 1.0 - self.cell_density / self.carrying_capacity
        proliferation = self.proliferation_rate * self.cell_density * np.maximum(0.0, growth_factor)
        return proliferation
    
    def _calculate_death(self) -> np.ndarray:
        """Calculate cell death term."""
        # Simple death rate
        death = self.death_rate * self.cell_density
        return death
    
    def _calculate_advection(self, velocity: np.ndarray) -> np.ndarray:
        """Calculate advection term ∇·(ρv)."""
        advection = np.zeros_like(self.cell_density)
        
        # Calculate divergence of (ρv)
        for i in range(self.ndim):
            flux_component = self.cell_density * velocity[..., i]
            
            if i == 0:  # x-direction
                flux_gradient = np.gradient(flux_component, self.grid_spacing[0], axis=0)
            elif i == 1:  # y-direction
                flux_gradient = np.gradient(flux_component, self.grid_spacing[1], axis=1)
            elif i == 2:  # z-direction
                flux_gradient = np.gradient(flux_component, self.grid_spacing[2], axis=2)
            
            advection += flux_gradient
        
        return advection


class StressFieldCalculator:
    """
    Calculator for stress fields and mechanical quantities.
    """
    
    def __init__(self, material_properties: MaterialProperties):
        """Initialize stress field calculator."""
        self.material_props = material_properties
    
    def calculate_von_mises_stress(self, stress_tensor: np.ndarray) -> np.ndarray:
        """
        Calculate von Mises stress from stress tensor.
        
        Args:
            stress_tensor: Stress tensor field
            
        Returns:
            von Mises stress field
        """
        if stress_tensor.shape[-1] == 2:  # 2D
            sigma_xx = stress_tensor[..., 0, 0]
            sigma_yy = stress_tensor[..., 1, 1]
            sigma_xy = stress_tensor[..., 0, 1]
            
            von_mises = np.sqrt(sigma_xx**2 + sigma_yy**2 - sigma_xx*sigma_yy + 3*sigma_xy**2)
            
        elif stress_tensor.shape[-1] == 3:  # 3D
            sigma_xx = stress_tensor[..., 0, 0]
            sigma_yy = stress_tensor[..., 1, 1]
            sigma_zz = stress_tensor[..., 2, 2]
            sigma_xy = stress_tensor[..., 0, 1]
            sigma_xz = stress_tensor[..., 0, 2]
            sigma_yz = stress_tensor[..., 1, 2]
            
            von_mises = np.sqrt(0.5 * ((sigma_xx - sigma_yy)**2 + 
                                      (sigma_yy - sigma_zz)**2 + 
                                      (sigma_zz - sigma_xx)**2 + 
                                      6 * (sigma_xy**2 + sigma_xz**2 + sigma_yz**2)))
        
        return von_mises
    
    def calculate_hydrostatic_pressure(self, stress_tensor: np.ndarray) -> np.ndarray:
        """Calculate hydrostatic pressure (negative mean stress)."""
        trace = np.trace(stress_tensor, axis1=-2, axis2=-1)
        pressure = -trace / stress_tensor.shape[-1]
        return pressure
    
    def calculate_deviatoric_stress(self, stress_tensor: np.ndarray) -> np.ndarray:
        """Calculate deviatoric stress tensor."""
        pressure = self.calculate_hydrostatic_pressure(stress_tensor)
        ndim = stress_tensor.shape[-1]
        
        # Create identity tensor
        identity = np.eye(ndim)
        identity_field = np.broadcast_to(identity, stress_tensor.shape)
        
        # Deviatoric stress = total stress - hydrostatic part
        pressure_field = pressure[..., np.newaxis, np.newaxis] * identity_field
        deviatoric_stress = stress_tensor - pressure_field
        
        return deviatoric_stress


class MechanicalSignalTransduction:
    """
    Model for mechanical signal transduction pathways.
    """
    
    def __init__(self):
        """Initialize mechanical signal transduction model."""
        self.mechanosensitive_channels = {}
        self.focal_adhesion_proteins = {}
        self.cytoskeletal_tension = {}
    
    def calculate_mechanosensitive_channel_activity(self, membrane_tension: np.ndarray,
                                                  channel_sensitivity: float = 1.0) -> np.ndarray:
        """
        Calculate mechanosensitive channel activity.
        
        Args:
            membrane_tension: Membrane tension field
            channel_sensitivity: Channel sensitivity parameter
            
        Returns:
            Channel activity (0-1)
        """
        # Sigmoid activation function
        activity = 1.0 / (1.0 + np.exp(-channel_sensitivity * (membrane_tension - 0.5)))
        return activity
    
    def calculate_focal_adhesion_maturation(self, traction_stress: np.ndarray,
                                          substrate_stiffness: np.ndarray) -> np.ndarray:
        """
        Calculate focal adhesion maturation based on mechanical forces.
        
        Args:
            traction_stress: Cellular traction stress
            substrate_stiffness: Substrate stiffness field
            
        Returns:
            Focal adhesion maturation level
        """
        # Focal adhesion maturation depends on both force and substrate stiffness
        optimal_stress = 1000.0  # Pa
        optimal_stiffness = 10000.0  # Pa
        
        stress_factor = traction_stress / (optimal_stress + traction_stress)
        stiffness_factor = substrate_stiffness / (optimal_stiffness + substrate_stiffness)
        
        maturation = stress_factor * stiffness_factor
        return maturation
    
    def calculate_cytoskeletal_remodeling(self, mechanical_stress: np.ndarray,
                                        current_fiber_density: np.ndarray,
                                        dt: float) -> np.ndarray:
        """
        Calculate cytoskeletal remodeling in response to mechanical stress.
        
        Args:
            mechanical_stress: Applied mechanical stress
            current_fiber_density: Current stress fiber density
            dt: Time step
            
        Returns:
            Updated fiber density
        """
        # Stress fiber formation rate depends on applied stress
        formation_threshold = 500.0  # Pa
        formation_rate = 0.1  # 1/s
        
        stress_magnitude = np.sqrt(np.sum(mechanical_stress**2, axis=(-2, -1)))
        
        # Formation occurs above threshold
        formation_stimulus = np.maximum(0.0, stress_magnitude - formation_threshold)
        fiber_formation = formation_rate * formation_stimulus * dt
        
        # Fiber disassembly (first-order decay)
        disassembly_rate = 0.05  # 1/s
        fiber_disassembly = disassembly_rate * current_fiber_density * dt
        
        # Update fiber density
        new_fiber_density = current_fiber_density + fiber_formation - fiber_disassembly
        new_fiber_density = np.maximum(0.0, np.minimum(1.0, new_fiber_density))
        
        return new_fiber_density


class ElasticViscousModel:
    """
    Elastic-viscous material model for tissue mechanics.
    """
    
    def __init__(self, material_properties: MaterialProperties):
        """Initialize elastic-viscous model."""
        self.material_props = material_properties
    
    def calculate_stress_rate(self, strain_rate: np.ndarray,
                            current_stress: np.ndarray) -> np.ndarray:
        """
        Calculate stress rate for viscoelastic material.
        
        Args:
            strain_rate: Strain rate tensor
            current_stress: Current stress tensor
            
        Returns:
            Stress rate tensor
        """
        # Maxwell model: σ̇ = E*ε̇ - σ/τ
        # where τ = η/E is the relaxation time
        
        E = self.material_props.youngs_modulus
        eta = self.material_props.viscosity
        tau = eta / E  # Relaxation time
        
        # Elastic stress rate
        elastic_stress_rate = E * strain_rate
        
        # Viscous relaxation
        viscous_relaxation = current_stress / tau
        
        stress_rate = elastic_stress_rate - viscous_relaxation
        
        return stress_rate
    
    def update_stress_viscoelastic(self, current_stress: np.ndarray,
                                 strain_rate: np.ndarray,
                                 dt: float) -> np.ndarray:
        """
        Update stress using viscoelastic constitutive law.
        
        Args:
            current_stress: Current stress tensor
            strain_rate: Strain rate tensor
            dt: Time step
            
        Returns:
            Updated stress tensor
        """
        stress_rate = self.calculate_stress_rate(strain_rate, current_stress)
        new_stress = current_stress + stress_rate * dt
        
        return new_stress


class CellMigrationMechanics:
    """
    Mechanics of cell migration including traction forces and protrusion dynamics.
    """
    
    def __init__(self, grid_shape: Tuple[int, ...]):
        """Initialize cell migration mechanics."""
        self.grid_shape = grid_shape
        self.ndim = len(grid_shape)
        
        # Initialize traction stress field
        self.traction_stress = np.zeros(grid_shape + (self.ndim, self.ndim))
        
        # Migration parameters
        self.max_traction_stress = 1000.0  # Pa
        self.protrusion_force = 100.0      # pN
        self.adhesion_strength = 50.0      # pN/μm²
    
    def calculate_traction_forces(self, cell_positions: np.ndarray,
                                cell_velocities: np.ndarray,
                                substrate_stiffness: np.ndarray) -> np.ndarray:
        """
        Calculate cellular traction forces.
        
        Args:
            cell_positions: Cell positions
            cell_velocities: Cell velocities
            substrate_stiffness: Substrate stiffness field
            
        Returns:
            Traction force field
        """
        traction_forces = np.zeros(self.grid_shape + (self.ndim,))
        
        # For each cell, calculate traction based on velocity and substrate
        for i, (pos, vel) in enumerate(zip(cell_positions, cell_velocities)):
            # Convert position to grid indices
            grid_idx = tuple(int(pos[j] / self.grid_shape[j]) for j in range(self.ndim))
            
            # Ensure indices are within bounds
            grid_idx = tuple(max(0, min(self.grid_shape[j]-1, grid_idx[j])) 
                           for j in range(self.ndim))
            
            # Calculate traction force based on velocity and substrate stiffness
            local_stiffness = substrate_stiffness[grid_idx] if substrate_stiffness.ndim > 0 else substrate_stiffness
            traction_magnitude = min(self.max_traction_stress, 
                                   local_stiffness * np.linalg.norm(vel))
            
            # Direction opposite to velocity (resistance)
            if np.linalg.norm(vel) > 1e-12:
                traction_direction = -vel / np.linalg.norm(vel)
            else:
                traction_direction = np.zeros(self.ndim)
            
            traction_force = traction_magnitude * traction_direction
            
            # Add to traction force field
            traction_forces[grid_idx] += traction_force
        
        return traction_forces
    
    def calculate_protrusion_dynamics(self, cell_positions: np.ndarray,
                                    chemical_gradients: Dict[str, np.ndarray],
                                    membrane_tension: np.ndarray) -> np.ndarray:
        """
        Calculate cell protrusion dynamics.
        
        Args:
            cell_positions: Cell positions
            chemical_gradients: Chemical gradient fields
            membrane_tension: Membrane tension field
            
        Returns:
            Protrusion force field
        """
        protrusion_forces = np.zeros(self.grid_shape + (self.ndim,))
        
        for i, pos in enumerate(cell_positions):
            # Convert position to grid indices
            grid_idx = tuple(int(pos[j] / self.grid_shape[j]) for j in range(self.ndim))
            grid_idx = tuple(max(0, min(self.grid_shape[j]-1, grid_idx[j])) 
                           for j in range(self.ndim))
            
            # Calculate protrusion force based on chemical gradients
            total_gradient = np.zeros(self.ndim)
            for species_name, gradient in chemical_gradients.items():
                if gradient.ndim == self.ndim + 1:  # Gradient field
                    local_gradient = gradient[grid_idx]
                    total_gradient += local_gradient
            
            # Protrusion force proportional to gradient
            protrusion_force = self.protrusion_force * total_gradient
            
            # Modulate by membrane tension (higher tension reduces protrusion)
            if membrane_tension.size > 1:
                local_tension = membrane_tension[grid_idx]
                tension_factor = 1.0 / (1.0 + local_tension)
                protrusion_force *= tension_factor
            
            protrusion_forces[grid_idx] += protrusion_force
        
        return protrusion_forces
