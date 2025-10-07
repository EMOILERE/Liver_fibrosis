"""
Stochastic Processes Module

Implements advanced stochastic algorithms for biological system modeling
including MCMC sampling, stochastic chemical reactions, and random processes.

Key Features:
- Metropolis-Hastings MCMC sampling
- Gillespie algorithm for stochastic chemical kinetics
- Ornstein-Uhlenbeck process simulation
- Wiener process generation
- Stochastic differential equation solvers
"""

import numpy as np
from scipy.stats import multivariate_normal, norm
from typing import Callable, Tuple, Optional, Dict, List
import warnings


class MetropolisHastingsSampler:
    """
    Metropolis-Hastings MCMC sampler for Bayesian inference and
    parameter estimation in biological systems.
    """
    
    def __init__(self, log_posterior: Callable, proposal_cov: np.ndarray,
                 initial_state: np.ndarray):
        """
        Initialize Metropolis-Hastings sampler.
        
        Args:
            log_posterior: Log posterior probability function
            proposal_cov: Proposal distribution covariance matrix
            initial_state: Initial parameter values
        """
        self.log_posterior = log_posterior
        self.proposal_cov = proposal_cov
        self.current_state = initial_state.copy()
        self.current_log_prob = log_posterior(initial_state)
        
        self.n_accepted = 0
        self.n_total = 0
        
    def sample(self, n_samples: int, burn_in: int = 1000,
               thin: int = 1) -> Tuple[np.ndarray, float]:
        """
        Generate MCMC samples.
        
        Args:
            n_samples: Number of samples to generate
            burn_in: Number of burn-in samples
            thin: Thinning interval
            
        Returns:
            samples: Array of samples (n_samples, n_params)
            acceptance_rate: Acceptance rate
        """
        total_iterations = burn_in + n_samples * thin
        samples = []
        
        for i in range(total_iterations):
            # Propose new state
            proposal = multivariate_normal.rvs(
                mean=self.current_state, 
                cov=self.proposal_cov
            )
            
            # Compute acceptance probability
            proposal_log_prob = self.log_posterior(proposal)
            log_alpha = proposal_log_prob - self.current_log_prob
            alpha = min(1.0, np.exp(log_alpha))
            
            # Accept or reject
            if np.random.rand() < alpha:
                self.current_state = proposal
                self.current_log_prob = proposal_log_prob
                self.n_accepted += 1
                
            self.n_total += 1
            
            # Store sample after burn-in and thinning
            if i >= burn_in and (i - burn_in) % thin == 0:
                samples.append(self.current_state.copy())
                
        acceptance_rate = self.n_accepted / self.n_total
        return np.array(samples), acceptance_rate
    
    def adaptive_proposal(self, target_acceptance: float = 0.44,
                         adaptation_rate: float = 0.01):
        """
        Adapt proposal covariance for optimal acceptance rate.
        
        Args:
            target_acceptance: Target acceptance rate
            adaptation_rate: Rate of covariance adaptation
        """
        current_acceptance = self.n_accepted / max(1, self.n_total)
        
        if current_acceptance > target_acceptance:
            # Increase proposal variance
            self.proposal_cov *= (1 + adaptation_rate)
        else:
            # Decrease proposal variance  
            self.proposal_cov *= (1 - adaptation_rate)


class GillespieSolver:
    """
    Gillespie algorithm for exact stochastic simulation of
    chemical reaction networks.
    """
    
    def __init__(self, reactions: List[Dict], initial_state: Dict[str, int]):
        """
        Initialize Gillespie solver.
        
        Args:
            reactions: List of reaction dictionaries with keys:
                      'reactants', 'products', 'rate_constant'
            initial_state: Initial molecule counts
        """
        self.reactions = reactions
        self.species = list(initial_state.keys())
        self.state = initial_state.copy()
        self.time = 0.0
        
        # Build stoichiometry matrix
        self.n_species = len(self.species)
        self.n_reactions = len(reactions)
        self.stoich_matrix = np.zeros((self.n_species, self.n_reactions))
        
        for j, reaction in enumerate(reactions):
            for species, coeff in reaction.get('reactants', {}).items():
                i = self.species.index(species)
                self.stoich_matrix[i, j] -= coeff
                
            for species, coeff in reaction.get('products', {}).items():
                i = self.species.index(species)
                self.stoich_matrix[i, j] += coeff
                
    def compute_propensities(self) -> np.ndarray:
        """
        Compute reaction propensities based on current state.
        
        Returns:
            Array of propensity values
        """
        propensities = np.zeros(self.n_reactions)
        
        for j, reaction in enumerate(self.reactions):
            rate = reaction['rate_constant']
            
            # Compute combinatorial factor
            for species, coeff in reaction.get('reactants', {}).items():
                count = self.state[species]
                for k in range(coeff):
                    rate *= max(0, count - k)
                    
            propensities[j] = rate
            
        return propensities
    
    def gillespie_step(self) -> Tuple[float, int]:
        """
        Perform one Gillespie step.
        
        Returns:
            tau: Time to next reaction
            reaction_index: Index of reaction that occurs
        """
        propensities = self.compute_propensities()
        total_propensity = np.sum(propensities)
        
        if total_propensity == 0:
            return np.inf, -1
            
        # Sample time to next reaction
        tau = np.random.exponential(1.0 / total_propensity)
        
        # Sample which reaction occurs
        cumulative_prop = np.cumsum(propensities)
        r2 = np.random.rand() * total_propensity
        reaction_index = np.searchsorted(cumulative_prop, r2)
        
        return tau, reaction_index
    
    def simulate(self, t_final: float, max_steps: int = 100000) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate the reaction network until final time.
        
        Args:
            t_final: Final simulation time
            max_steps: Maximum number of steps
            
        Returns:
            times: Array of time points
            populations: Array of population trajectories
        """
        times = [self.time]
        populations = [list(self.state.values())]
        
        step = 0
        while self.time < t_final and step < max_steps:
            tau, reaction_idx = self.gillespie_step()
            
            if reaction_idx == -1:  # No more reactions possible
                break
                
            # Update time and state
            self.time += tau
            
            # Apply reaction
            for i, species in enumerate(self.species):
                self.state[species] += int(self.stoich_matrix[i, reaction_idx])
                self.state[species] = max(0, self.state[species])  # Ensure non-negative
                
            times.append(self.time)
            populations.append(list(self.state.values()))
            step += 1
            
        return np.array(times), np.array(populations)


class OrnsteinUhlenbeckProcess:
    """
    Ornstein-Uhlenbeck process for modeling mean-reverting stochastic processes
    in biological systems.
    
    dX_t = -θ(X_t - μ)dt + σ dW_t
    """
    
    def __init__(self, theta: float, mu: float, sigma: float):
        """
        Initialize Ornstein-Uhlenbeck process.
        
        Args:
            theta: Mean reversion rate
            mu: Long-term mean
            sigma: Volatility parameter
        """
        self.theta = theta
        self.mu = mu
        self.sigma = sigma
        
    def simulate(self, x0: float, t_span: Tuple[float, float], 
                dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate OU process trajectory.
        
        Args:
            x0: Initial value
            t_span: (t_start, t_end) time interval
            dt: Time step
            
        Returns:
            times: Time points
            trajectory: Process values
        """
        t_start, t_end = t_span
        times = np.arange(t_start, t_end + dt, dt)
        n_steps = len(times)
        
        trajectory = np.zeros(n_steps)
        trajectory[0] = x0
        
        # Exact solution for discrete time steps
        exp_theta_dt = np.exp(-self.theta * dt)
        variance = self.sigma**2 * (1 - exp_theta_dt**2) / (2 * self.theta)
        
        for i in range(1, n_steps):
            mean = self.mu + (trajectory[i-1] - self.mu) * exp_theta_dt
            trajectory[i] = norm.rvs(loc=mean, scale=np.sqrt(variance))
            
        return times, trajectory
    
    def analytical_moments(self, x0: float, t: float) -> Tuple[float, float]:
        """
        Compute analytical mean and variance at time t.
        
        Args:
            x0: Initial value
            t: Time
            
        Returns:
            mean: Expected value at time t
            variance: Variance at time t
        """
        exp_theta_t = np.exp(-self.theta * t)
        
        mean = self.mu + (x0 - self.mu) * exp_theta_t
        variance = self.sigma**2 * (1 - exp_theta_t**2) / (2 * self.theta)
        
        return mean, variance


class WienerProcessGenerator:
    """
    Generator for Wiener processes (Brownian motion) and related
    stochastic processes.
    """
    
    @staticmethod
    def standard_brownian(t_span: Tuple[float, float], 
                         dt: float, n_paths: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate standard Brownian motion paths.
        
        Args:
            t_span: (t_start, t_end) time interval
            dt: Time step
            n_paths: Number of paths to generate
            
        Returns:
            times: Time points
            paths: Brownian motion paths (n_steps, n_paths)
        """
        t_start, t_end = t_span
        times = np.arange(t_start, t_end + dt, dt)
        n_steps = len(times)
        
        # Generate increments
        dW = norm.rvs(loc=0, scale=np.sqrt(dt), size=(n_steps-1, n_paths))
        
        # Cumulative sum to get Brownian motion
        paths = np.zeros((n_steps, n_paths))
        paths[1:] = np.cumsum(dW, axis=0)
        
        return times, paths
    
    @staticmethod
    def geometric_brownian(s0: float, mu: float, sigma: float,
                          t_span: Tuple[float, float], dt: float,
                          n_paths: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate geometric Brownian motion paths.
        
        dS_t = μS_t dt + σS_t dW_t
        
        Args:
            s0: Initial value
            mu: Drift parameter
            sigma: Volatility parameter
            t_span: Time interval
            dt: Time step
            n_paths: Number of paths
            
        Returns:
            times: Time points
            paths: GBM paths
        """
        times, W = WienerProcessGenerator.standard_brownian(t_span, dt, n_paths)
        
        # Exact solution: S_t = S_0 * exp((μ - σ²/2)t + σW_t)
        drift_term = (mu - 0.5 * sigma**2) * times[:, np.newaxis]
        diffusion_term = sigma * W
        
        paths = s0 * np.exp(drift_term + diffusion_term)
        
        return times, paths
    
    @staticmethod
    def fractional_brownian(hurst: float, t_span: Tuple[float, float],
                           dt: float, n_paths: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate fractional Brownian motion with Hurst parameter H.
        
        Args:
            hurst: Hurst parameter (0 < H < 1)
            t_span: Time interval
            dt: Time step
            n_paths: Number of paths
            
        Returns:
            times: Time points
            paths: Fractional Brownian motion paths
        """
        if not 0 < hurst < 1:
            raise ValueError("Hurst parameter must be between 0 and 1")
            
        t_start, t_end = t_span
        times = np.arange(t_start, t_end + dt, dt)
        n_steps = len(times)
        
        # Covariance matrix for fractional Brownian motion
        def fbm_covariance(s, t, H):
            return 0.5 * (s**(2*H) + t**(2*H) - np.abs(t-s)**(2*H))
        
        # Build covariance matrix
        cov_matrix = np.zeros((n_steps, n_steps))
        for i in range(n_steps):
            for j in range(n_steps):
                if times[i] > 0 and times[j] > 0:
                    cov_matrix[i, j] = fbm_covariance(times[i], times[j], hurst)
        
        # Generate correlated Gaussian variables
        paths = np.zeros((n_steps, n_paths))
        for path in range(n_paths):
            paths[:, path] = multivariate_normal.rvs(
                mean=np.zeros(n_steps), cov=cov_matrix)
        
        return times, paths


class StochasticDifferentialEquationSolver:
    """
    Numerical solver for stochastic differential equations using
    Euler-Maruyama and Milstein schemes.
    """
    
    def __init__(self, method: str = "euler_maruyama"):
        """
        Initialize SDE solver.
        
        Args:
            method: "euler_maruyama" or "milstein"
        """
        self.method = method
        
    def solve(self, drift: Callable, diffusion: Callable,
             y0: np.ndarray, t_span: Tuple[float, float], dt: float,
             diffusion_derivative: Optional[Callable] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve SDE: dY_t = f(t, Y_t)dt + g(t, Y_t)dW_t
        
        Args:
            drift: Drift function f(t, y)
            diffusion: Diffusion function g(t, y)  
            y0: Initial conditions
            t_span: Time interval
            dt: Time step
            diffusion_derivative: Derivative of diffusion function (for Milstein)
            
        Returns:
            times: Time points
            solution: Solution trajectory
        """
        t_start, t_end = t_span
        times = np.arange(t_start, t_end + dt, dt)
        n_steps = len(times)
        
        solution = np.zeros((n_steps, len(y0)))
        solution[0] = y0
        
        sqrt_dt = np.sqrt(dt)
        
        for i in range(1, n_steps):
            t = times[i-1]
            y = solution[i-1]
            
            # Generate Wiener increment
            dW = norm.rvs(loc=0, scale=sqrt_dt, size=len(y0))
            
            if self.method == "euler_maruyama":
                # Euler-Maruyama scheme
                solution[i] = y + drift(t, y) * dt + diffusion(t, y) * dW
                
            elif self.method == "milstein":
                # Milstein scheme (higher order)
                if diffusion_derivative is None:
                    raise ValueError("Milstein method requires diffusion derivative")
                    
                g = diffusion(t, y)
                g_prime = diffusion_derivative(t, y)
                
                solution[i] = (y + drift(t, y) * dt + g * dW + 
                             0.5 * g * g_prime * (dW**2 - dt))
            else:
                raise ValueError(f"Unknown method: {self.method}")
                
        return times, solution
    
    def strong_convergence_test(self, drift: Callable, diffusion: Callable,
                               y0: np.ndarray, t_span: Tuple[float, float],
                               reference_solution: Callable,
                               dt_values: List[float]) -> np.ndarray:
        """
        Test strong convergence of the numerical scheme.
        
        Args:
            drift: Drift function
            diffusion: Diffusion function
            y0: Initial conditions
            t_span: Time interval
            reference_solution: Analytical solution (if available)
            dt_values: List of time steps to test
            
        Returns:
            errors: Strong errors for each time step
        """
        errors = np.zeros(len(dt_values))
        
        # Use same random seed for fair comparison
        np.random.seed(42)
        
        for i, dt in enumerate(dt_values):
            times, numerical_sol = self.solve(drift, diffusion, y0, t_span, dt)
            
            # Compare with reference at final time
            analytical_val = reference_solution(times[-1])
            numerical_val = numerical_sol[-1]
            
            errors[i] = np.linalg.norm(numerical_val - analytical_val)
            
        return errors
