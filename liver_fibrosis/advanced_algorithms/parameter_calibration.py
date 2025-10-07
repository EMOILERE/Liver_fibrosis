"""
Parameter Calibration Module

Implements advanced parameter calibration and uncertainty quantification
techniques for biological system modeling including Bayesian inference,
MCMC sampling, and automatic differentiation.

Key Features:
- Bayesian parameter inference with prior specification
- Markov Chain Monte Carlo sampling algorithms
- Automatic differentiation for gradient-based optimization
- Comprehensive uncertainty quantification
- Model selection and comparison metrics
- Global sensitivity analysis with Sobol indices
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable, Any, Union
from dataclasses import dataclass
from enum import Enum
import warnings
from scipy import stats, optimize
from scipy.linalg import cholesky, solve_triangular
import time


class PriorType(Enum):
    """Types of prior distributions."""
    UNIFORM = "uniform"
    NORMAL = "normal"
    LOGNORMAL = "lognormal"
    GAMMA = "gamma"
    BETA = "beta"
    EXPONENTIAL = "exponential"


class MCMCAlgorithm(Enum):
    """MCMC sampling algorithms."""
    METROPOLIS_HASTINGS = "metropolis_hastings"
    HAMILTONIAN_MC = "hamiltonian_mc"
    NUTS = "nuts"
    GIBBS = "gibbs"
    ADAPTIVE_METROPOLIS = "adaptive_metropolis"


@dataclass
class PriorDistribution:
    """Container for prior distribution specification."""
    distribution_type: PriorType
    parameters: Dict[str, float]
    bounds: Optional[Tuple[float, float]] = None
    
    def sample(self, size: int = 1) -> np.ndarray:
        """Sample from the prior distribution."""
        if self.distribution_type == PriorType.UNIFORM:
            return np.random.uniform(
                self.parameters['low'], 
                self.parameters['high'], 
                size
            )
        elif self.distribution_type == PriorType.NORMAL:
            return np.random.normal(
                self.parameters['loc'], 
                self.parameters['scale'], 
                size
            )
        elif self.distribution_type == PriorType.LOGNORMAL:
            return np.random.lognormal(
                self.parameters['mean'], 
                self.parameters['sigma'], 
                size
            )
        elif self.distribution_type == PriorType.GAMMA:
            return np.random.gamma(
                self.parameters['shape'], 
                self.parameters['scale'], 
                size
            )
        elif self.distribution_type == PriorType.BETA:
            return np.random.beta(
                self.parameters['alpha'], 
                self.parameters['beta'], 
                size
            )
        elif self.distribution_type == PriorType.EXPONENTIAL:
            return np.random.exponential(
                self.parameters['scale'], 
                size
            )
        else:
            raise ValueError(f"Unknown prior type: {self.distribution_type}")
    
    def log_pdf(self, x: np.ndarray) -> np.ndarray:
        """Calculate log probability density."""
        if self.distribution_type == PriorType.UNIFORM:
            low, high = self.parameters['low'], self.parameters['high']
            return np.where(
                (x >= low) & (x <= high),
                -np.log(high - low),
                -np.inf
            )
        elif self.distribution_type == PriorType.NORMAL:
            return stats.norm.logpdf(
                x, 
                self.parameters['loc'], 
                self.parameters['scale']
            )
        elif self.distribution_type == PriorType.LOGNORMAL:
            return stats.lognorm.logpdf(
                x, 
                self.parameters['sigma'], 
                scale=np.exp(self.parameters['mean'])
            )
        elif self.distribution_type == PriorType.GAMMA:
            return stats.gamma.logpdf(
                x, 
                self.parameters['shape'], 
                scale=self.parameters['scale']
            )
        elif self.distribution_type == PriorType.BETA:
            return stats.beta.logpdf(
                x, 
                self.parameters['alpha'], 
                self.parameters['beta']
            )
        elif self.distribution_type == PriorType.EXPONENTIAL:
            return stats.expon.logpdf(x, scale=self.parameters['scale'])
        else:
            raise ValueError(f"Unknown prior type: {self.distribution_type}")


@dataclass
class CalibrationResult:
    """Container for calibration results."""
    method: str
    parameters: Dict[str, Any]
    posterior_samples: Optional[np.ndarray] = None
    parameter_names: Optional[List[str]] = None
    log_likelihood_trace: Optional[np.ndarray] = None
    acceptance_rate: Optional[float] = None
    convergence_diagnostics: Optional[Dict[str, float]] = None
    uncertainty_metrics: Optional[Dict[str, Any]] = None
    computation_time: float = 0.0


class BayesianParameterInference:
    """
    Bayesian parameter inference framework with multiple sampling algorithms.
    """
    
    def __init__(self, model_function: Callable, 
                 parameter_names: List[str],
                 priors: Dict[str, PriorDistribution]):
        """
        Initialize Bayesian parameter inference.
        
        Args:
            model_function: Function that takes parameters and returns model predictions
            parameter_names: Names of parameters to calibrate
            priors: Prior distributions for each parameter
        """
        self.model_function = model_function
        self.parameter_names = parameter_names
        self.priors = priors
        
        # Validate that all parameters have priors
        for param_name in parameter_names:
            if param_name not in priors:
                raise ValueError(f"No prior specified for parameter: {param_name}")
    
    def log_prior(self, parameters: np.ndarray) -> float:
        """Calculate log prior probability."""
        log_prob = 0.0
        
        for i, param_name in enumerate(self.parameter_names):
            prior = self.priors[param_name]
            log_prob += prior.log_pdf(parameters[i])
        
        return log_prob
    
    def log_likelihood(self, parameters: np.ndarray, 
                      observed_data: np.ndarray,
                      observation_noise: float = 0.1) -> float:
        """
        Calculate log likelihood of parameters given observed data.
        
        Args:
            parameters: Parameter values
            observed_data: Observed data for comparison
            observation_noise: Observation noise standard deviation
            
        Returns:
            Log likelihood value
        """
        try:
            # Get model predictions
            predictions = self.model_function(parameters)
            
            # Ensure predictions and data have same shape
            if predictions.shape != observed_data.shape:
                return -np.inf
            
            # Calculate log likelihood assuming Gaussian noise
            residuals = predictions - observed_data
            log_likelihood = -0.5 * np.sum(residuals**2) / (observation_noise**2)
            log_likelihood -= 0.5 * len(observed_data) * np.log(2 * np.pi * observation_noise**2)
            
            return log_likelihood
            
        except Exception as e:
            warnings.warn(f"Model evaluation failed: {str(e)}")
            return -np.inf
    
    def log_posterior(self, parameters: np.ndarray, 
                     observed_data: np.ndarray,
                     observation_noise: float = 0.1) -> float:
        """Calculate log posterior probability."""
        log_prior_val = self.log_prior(parameters)
        
        if np.isfinite(log_prior_val):
            log_likelihood_val = self.log_likelihood(parameters, observed_data, observation_noise)
            return log_prior_val + log_likelihood_val
        else:
            return -np.inf
    
    def maximum_a_posteriori(self, observed_data: np.ndarray,
                           observation_noise: float = 0.1,
                           initial_guess: Optional[np.ndarray] = None) -> CalibrationResult:
        """
        Find maximum a posteriori (MAP) estimate.
        
        Args:
            observed_data: Observed data for calibration
            observation_noise: Observation noise level
            initial_guess: Initial parameter guess
            
        Returns:
            MAP calibration result
        """
        start_time = time.time()
        
        # Initial guess from prior means if not provided
        if initial_guess is None:
            initial_guess = np.array([
                self.priors[name].sample(1)[0] for name in self.parameter_names
            ])
        
        # Define objective function (negative log posterior)
        def objective(params):
            return -self.log_posterior(params, observed_data, observation_noise)
        
        # Optimize
        result = optimize.minimize(
            objective,
            initial_guess,
            method='L-BFGS-B',
            options={'maxiter': 1000}
        )
        
        computation_time = time.time() - start_time
        
        return CalibrationResult(
            method="MAP",
            parameters={
                'map_estimate': result.x,
                'optimization_result': result
            },
            parameter_names=self.parameter_names,
            computation_time=computation_time
        )


class MarkovChainMonteCarlo:
    """
    Markov Chain Monte Carlo sampling algorithms for Bayesian inference.
    """
    
    def __init__(self, inference_framework: BayesianParameterInference):
        """Initialize MCMC sampler."""
        self.inference = inference_framework
        
    def metropolis_hastings(self, observed_data: np.ndarray,
                          n_samples: int = 10000,
                          burn_in: int = 1000,
                          thin: int = 1,
                          proposal_cov: Optional[np.ndarray] = None,
                          observation_noise: float = 0.1) -> CalibrationResult:
        """
        Metropolis-Hastings MCMC sampling.
        
        Args:
            observed_data: Observed data for calibration
            n_samples: Number of samples to generate
            burn_in: Number of burn-in samples
            thin: Thinning interval
            proposal_cov: Proposal covariance matrix
            observation_noise: Observation noise level
            
        Returns:
            MCMC calibration result
        """
        start_time = time.time()
        
        n_params = len(self.inference.parameter_names)
        
        # Initialize proposal covariance
        if proposal_cov is None:
            proposal_cov = 0.1 * np.eye(n_params)
        
        # Initialize chain
        current_params = np.array([
            self.inference.priors[name].sample(1)[0] 
            for name in self.inference.parameter_names
        ])
        
        current_log_posterior = self.inference.log_posterior(
            current_params, observed_data, observation_noise
        )
        
        # Storage
        total_iterations = burn_in + n_samples * thin
        samples = []
        log_posterior_trace = []
        n_accepted = 0
        
        # MCMC loop
        for i in range(total_iterations):
            # Propose new parameters
            proposal = np.random.multivariate_normal(current_params, proposal_cov)
            
            # Calculate acceptance probability
            proposal_log_posterior = self.inference.log_posterior(
                proposal, observed_data, observation_noise
            )
            
            log_alpha = proposal_log_posterior - current_log_posterior
            alpha = min(1.0, np.exp(log_alpha))
            
            # Accept or reject
            if np.random.rand() < alpha:
                current_params = proposal
                current_log_posterior = proposal_log_posterior
                n_accepted += 1
            
            # Store sample after burn-in and thinning
            if i >= burn_in and (i - burn_in) % thin == 0:
                samples.append(current_params.copy())
                log_posterior_trace.append(current_log_posterior)
        
        samples = np.array(samples)
        log_posterior_trace = np.array(log_posterior_trace)
        acceptance_rate = n_accepted / total_iterations
        
        # Convergence diagnostics
        convergence_diagnostics = self._calculate_convergence_diagnostics(samples)
        
        # Uncertainty quantification
        uncertainty_metrics = self._calculate_uncertainty_metrics(samples)
        
        computation_time = time.time() - start_time
        
        return CalibrationResult(
            method="Metropolis-Hastings",
            parameters={
                'n_samples': n_samples,
                'burn_in': burn_in,
                'thin': thin,
                'proposal_cov': proposal_cov
            },
            posterior_samples=samples,
            parameter_names=self.inference.parameter_names,
            log_likelihood_trace=log_posterior_trace,
            acceptance_rate=acceptance_rate,
            convergence_diagnostics=convergence_diagnostics,
            uncertainty_metrics=uncertainty_metrics,
            computation_time=computation_time
        )
    
    def adaptive_metropolis(self, observed_data: np.ndarray,
                          n_samples: int = 10000,
                          burn_in: int = 1000,
                          adaptation_interval: int = 100,
                          target_acceptance: float = 0.44) -> CalibrationResult:
        """
        Adaptive Metropolis algorithm with covariance adaptation.
        
        Args:
            observed_data: Observed data for calibration
            n_samples: Number of samples to generate
            burn_in: Number of burn-in samples
            adaptation_interval: Interval for covariance adaptation
            target_acceptance: Target acceptance rate
            
        Returns:
            Adaptive MCMC calibration result
        """
        start_time = time.time()
        
        n_params = len(self.inference.parameter_names)
        
        # Initialize
        current_params = np.array([
            self.inference.priors[name].sample(1)[0] 
            for name in self.inference.parameter_names
        ])
        
        proposal_cov = 0.1 * np.eye(n_params)
        current_log_posterior = self.inference.log_posterior(
            current_params, observed_data
        )
        
        # Storage
        total_iterations = burn_in + n_samples
        samples = []
        log_posterior_trace = []
        n_accepted = 0
        adaptation_history = []
        
        # Running statistics for adaptation
        sample_mean = current_params.copy()
        sample_cov = proposal_cov.copy()
        
        for i in range(total_iterations):
            # Propose new parameters
            proposal = np.random.multivariate_normal(current_params, proposal_cov)
            
            # Calculate acceptance probability
            proposal_log_posterior = self.inference.log_posterior(
                proposal, observed_data
            )
            
            log_alpha = proposal_log_posterior - current_log_posterior
            alpha = min(1.0, np.exp(log_alpha))
            
            # Accept or reject
            if np.random.rand() < alpha:
                current_params = proposal
                current_log_posterior = proposal_log_posterior
                n_accepted += 1
            
            # Adapt proposal covariance
            if i > 0 and i % adaptation_interval == 0:
                current_acceptance = n_accepted / (i + 1)
                
                # Update running statistics
                if i > adaptation_interval:
                    # Update sample mean and covariance
                    n_samples_so_far = len(samples)
                    if n_samples_so_far > 1:
                        sample_mean = np.mean(samples, axis=0)
                        sample_cov = np.cov(samples, rowvar=False)
                        
                        # Regularize covariance matrix
                        sample_cov += 1e-6 * np.eye(n_params)
                        
                        # Adaptive scaling
                        if current_acceptance > target_acceptance:
                            scale_factor = 1.1
                        else:
                            scale_factor = 0.9
                        
                        proposal_cov = scale_factor * sample_cov
                
                adaptation_history.append({
                    'iteration': i,
                    'acceptance_rate': current_acceptance,
                    'proposal_cov_trace': np.trace(proposal_cov)
                })
            
            # Store sample after burn-in
            if i >= burn_in:
                samples.append(current_params.copy())
                log_posterior_trace.append(current_log_posterior)
        
        samples = np.array(samples)
        log_posterior_trace = np.array(log_posterior_trace)
        acceptance_rate = n_accepted / total_iterations
        
        # Convergence diagnostics
        convergence_diagnostics = self._calculate_convergence_diagnostics(samples)
        
        # Uncertainty quantification
        uncertainty_metrics = self._calculate_uncertainty_metrics(samples)
        
        computation_time = time.time() - start_time
        
        return CalibrationResult(
            method="Adaptive Metropolis",
            parameters={
                'n_samples': n_samples,
                'burn_in': burn_in,
                'adaptation_interval': adaptation_interval,
                'target_acceptance': target_acceptance,
                'adaptation_history': adaptation_history
            },
            posterior_samples=samples,
            parameter_names=self.inference.parameter_names,
            log_likelihood_trace=log_posterior_trace,
            acceptance_rate=acceptance_rate,
            convergence_diagnostics=convergence_diagnostics,
            uncertainty_metrics=uncertainty_metrics,
            computation_time=computation_time
        )
    
    def _calculate_convergence_diagnostics(self, samples: np.ndarray) -> Dict[str, float]:
        """Calculate convergence diagnostics for MCMC chains."""
        n_samples, n_params = samples.shape
        
        diagnostics = {}
        
        # Effective sample size (simplified)
        for i, param_name in enumerate(self.inference.parameter_names):
            param_samples = samples[:, i]
            
            # Autocorrelation function
            autocorr = np.correlate(param_samples, param_samples, mode='full')
            autocorr = autocorr[autocorr.size // 2:]
            autocorr = autocorr / autocorr[0]
            
            # Find first negative autocorrelation
            first_negative = np.where(autocorr < 0)[0]
            if len(first_negative) > 0:
                tau_int = 1 + 2 * np.sum(autocorr[1:first_negative[0]])
            else:
                tau_int = n_samples  # Conservative estimate
            
            eff_sample_size = n_samples / (2 * tau_int + 1)
            diagnostics[f'{param_name}_eff_sample_size'] = eff_sample_size
        
        # Gelman-Rubin statistic (simplified, single chain)
        # Split chain into two halves
        mid_point = n_samples // 2
        first_half = samples[:mid_point]
        second_half = samples[mid_point:]
        
        for i, param_name in enumerate(self.inference.parameter_names):
            mean1 = np.mean(first_half[:, i])
            mean2 = np.mean(second_half[:, i])
            var1 = np.var(first_half[:, i])
            var2 = np.var(second_half[:, i])
            
            # Between-chain variance
            B = mid_point * (mean1 - mean2)**2 / 2
            
            # Within-chain variance
            W = (var1 + var2) / 2
            
            # Potential scale reduction factor
            if W > 0:
                R_hat = np.sqrt((mid_point - 1) / mid_point + B / (mid_point * W))
            else:
                R_hat = 1.0
            
            diagnostics[f'{param_name}_R_hat'] = R_hat
        
        return diagnostics
    
    def _calculate_uncertainty_metrics(self, samples: np.ndarray) -> Dict[str, Any]:
        """Calculate uncertainty quantification metrics."""
        metrics = {}
        
        for i, param_name in enumerate(self.inference.parameter_names):
            param_samples = samples[:, i]
            
            metrics[param_name] = {
                'mean': np.mean(param_samples),
                'std': np.std(param_samples),
                'median': np.median(param_samples),
                'q025': np.percentile(param_samples, 2.5),
                'q975': np.percentile(param_samples, 97.5),
                'mode': self._calculate_mode(param_samples)
            }
        
        # Correlation matrix
        correlation_matrix = np.corrcoef(samples, rowvar=False)
        metrics['correlation_matrix'] = correlation_matrix
        
        return metrics
    
    def _calculate_mode(self, samples: np.ndarray, n_bins: int = 50) -> float:
        """Calculate mode using histogram."""
        hist, bin_edges = np.histogram(samples, bins=n_bins)
        max_bin_idx = np.argmax(hist)
        mode = (bin_edges[max_bin_idx] + bin_edges[max_bin_idx + 1]) / 2
        return mode


class AutomaticDifferentiation:
    """
    Automatic differentiation for gradient-based parameter optimization.
    """
    
    def __init__(self, function: Callable):
        """
        Initialize automatic differentiation.
        
        Args:
            function: Function to differentiate
        """
        self.function = function
        
    def finite_difference_gradient(self, x: np.ndarray, 
                                 h: float = 1e-8) -> np.ndarray:
        """
        Calculate gradient using finite differences.
        
        Args:
            x: Point at which to calculate gradient
            h: Step size for finite differences
            
        Returns:
            Gradient vector
        """
        n = len(x)
        gradient = np.zeros(n)
        
        f_x = self.function(x)
        
        for i in range(n):
            x_plus = x.copy()
            x_plus[i] += h
            
            f_x_plus = self.function(x_plus)
            gradient[i] = (f_x_plus - f_x) / h
        
        return gradient
    
    def central_difference_gradient(self, x: np.ndarray,
                                  h: float = 1e-6) -> np.ndarray:
        """
        Calculate gradient using central differences (more accurate).
        
        Args:
            x: Point at which to calculate gradient
            h: Step size for finite differences
            
        Returns:
            Gradient vector
        """
        n = len(x)
        gradient = np.zeros(n)
        
        for i in range(n):
            x_plus = x.copy()
            x_minus = x.copy()
            x_plus[i] += h
            x_minus[i] -= h
            
            f_x_plus = self.function(x_plus)
            f_x_minus = self.function(x_minus)
            
            gradient[i] = (f_x_plus - f_x_minus) / (2 * h)
        
        return gradient
    
    def hessian_finite_difference(self, x: np.ndarray,
                                h: float = 1e-6) -> np.ndarray:
        """
        Calculate Hessian matrix using finite differences.
        
        Args:
            x: Point at which to calculate Hessian
            h: Step size for finite differences
            
        Returns:
            Hessian matrix
        """
        n = len(x)
        hessian = np.zeros((n, n))
        
        f_x = self.function(x)
        
        # Diagonal elements
        for i in range(n):
            x_plus = x.copy()
            x_minus = x.copy()
            x_plus[i] += h
            x_minus[i] -= h
            
            f_x_plus = self.function(x_plus)
            f_x_minus = self.function(x_minus)
            
            hessian[i, i] = (f_x_plus - 2 * f_x + f_x_minus) / (h**2)
        
        # Off-diagonal elements
        for i in range(n):
            for j in range(i + 1, n):
                x_pp = x.copy()
                x_pm = x.copy()
                x_mp = x.copy()
                x_mm = x.copy()
                
                x_pp[i] += h
                x_pp[j] += h
                
                x_pm[i] += h
                x_pm[j] -= h
                
                x_mp[i] -= h
                x_mp[j] += h
                
                x_mm[i] -= h
                x_mm[j] -= h
                
                f_pp = self.function(x_pp)
                f_pm = self.function(x_pm)
                f_mp = self.function(x_mp)
                f_mm = self.function(x_mm)
                
                hessian[i, j] = (f_pp - f_pm - f_mp + f_mm) / (4 * h**2)
                hessian[j, i] = hessian[i, j]  # Symmetry
        
        return hessian


class UncertaintyQuantification:
    """
    Comprehensive uncertainty quantification for model parameters and predictions.
    """
    
    def __init__(self):
        """Initialize uncertainty quantification."""
        pass
    
    def propagate_uncertainty(self, model_function: Callable,
                            parameter_samples: np.ndarray,
                            n_predictions: int = 1000) -> Dict[str, Any]:
        """
        Propagate parameter uncertainty to model predictions.
        
        Args:
            model_function: Model function
            parameter_samples: Samples from parameter posterior
            n_predictions: Number of prediction samples
            
        Returns:
            Uncertainty propagation results
        """
        # Select random subset of parameter samples
        n_param_samples = len(parameter_samples)
        selected_indices = np.random.choice(
            n_param_samples, 
            min(n_predictions, n_param_samples), 
            replace=False
        )
        
        selected_samples = parameter_samples[selected_indices]
        
        # Generate predictions
        predictions = []
        for params in selected_samples:
            try:
                pred = model_function(params)
                predictions.append(pred)
            except Exception as e:
                warnings.warn(f"Prediction failed for parameters {params}: {str(e)}")
        
        if not predictions:
            return {'error': 'No successful predictions'}
        
        predictions = np.array(predictions)
        
        # Calculate prediction statistics
        results = {
            'predictions': predictions,
            'mean_prediction': np.mean(predictions, axis=0),
            'std_prediction': np.std(predictions, axis=0),
            'median_prediction': np.median(predictions, axis=0),
            'q025_prediction': np.percentile(predictions, 2.5, axis=0),
            'q975_prediction': np.percentile(predictions, 97.5, axis=0),
            'prediction_intervals': {
                '50%': (np.percentile(predictions, 25, axis=0),
                       np.percentile(predictions, 75, axis=0)),
                '95%': (np.percentile(predictions, 2.5, axis=0),
                       np.percentile(predictions, 97.5, axis=0))
            }
        }
        
        return results
    
    def sensitivity_analysis_morris(self, model_function: Callable,
                                  parameter_bounds: List[Tuple[float, float]],
                                  n_trajectories: int = 100,
                                  n_levels: int = 4) -> Dict[str, Any]:
        """
        Morris sensitivity analysis (elementary effects method).
        
        Args:
            model_function: Model function to analyze
            parameter_bounds: Parameter bounds for sampling
            n_trajectories: Number of trajectories
            n_levels: Number of levels for parameter grid
            
        Returns:
            Morris sensitivity analysis results
        """
        n_params = len(parameter_bounds)
        
        # Generate Morris sample
        elementary_effects = {i: [] for i in range(n_params)}
        
        for trajectory in range(n_trajectories):
            # Generate base point
            base_point = np.array([
                np.random.uniform(bounds[0], bounds[1]) 
                for bounds in parameter_bounds
            ])
            
            # Evaluate at base point
            base_output = model_function(base_point)
            
            # Generate trajectory
            current_point = base_point.copy()
            
            # Random permutation of parameters
            param_order = np.random.permutation(n_params)
            
            for param_idx in param_order:
                # Calculate step size
                bounds = parameter_bounds[param_idx]
                step_size = (bounds[1] - bounds[0]) / (n_levels - 1)
                
                # Move in parameter space
                new_point = current_point.copy()
                if np.random.rand() > 0.5:
                    new_point[param_idx] = min(bounds[1], current_point[param_idx] + step_size)
                else:
                    new_point[param_idx] = max(bounds[0], current_point[param_idx] - step_size)
                
                # Evaluate at new point
                try:
                    new_output = model_function(new_point)
                    
                    # Calculate elementary effect
                    if isinstance(base_output, np.ndarray):
                        effect = np.mean(np.abs(new_output - base_output))
                    else:
                        effect = abs(new_output - base_output)
                    
                    elementary_effects[param_idx].append(effect)
                    
                    # Update for next iteration
                    current_point = new_point
                    base_output = new_output
                    
                except Exception as e:
                    warnings.warn(f"Morris analysis failed at point {new_point}: {str(e)}")
        
        # Calculate Morris statistics
        results = {}
        for param_idx in range(n_params):
            effects = elementary_effects[param_idx]
            if effects:
                results[f'parameter_{param_idx}'] = {
                    'mean_elementary_effect': np.mean(effects),
                    'std_elementary_effect': np.std(effects),
                    'mean_absolute_effect': np.mean(np.abs(effects))
                }
        
        return results


class ModelSelection:
    """
    Model selection and comparison framework.
    """
    
    def __init__(self):
        """Initialize model selection."""
        pass
    
    def information_criteria(self, log_likelihood: float,
                           n_parameters: int,
                           n_observations: int) -> Dict[str, float]:
        """
        Calculate information criteria for model selection.
        
        Args:
            log_likelihood: Maximum log likelihood
            n_parameters: Number of model parameters
            n_observations: Number of observations
            
        Returns:
            Information criteria values
        """
        # Akaike Information Criterion
        aic = 2 * n_parameters - 2 * log_likelihood
        
        # Corrected AIC for small samples
        if n_observations / n_parameters < 40:
            aicc = aic + (2 * n_parameters * (n_parameters + 1)) / (n_observations - n_parameters - 1)
        else:
            aicc = aic
        
        # Bayesian Information Criterion
        bic = n_parameters * np.log(n_observations) - 2 * log_likelihood
        
        # Deviance Information Criterion (simplified)
        dic = -2 * log_likelihood + 2 * n_parameters
        
        return {
            'aic': aic,
            'aicc': aicc,
            'bic': bic,
            'dic': dic,
            'log_likelihood': log_likelihood
        }
    
    def cross_validation_score(self, model_function: Callable,
                             parameters: np.ndarray,
                             data: np.ndarray,
                             k_folds: int = 5) -> Dict[str, float]:
        """
        Calculate cross-validation score for model assessment.
        
        Args:
            model_function: Model function
            parameters: Model parameters
            data: Observed data
            k_folds: Number of cross-validation folds
            
        Returns:
            Cross-validation results
        """
        n_observations = len(data)
        fold_size = n_observations // k_folds
        
        cv_scores = []
        
        for fold in range(k_folds):
            # Split data
            start_idx = fold * fold_size
            end_idx = start_idx + fold_size if fold < k_folds - 1 else n_observations
            
            test_indices = list(range(start_idx, end_idx))
            train_indices = [i for i in range(n_observations) if i not in test_indices]
            
            train_data = data[train_indices]
            test_data = data[test_indices]
            
            try:
                # Generate predictions for test data
                predictions = model_function(parameters)
                
                if len(predictions) >= len(test_data):
                    test_predictions = predictions[test_indices]
                    
                    # Calculate mean squared error
                    mse = np.mean((test_predictions - test_data)**2)
                    cv_scores.append(mse)
                    
            except Exception as e:
                warnings.warn(f"Cross-validation failed for fold {fold}: {str(e)}")
        
        if cv_scores:
            return {
                'cv_score': np.mean(cv_scores),
                'cv_std': np.std(cv_scores),
                'individual_scores': cv_scores
            }
        else:
            return {'error': 'Cross-validation failed for all folds'}
    
    def model_comparison(self, models: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        Compare multiple models using various criteria.
        
        Args:
            models: Dictionary of models with their statistics
            
        Returns:
            Model comparison results
        """
        comparison_results = {
            'model_rankings': {},
            'criteria_comparison': {},
            'best_model_by_criterion': {}
        }
        
        criteria = ['aic', 'bic', 'cv_score']
        
        for criterion in criteria:
            if all(criterion in model_stats for model_stats in models.values()):
                # Extract criterion values
                criterion_values = {
                    model_name: model_stats[criterion] 
                    for model_name, model_stats in models.items()
                }
                
                # Rank models (lower is better for AIC, BIC, CV score)
                sorted_models = sorted(criterion_values.items(), key=lambda x: x[1])
                
                comparison_results['criteria_comparison'][criterion] = criterion_values
                comparison_results['model_rankings'][criterion] = [
                    model_name for model_name, _ in sorted_models
                ]
                comparison_results['best_model_by_criterion'][criterion] = sorted_models[0][0]
        
        return comparison_results


class SensitivityAnalyzer:
    """
    Global sensitivity analysis using Sobol indices and other methods.
    """
    
    def __init__(self):
        """Initialize sensitivity analyzer."""
        pass
    
    def sobol_indices(self, model_function: Callable,
                     parameter_bounds: List[Tuple[float, float]],
                     n_samples: int = 1000) -> Dict[str, Any]:
        """
        Calculate Sobol sensitivity indices.
        
        Args:
            model_function: Model function to analyze
            parameter_bounds: Parameter bounds
            n_samples: Number of samples for estimation
            
        Returns:
            Sobol sensitivity indices
        """
        n_params = len(parameter_bounds)
        
        # Generate sample matrices using Sobol sequences (simplified)
        # In practice, would use proper Sobol sequence generator
        
        # Matrix A
        A = np.random.uniform(0, 1, (n_samples, n_params))
        for i in range(n_params):
            bounds = parameter_bounds[i]
            A[:, i] = bounds[0] + A[:, i] * (bounds[1] - bounds[0])
        
        # Matrix B
        B = np.random.uniform(0, 1, (n_samples, n_params))
        for i in range(n_params):
            bounds = parameter_bounds[i]
            B[:, i] = bounds[0] + B[:, i] * (bounds[1] - bounds[0])
        
        # Evaluate model at A and B
        f_A = np.array([model_function(A[i]) for i in range(n_samples)])
        f_B = np.array([model_function(B[i]) for i in range(n_samples)])
        
        # Calculate total variance
        f_total = np.concatenate([f_A, f_B])
        V = np.var(f_total)
        
        # First-order indices
        S1 = np.zeros(n_params)
        
        for i in range(n_params):
            # Create C_i matrix (A with i-th column from B)
            C_i = A.copy()
            C_i[:, i] = B[:, i]
            
            f_C_i = np.array([model_function(C_i[j]) for j in range(n_samples)])
            
            # First-order index
            if V > 0:
                S1[i] = np.mean(f_B * (f_C_i - f_A)) / V
            else:
                S1[i] = 0.0
        
        # Total-order indices
        ST = np.zeros(n_params)
        
        for i in range(n_params):
            # Create C_i matrix (B with i-th column from A)
            C_i = B.copy()
            C_i[:, i] = A[:, i]
            
            f_C_i = np.array([model_function(C_i[j]) for j in range(n_samples)])
            
            # Total-order index
            if V > 0:
                ST[i] = 1 - np.mean(f_A * (f_C_i - f_B)) / V
            else:
                ST[i] = 0.0
        
        return {
            'first_order_indices': S1,
            'total_order_indices': ST,
            'total_variance': V,
            'parameter_bounds': parameter_bounds
        }
