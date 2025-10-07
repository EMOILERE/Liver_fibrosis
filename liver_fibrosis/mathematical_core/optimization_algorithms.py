"""
Optimization Algorithms Module

Implements advanced optimization methods for parameter estimation,
experimental design, and multi-objective optimization in biological systems.

Key Features:
- Bayesian optimization with Gaussian processes
- NSGA-II multi-objective optimization
- Response surface methodology
- Pareto front calculation
- Sobol sensitivity analysis
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.stats import norm, qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, WhiteKernel
from typing import Callable, Tuple, List, Dict, Optional, Union
import warnings


class BayesianOptimizer:
    """
    Bayesian optimization using Gaussian processes for expensive
    black-box function optimization.
    """
    
    def __init__(self, bounds: List[Tuple[float, float]], 
                 kernel: Optional[object] = None,
                 acquisition_function: str = "expected_improvement",
                 xi: float = 0.01, kappa: float = 2.576):
        """
        Initialize Bayesian optimizer.
        
        Args:
            bounds: List of (min, max) bounds for each parameter
            kernel: GP kernel (default: Matern)
            acquisition_function: "expected_improvement", "upper_confidence_bound", or "probability_improvement"
            xi: Exploration parameter for EI
            kappa: Exploration parameter for UCB
        """
        self.bounds = np.array(bounds)
        self.dim = len(bounds)
        self.xi = xi
        self.kappa = kappa
        self.acquisition_func = acquisition_function
        
        # Default kernel: Matern with automatic relevance determination
        if kernel is None:
            kernel = Matern(length_scale=np.ones(self.dim), nu=2.5) + WhiteKernel()
        
        self.gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=1e-6,
            normalize_y=True,
            n_restarts_optimizer=10
        )
        
        self.X_observed = []
        self.y_observed = []
        
    def _normalize_bounds(self, X: np.ndarray) -> np.ndarray:
        """Normalize parameters to [0, 1] range."""
        return (X - self.bounds[:, 0]) / (self.bounds[:, 1] - self.bounds[:, 0])
    
    def _denormalize_bounds(self, X_norm: np.ndarray) -> np.ndarray:
        """Denormalize parameters from [0, 1] to original range."""
        return X_norm * (self.bounds[:, 1] - self.bounds[:, 0]) + self.bounds[:, 0]
    
    def _expected_improvement(self, X: np.ndarray) -> np.ndarray:
        """
        Compute Expected Improvement acquisition function.
        
        Args:
            X: Candidate points (n_points, n_dims)
            
        Returns:
            EI values for each point
        """
        X_norm = self._normalize_bounds(X)
        
        if len(self.X_observed) == 0:
            return np.ones(X.shape[0])
            
        mu, sigma = self.gp.predict(X_norm, return_std=True)
        sigma = sigma.reshape(-1, 1)
        
        # Current best observation
        f_best = np.max(self.y_observed)
        
        # Expected improvement
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            improvement = mu - f_best - self.xi
            Z = improvement / sigma
            ei = improvement * norm.cdf(Z) + sigma * norm.pdf(Z)
            ei[sigma == 0.0] = 0.0
            
        return ei.flatten()
    
    def _upper_confidence_bound(self, X: np.ndarray) -> np.ndarray:
        """
        Compute Upper Confidence Bound acquisition function.
        
        Args:
            X: Candidate points
            
        Returns:
            UCB values
        """
        X_norm = self._normalize_bounds(X)
        
        if len(self.X_observed) == 0:
            return np.ones(X.shape[0])
            
        mu, sigma = self.gp.predict(X_norm, return_std=True)
        return mu + self.kappa * sigma
    
    def _probability_improvement(self, X: np.ndarray) -> np.ndarray:
        """
        Compute Probability of Improvement acquisition function.
        
        Args:
            X: Candidate points
            
        Returns:
            PI values
        """
        X_norm = self._normalize_bounds(X)
        
        if len(self.X_observed) == 0:
            return np.ones(X.shape[0])
            
        mu, sigma = self.gp.predict(X_norm, return_std=True)
        f_best = np.max(self.y_observed)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            Z = (mu - f_best - self.xi) / sigma
            pi = norm.cdf(Z)
            
        return pi
    
    def acquisition(self, X: np.ndarray) -> np.ndarray:
        """
        Compute acquisition function values.
        
        Args:
            X: Candidate points
            
        Returns:
            Acquisition values
        """
        if self.acquisition_func == "expected_improvement":
            return self._expected_improvement(X)
        elif self.acquisition_func == "upper_confidence_bound":
            return self._upper_confidence_bound(X)
        elif self.acquisition_func == "probability_improvement":
            return self._probability_improvement(X)
        else:
            raise ValueError(f"Unknown acquisition function: {self.acquisition_func}")
    
    def suggest_next_point(self) -> np.ndarray:
        """
        Suggest next point to evaluate using acquisition function optimization.
        
        Returns:
            Next point to evaluate
        """
        def objective(x):
            return -self.acquisition(x.reshape(1, -1))[0]
        
        # Multi-start optimization
        best_x = None
        best_acq = -np.inf
        
        for _ in range(10):
            # Random starting point
            x0 = np.random.uniform(self.bounds[:, 0], self.bounds[:, 1])
            
            result = minimize(
                objective, x0,
                bounds=self.bounds,
                method='L-BFGS-B'
            )
            
            if -result.fun > best_acq:
                best_acq = -result.fun
                best_x = result.x
                
        return best_x
    
    def update(self, X_new: np.ndarray, y_new: float):
        """
        Update GP with new observation.
        
        Args:
            X_new: New parameter vector
            y_new: New objective value
        """
        self.X_observed.append(X_new)
        self.y_observed.append(y_new)
        
        # Fit GP to all observations
        X_norm = self._normalize_bounds(np.array(self.X_observed))
        self.gp.fit(X_norm, np.array(self.y_observed))
    
    def optimize(self, objective_func: Callable, n_iterations: int = 50,
                n_initial: int = 5) -> Tuple[np.ndarray, float]:
        """
        Run Bayesian optimization.
        
        Args:
            objective_func: Function to optimize
            n_iterations: Number of optimization iterations
            n_initial: Number of initial random samples
            
        Returns:
            best_x: Best parameter vector found
            best_y: Best objective value found
        """
        # Initial random sampling
        for _ in range(n_initial):
            x = np.random.uniform(self.bounds[:, 0], self.bounds[:, 1])
            y = objective_func(x)
            self.update(x, y)
        
        # Bayesian optimization loop
        for iteration in range(n_iterations):
            # Suggest next point
            x_next = self.suggest_next_point()
            
            # Evaluate objective
            y_next = objective_func(x_next)
            
            # Update GP
            self.update(x_next, y_next)
            
        # Return best point found
        best_idx = np.argmax(self.y_observed)
        return self.X_observed[best_idx], self.y_observed[best_idx]


class NSGA2Optimizer:
    """
    NSGA-II (Non-dominated Sorting Genetic Algorithm II) for
    multi-objective optimization problems.
    """
    
    def __init__(self, bounds: List[Tuple[float, float]], 
                 population_size: int = 100,
                 crossover_prob: float = 0.9,
                 mutation_prob: float = 0.1,
                 eta_c: float = 20.0,
                 eta_m: float = 20.0):
        """
        Initialize NSGA-II optimizer.
        
        Args:
            bounds: Parameter bounds
            population_size: Size of population
            crossover_prob: Crossover probability
            mutation_prob: Mutation probability  
            eta_c: Crossover distribution index
            eta_m: Mutation distribution index
        """
        self.bounds = np.array(bounds)
        self.dim = len(bounds)
        self.pop_size = population_size
        self.pc = crossover_prob
        self.pm = mutation_prob
        self.eta_c = eta_c
        self.eta_m = eta_m
        
    def _initialize_population(self) -> np.ndarray:
        """Initialize random population."""
        return np.random.uniform(
            self.bounds[:, 0], 
            self.bounds[:, 1],
            (self.pop_size, self.dim)
        )
    
    def _fast_non_dominated_sort(self, objectives: np.ndarray) -> List[List[int]]:
        """
        Fast non-dominated sorting algorithm.
        
        Args:
            objectives: Objective values (pop_size, n_objectives)
            
        Returns:
            List of fronts (each front is list of individual indices)
        """
        n = len(objectives)
        domination_count = np.zeros(n, dtype=int)
        dominated_solutions = [[] for _ in range(n)]
        fronts = [[]]
        
        # Find domination relationships
        for i in range(n):
            for j in range(n):
                if i != j:
                    if self._dominates(objectives[i], objectives[j]):
                        dominated_solutions[i].append(j)
                    elif self._dominates(objectives[j], objectives[i]):
                        domination_count[i] += 1
            
            if domination_count[i] == 0:
                fronts[0].append(i)
        
        # Build subsequent fronts
        front_idx = 0
        while len(fronts[front_idx]) > 0:
            next_front = []
            for i in fronts[front_idx]:
                for j in dominated_solutions[i]:
                    domination_count[j] -= 1
                    if domination_count[j] == 0:
                        next_front.append(j)
            
            if next_front:
                fronts.append(next_front)
            front_idx += 1
            
        return fronts[:-1]  # Remove empty last front
    
    def _dominates(self, obj1: np.ndarray, obj2: np.ndarray) -> bool:
        """Check if obj1 dominates obj2 (assuming minimization)."""
        return np.all(obj1 <= obj2) and np.any(obj1 < obj2)
    
    def _crowding_distance(self, objectives: np.ndarray, front: List[int]) -> np.ndarray:
        """
        Calculate crowding distance for individuals in a front.
        
        Args:
            objectives: Objective values
            front: Indices of individuals in the front
            
        Returns:
            Crowding distances
        """
        n = len(front)
        if n <= 2:
            return np.full(n, np.inf)
            
        distances = np.zeros(n)
        n_objectives = objectives.shape[1]
        
        for m in range(n_objectives):
            # Sort by m-th objective
            sorted_indices = np.argsort([objectives[front[i], m] for i in range(n)])
            
            # Boundary points get infinite distance
            distances[sorted_indices[0]] = np.inf
            distances[sorted_indices[-1]] = np.inf
            
            # Calculate distances for interior points
            obj_range = objectives[front[sorted_indices[-1]], m] - objectives[front[sorted_indices[0]], m]
            
            if obj_range > 0:
                for i in range(1, n-1):
                    distances[sorted_indices[i]] += (
                        objectives[front[sorted_indices[i+1]], m] - 
                        objectives[front[sorted_indices[i-1]], m]
                    ) / obj_range
                    
        return distances
    
    def _tournament_selection(self, population: np.ndarray, 
                            fronts: List[List[int]], 
                            crowding_distances: Dict[int, float]) -> np.ndarray:
        """Tournament selection based on dominance and crowding distance."""
        selected = []
        
        for _ in range(self.pop_size):
            # Select two random individuals
            i, j = np.random.choice(len(population), 2, replace=False)
            
            # Find which fronts they belong to
            front_i = front_j = -1
            for f_idx, front in enumerate(fronts):
                if i in front:
                    front_i = f_idx
                if j in front:
                    front_j = f_idx
            
            # Selection based on dominance rank and crowding distance
            if front_i < front_j:
                selected.append(population[i])
            elif front_j < front_i:
                selected.append(population[j])
            else:
                # Same front, use crowding distance
                if crowding_distances[i] > crowding_distances[j]:
                    selected.append(population[i])
                else:
                    selected.append(population[j])
                    
        return np.array(selected)
    
    def _simulated_binary_crossover(self, parent1: np.ndarray, 
                                   parent2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Simulated binary crossover operator."""
        child1 = parent1.copy()
        child2 = parent2.copy()
        
        for i in range(self.dim):
            if np.random.rand() <= 0.5:
                if abs(parent1[i] - parent2[i]) > 1e-14:
                    y1 = min(parent1[i], parent2[i])
                    y2 = max(parent1[i], parent2[i])
                    
                    lb = self.bounds[i, 0]
                    ub = self.bounds[i, 1]
                    
                    rand = np.random.rand()
                    
                    beta1 = 1.0 + (2.0 * (y1 - lb) / (y2 - y1))
                    beta2 = 1.0 + (2.0 * (ub - y2) / (y2 - y1))
                    
                    alpha1 = 2.0 - beta1**(-self.eta_c - 1.0)
                    alpha2 = 2.0 - beta2**(-self.eta_c - 1.0)
                    
                    if rand <= 1.0 / alpha1:
                        betaq1 = (rand * alpha1)**(1.0 / (self.eta_c + 1.0))
                    else:
                        betaq1 = (1.0 / (2.0 - rand * alpha1))**(1.0 / (self.eta_c + 1.0))
                        
                    if rand <= 1.0 / alpha2:
                        betaq2 = (rand * alpha2)**(1.0 / (self.eta_c + 1.0))
                    else:
                        betaq2 = (1.0 / (2.0 - rand * alpha2))**(1.0 / (self.eta_c + 1.0))
                    
                    child1[i] = 0.5 * ((y1 + y2) - betaq1 * (y2 - y1))
                    child2[i] = 0.5 * ((y1 + y2) + betaq2 * (y2 - y1))
                    
                    # Ensure bounds
                    child1[i] = np.clip(child1[i], lb, ub)
                    child2[i] = np.clip(child2[i], lb, ub)
                    
        return child1, child2
    
    def _polynomial_mutation(self, individual: np.ndarray) -> np.ndarray:
        """Polynomial mutation operator."""
        mutated = individual.copy()
        
        for i in range(self.dim):
            if np.random.rand() <= self.pm:
                y = mutated[i]
                lb = self.bounds[i, 0]
                ub = self.bounds[i, 1]
                
                delta1 = (y - lb) / (ub - lb)
                delta2 = (ub - y) / (ub - lb)
                
                rand = np.random.rand()
                mut_pow = 1.0 / (self.eta_m + 1.0)
                
                if rand <= 0.5:
                    xy = 1.0 - delta1
                    val = 2.0 * rand + (1.0 - 2.0 * rand) * xy**(self.eta_m + 1.0)
                    deltaq = val**mut_pow - 1.0
                else:
                    xy = 1.0 - delta2
                    val = 2.0 * (1.0 - rand) + 2.0 * (rand - 0.5) * xy**(self.eta_m + 1.0)
                    deltaq = 1.0 - val**mut_pow
                
                y = y + deltaq * (ub - lb)
                mutated[i] = np.clip(y, lb, ub)
                
        return mutated
    
    def optimize(self, objective_functions: List[Callable], 
                n_generations: int = 100) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run NSGA-II optimization.
        
        Args:
            objective_functions: List of objective functions to minimize
            n_generations: Number of generations
            
        Returns:
            pareto_front: Pareto optimal solutions
            pareto_objectives: Corresponding objective values
        """
        # Initialize population
        population = self._initialize_population()
        
        for generation in range(n_generations):
            # Evaluate objectives
            objectives = np.zeros((len(population), len(objective_functions)))
            for i, individual in enumerate(population):
                for j, obj_func in enumerate(objective_functions):
                    objectives[i, j] = obj_func(individual)
            
            # Non-dominated sorting
            fronts = self._fast_non_dominated_sort(objectives)
            
            # Calculate crowding distances
            crowding_distances = {}
            for front in fronts:
                distances = self._crowding_distance(objectives, front)
                for i, idx in enumerate(front):
                    crowding_distances[idx] = distances[i]
            
            # Selection for next generation
            new_population = []
            front_idx = 0
            
            while len(new_population) + len(fronts[front_idx]) <= self.pop_size:
                new_population.extend([population[i] for i in fronts[front_idx]])
                front_idx += 1
                
            # Fill remaining slots using crowding distance
            if len(new_population) < self.pop_size:
                remaining = self.pop_size - len(new_population)
                last_front = fronts[front_idx]
                
                # Sort by crowding distance
                sorted_indices = sorted(last_front, 
                                      key=lambda x: crowding_distances[x], 
                                      reverse=True)
                
                new_population.extend([population[i] for i in sorted_indices[:remaining]])
            
            population = np.array(new_population)
            
            # Generate offspring through crossover and mutation
            offspring = []
            for i in range(0, self.pop_size, 2):
                parent1 = population[i]
                parent2 = population[(i + 1) % self.pop_size]
                
                if np.random.rand() <= self.pc:
                    child1, child2 = self._simulated_binary_crossover(parent1, parent2)
                else:
                    child1, child2 = parent1.copy(), parent2.copy()
                
                child1 = self._polynomial_mutation(child1)
                child2 = self._polynomial_mutation(child2)
                
                offspring.extend([child1, child2])
            
            # Combine parent and offspring populations
            combined_pop = np.vstack([population, offspring[:self.pop_size]])
            population = combined_pop
        
        # Final evaluation and return Pareto front
        final_objectives = np.zeros((len(population), len(objective_functions)))
        for i, individual in enumerate(population):
            for j, obj_func in enumerate(objective_functions):
                final_objectives[i, j] = obj_func(individual)
        
        fronts = self._fast_non_dominated_sort(final_objectives)
        pareto_indices = fronts[0]
        
        return population[pareto_indices], final_objectives[pareto_indices]


class ResponseSurfaceMethod:
    """
    Response Surface Methodology for experimental design and optimization.
    """
    
    def __init__(self, bounds: List[Tuple[float, float]]):
        """
        Initialize response surface method.
        
        Args:
            bounds: Parameter bounds
        """
        self.bounds = np.array(bounds)
        self.dim = len(bounds)
        
    def central_composite_design(self, alpha: Optional[float] = None) -> np.ndarray:
        """
        Generate Central Composite Design points.
        
        Args:
            alpha: Axial distance (default: face-centered)
            
        Returns:
            Design matrix
        """
        if alpha is None:
            alpha = 1.0  # Face-centered design
            
        # Factorial points (2^k)
        factorial_points = []
        for i in range(2**self.dim):
            point = []
            for j in range(self.dim):
                if (i >> j) & 1:
                    point.append(1.0)
                else:
                    point.append(-1.0)
            factorial_points.append(point)
        
        # Axial points (2*k)
        axial_points = []
        for i in range(self.dim):
            # Positive axial point
            point_pos = [0.0] * self.dim
            point_pos[i] = alpha
            axial_points.append(point_pos)
            
            # Negative axial point
            point_neg = [0.0] * self.dim
            point_neg[i] = -alpha
            axial_points.append(point_neg)
        
        # Center point
        center_point = [[0.0] * self.dim]
        
        # Combine all points
        design_points = np.array(factorial_points + axial_points + center_point)
        
        # Transform to actual bounds
        design_actual = np.zeros_like(design_points)
        for i in range(self.dim):
            center = (self.bounds[i, 1] + self.bounds[i, 0]) / 2
            range_half = (self.bounds[i, 1] - self.bounds[i, 0]) / 2
            design_actual[:, i] = center + design_points[:, i] * range_half
            
        return design_actual
    
    def fit_quadratic_model(self, X: np.ndarray, y: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Fit quadratic response surface model.
        
        Args:
            X: Design points
            y: Response values
            
        Returns:
            Model coefficients
        """
        n, p = X.shape
        
        # Build design matrix for quadratic model
        # y = β₀ + Σβᵢxᵢ + Σβᵢᵢxᵢ² + ΣΣβᵢⱼxᵢxⱼ
        
        # Number of terms: 1 + p + p + p(p-1)/2
        n_terms = 1 + p + p + p * (p - 1) // 2
        design_matrix = np.zeros((n, n_terms))
        
        col = 0
        
        # Intercept
        design_matrix[:, col] = 1.0
        col += 1
        
        # Linear terms
        for i in range(p):
            design_matrix[:, col] = X[:, i]
            col += 1
        
        # Quadratic terms
        for i in range(p):
            design_matrix[:, col] = X[:, i]**2
            col += 1
        
        # Interaction terms
        for i in range(p):
            for j in range(i + 1, p):
                design_matrix[:, col] = X[:, i] * X[:, j]
                col += 1
        
        # Least squares fit
        coefficients = np.linalg.lstsq(design_matrix, y, rcond=None)[0]
        
        return {
            'coefficients': coefficients,
            'design_matrix': design_matrix,
            'r_squared': 1 - np.sum((y - design_matrix @ coefficients)**2) / np.sum((y - np.mean(y))**2)
        }
    
    def predict(self, model: Dict[str, np.ndarray], X_new: np.ndarray) -> np.ndarray:
        """
        Predict response at new points using fitted model.
        
        Args:
            model: Fitted model from fit_quadratic_model
            X_new: New points to predict
            
        Returns:
            Predicted responses
        """
        n, p = X_new.shape
        n_terms = len(model['coefficients'])
        
        # Build design matrix for new points
        design_matrix = np.zeros((n, n_terms))
        col = 0
        
        # Intercept
        design_matrix[:, col] = 1.0
        col += 1
        
        # Linear terms
        for i in range(p):
            design_matrix[:, col] = X_new[:, i]
            col += 1
        
        # Quadratic terms
        for i in range(p):
            design_matrix[:, col] = X_new[:, i]**2
            col += 1
        
        # Interaction terms
        for i in range(p):
            for j in range(i + 1, p):
                design_matrix[:, col] = X_new[:, i] * X_new[:, j]
                col += 1
        
        return design_matrix @ model['coefficients']


class ParetoFrontCalculator:
    """
    Utilities for Pareto front analysis and visualization.
    """
    
    @staticmethod
    def is_pareto_efficient(objectives: np.ndarray) -> np.ndarray:
        """
        Find Pareto efficient points.
        
        Args:
            objectives: Objective values (n_points, n_objectives)
            
        Returns:
            Boolean array indicating Pareto efficient points
        """
        n_points = objectives.shape[0]
        is_efficient = np.ones(n_points, dtype=bool)
        
        for i in range(n_points):
            if is_efficient[i]:
                # Check if any other point dominates point i
                dominated = np.all(objectives <= objectives[i], axis=1) & \
                           np.any(objectives < objectives[i], axis=1)
                is_efficient[dominated] = False
                
        return is_efficient
    
    @staticmethod
    def hypervolume(pareto_front: np.ndarray, reference_point: np.ndarray) -> float:
        """
        Calculate hypervolume indicator.
        
        Args:
            pareto_front: Pareto front points
            reference_point: Reference point for hypervolume calculation
            
        Returns:
            Hypervolume value
        """
        if pareto_front.shape[1] != 2:
            raise NotImplementedError("Hypervolume calculation only implemented for 2D")
        
        # Sort by first objective
        sorted_indices = np.argsort(pareto_front[:, 0])
        sorted_front = pareto_front[sorted_indices]
        
        hypervolume = 0.0
        prev_x = reference_point[0]
        
        for point in sorted_front:
            if point[0] > prev_x and point[1] < reference_point[1]:
                width = point[0] - prev_x
                height = reference_point[1] - point[1]
                hypervolume += width * height
                prev_x = point[0]
                
        return hypervolume
    
    @staticmethod
    def crowding_distance(pareto_front: np.ndarray) -> np.ndarray:
        """
        Calculate crowding distance for Pareto front points.
        
        Args:
            pareto_front: Pareto front points
            
        Returns:
            Crowding distances
        """
        n_points, n_objectives = pareto_front.shape
        distances = np.zeros(n_points)
        
        for m in range(n_objectives):
            # Sort by m-th objective
            sorted_indices = np.argsort(pareto_front[:, m])
            
            # Boundary points get infinite distance
            distances[sorted_indices[0]] = np.inf
            distances[sorted_indices[-1]] = np.inf
            
            # Calculate distances for interior points
            obj_range = pareto_front[sorted_indices[-1], m] - pareto_front[sorted_indices[0], m]
            
            if obj_range > 0:
                for i in range(1, n_points - 1):
                    distances[sorted_indices[i]] += (
                        pareto_front[sorted_indices[i + 1], m] - 
                        pareto_front[sorted_indices[i - 1], m]
                    ) / obj_range
                    
        return distances


class SobolSensitivityAnalyzer:
    """
    Sobol sensitivity analysis for global sensitivity analysis.
    """
    
    def __init__(self, bounds: List[Tuple[float, float]]):
        """
        Initialize Sobol analyzer.
        
        Args:
            bounds: Parameter bounds
        """
        self.bounds = np.array(bounds)
        self.dim = len(bounds)
        
    def generate_samples(self, n_samples: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate Sobol sample matrices.
        
        Args:
            n_samples: Number of samples
            
        Returns:
            A, B, C matrices for Sobol analysis
        """
        # Generate Sobol sequences
        sampler = qmc.Sobol(d=self.dim, scramble=True)
        
        # Generate base samples
        sobol_A = sampler.random(n_samples)
        sobol_B = sampler.random(n_samples)
        
        # Scale to bounds
        A = qmc.scale(sobol_A, self.bounds[:, 0], self.bounds[:, 1])
        B = qmc.scale(sobol_B, self.bounds[:, 0], self.bounds[:, 1])
        
        # Generate C matrices (A with i-th column from B)
        C = np.zeros((self.dim, n_samples, self.dim))
        for i in range(self.dim):
            C[i] = A.copy()
            C[i][:, i] = B[:, i]
            
        return A, B, C
    
    def compute_indices(self, model_func: Callable, n_samples: int = 1000) -> Dict[str, np.ndarray]:
        """
        Compute Sobol sensitivity indices.
        
        Args:
            model_func: Model function to analyze
            n_samples: Number of samples for estimation
            
        Returns:
            Dictionary with first-order and total-order indices
        """
        # Generate sample matrices
        A, B, C = self.generate_samples(n_samples)
        
        # Evaluate model
        f_A = np.array([model_func(x) for x in A])
        f_B = np.array([model_func(x) for x in B])
        f_C = np.zeros((self.dim, n_samples))
        
        for i in range(self.dim):
            f_C[i] = np.array([model_func(x) for x in C[i]])
        
        # Compute variance
        f_total = np.concatenate([f_A, f_B] + [f_C[i] for i in range(self.dim)])
        V = np.var(f_total)
        
        # First-order indices
        S1 = np.zeros(self.dim)
        for i in range(self.dim):
            S1[i] = np.mean(f_B * (f_C[i] - f_A)) / V
        
        # Total-order indices
        ST = np.zeros(self.dim)
        for i in range(self.dim):
            ST[i] = 1 - np.mean(f_A * (f_C[i] - f_B)) / V
        
        return {
            'first_order': S1,
            'total_order': ST,
            'variance': V
        }
