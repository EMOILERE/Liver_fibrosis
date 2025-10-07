"""
Mathematical Core Module for Liver Fibrosis Simulation Platform

This module provides advanced mathematical algorithms and computational frameworks
for multi-scale biological system modeling, including:
- Multi-scale numerical solvers
- Stochastic process computations  
- Optimization algorithms
- Error estimation and adaptive control

Author: Liver Fibrosis Simulation Team
License: MIT
"""

from .multiscale_solvers import (
    ADISolver,
    RungeKuttaFehlbergSolver, 
    OperatorSplittingSolver,
    MultipleTimeSteppingSolver,
    AdaptiveTimeStepController
)

from .stochastic_processes import (
    MetropolisHastingsSampler,
    GillespieSolver,
    OrnsteinUhlenbeckProcess,
    WienerProcessGenerator,
    StochasticDifferentialEquationSolver
)

from .optimization_algorithms import (
    BayesianOptimizer,
    NSGA2Optimizer,
    ResponseSurfaceMethod,
    ParetoFrontCalculator,
    SobolSensitivityAnalyzer
)

__version__ = "1.0.0"
__author__ = "Liver Fibrosis Simulation Team"

__all__ = [
    # Solvers
    "ADISolver",
    "RungeKuttaFehlbergSolver",
    "OperatorSplittingSolver", 
    "MultipleTimeSteppingSolver",
    "AdaptiveTimeStepController",
    
    # Stochastic Processes
    "MetropolisHastingsSampler",
    "GillespieSolver",
    "OrnsteinUhlenbeckProcess", 
    "WienerProcessGenerator",
    "StochasticDifferentialEquationSolver",
    
    # Optimization
    "BayesianOptimizer",
    "NSGA2Optimizer",
    "ResponseSurfaceMethod",
    "ParetoFrontCalculator", 
    "SobolSensitivityAnalyzer"
]
