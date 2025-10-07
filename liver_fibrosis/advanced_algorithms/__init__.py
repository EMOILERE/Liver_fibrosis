"""
Advanced Algorithms Module for Liver Fibrosis Simulation Platform

This module implements advanced computational algorithms for high-performance
biological system modeling including parallel computing, data analysis,
and parameter calibration techniques.

Key Features:
- Parallel computing frameworks with load balancing
- High-dimensional data analysis and dimensionality reduction
- Intelligent parameter calibration and optimization
- Machine learning integration for biological modeling
- Advanced statistical analysis and uncertainty quantification
- High-performance computing optimizations

Author: Liver Fibrosis Simulation Team
License: MIT
"""

from .parallel_computing import (
    ParallelComputingFramework,
    LoadBalancer,
    DistributedMemoryManager,
    CommunicationOptimizer,
    TaskScheduler,
    PerformanceMonitor
)

from .data_analysis import (
    DimensionalityReducer,
    TSNEAnalyzer,
    PCAAnalyzer,
    GaussianMixtureAnalyzer,
    SpatialCorrelationAnalyzer,
    TimeSeriesAnalyzer,
    StatisticalTestSuite
)

from .parameter_calibration import (
    BayesianParameterInference,
    MarkovChainMonteCarlo,
    AutomaticDifferentiation,
    UncertaintyQuantification,
    ModelSelection,
    SensitivityAnalyzer
)

__version__ = "1.0.0"
__author__ = "Liver Fibrosis Simulation Team"

__all__ = [
    # Parallel Computing
    "ParallelComputingFramework",
    "LoadBalancer",
    "DistributedMemoryManager",
    "CommunicationOptimizer",
    "TaskScheduler",
    "PerformanceMonitor",
    
    # Data Analysis
    "DimensionalityReducer",
    "TSNEAnalyzer",
    "PCAAnalyzer",
    "GaussianMixtureAnalyzer",
    "SpatialCorrelationAnalyzer",
    "TimeSeriesAnalyzer",
    "StatisticalTestSuite",
    
    # Parameter Calibration
    "BayesianParameterInference",
    "MarkovChainMonteCarlo",
    "AutomaticDifferentiation",
    "UncertaintyQuantification",
    "ModelSelection",
    "SensitivityAnalyzer"
]
