"""
Multi-Scale Integration Module for Liver Fibrosis Simulation Platform

This module implements multi-scale coupling systems that integrate molecular,
cellular, and tissue-level processes for comprehensive biological modeling.

Key Features:
- Scale coupling frameworks for molecular-cellular-tissue integration
- Reaction-diffusion systems for molecular transport
- Tissue mechanics and mechanical signal transduction
- Multi-physics field coupling
- Scale separation and homogenization techniques
- Cross-scale information transfer mechanisms

Author: Liver Fibrosis Simulation Team
License: MIT
"""

from .scale_coupling import (
    MultiScaleCouplingFramework,
    ScaleSeparationHandler,
    CrossScaleInformationTransfer,
    HomogenizationTechniques,
    UpscalingOperator,
    DownscalingOperator
)

from .reaction_diffusion import (
    ReactionDiffusionSolver,
    MolecularSpeciesManager,
    DiffusionTensorCalculator,
    ReactionNetworkSolver,
    BoundaryConditionHandler,
    ConcentrationFieldUpdater
)

from .tissue_mechanics import (
    ContinuumMechanicsModel,
    CellDensityFieldModel,
    StressFieldCalculator,
    MechanicalSignalTransduction,
    ElasticViscousModel,
    CellMigrationMechanics
)

__version__ = "1.0.0"
__author__ = "Liver Fibrosis Simulation Team"

__all__ = [
    # Scale Coupling
    "MultiScaleCouplingFramework",
    "ScaleSeparationHandler",
    "CrossScaleInformationTransfer",
    "HomogenizationTechniques",
    "UpscalingOperator",
    "DownscalingOperator",
    
    # Reaction-Diffusion
    "ReactionDiffusionSolver",
    "MolecularSpeciesManager",
    "DiffusionTensorCalculator",
    "ReactionNetworkSolver",
    "BoundaryConditionHandler",
    "ConcentrationFieldUpdater",
    
    # Tissue Mechanics
    "ContinuumMechanicsModel",
    "CellDensityFieldModel", 
    "StressFieldCalculator",
    "MechanicalSignalTransduction",
    "ElasticViscousModel",
    "CellMigrationMechanics"
]
