"""
Biological Models Module for Liver Fibrosis Simulation Platform

This module implements advanced biological mechanism models for multi-scale
cellular and molecular system modeling, including:
- Cellular dynamics and behavior modules
- Molecular interaction networks  
- Membrane transport systems
- Signal transduction pathways
- Gene regulatory networks

Author: Liver Fibrosis Simulation Team
License: MIT
"""

from .cellular_dynamics import (
    CellCycleModel,
    MetabolismModel,
    StressResponseModel,
    ReceptorDynamicsModel,
    CytoskeletonModel,
    MigrationModel,
    ApoptosisModel,
    SenescenceModel,
    CellularBehaviorManager
)

from .molecular_networks import (
    TGFBetaSignalingNetwork,
    MiRNARegulationNetwork,
    BooleanContinuousModel,
    SynergyEffectCalculator,
    MasterEquationSolver,
    SignalTransductionPathway
)

from .membrane_transport import (
    EndocytosisModel,
    ExocytosisModel,
    IntracellularTraffickingModel,
    MembraneReceptorModel,
    VesicleTransportModel,
    AutophagyModel
)

__version__ = "1.0.0"
__author__ = "Liver Fibrosis Simulation Team"

__all__ = [
    # Cellular Dynamics
    "CellCycleModel",
    "MetabolismModel", 
    "StressResponseModel",
    "ReceptorDynamicsModel",
    "CytoskeletonModel",
    "MigrationModel",
    "ApoptosisModel",
    "SenescenceModel",
    "CellularBehaviorManager",
    
    # Molecular Networks
    "TGFBetaSignalingNetwork",
    "MiRNARegulationNetwork",
    "BooleanContinuousModel",
    "SynergyEffectCalculator", 
    "MasterEquationSolver",
    "SignalTransductionPathway",
    
    # Membrane Transport
    "EndocytosisModel",
    "ExocytosisModel",
    "IntracellularTraffickingModel",
    "MembraneReceptorModel",
    "VesicleTransportModel",
    "AutophagyModel"
]
