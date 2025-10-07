"""
Liver Fibrosis Simulation Package
=================================

A comprehensive PhysiCell-based simulation platform for studying hepatic stellate cell 
activation and miRNA-based therapeutic interventions in liver fibrosis.

Main modules:
- simulation: Core PhysiCell simulation interface
- analysis: Data analysis and statistical processing
- visualization: Professional scientific visualization
- data_generation: Synthetic data generation utilities
- config: Configuration management

Example usage:
    >>> from liver_fibrosis_sim import LiverFibrosisSimulation
    >>> sim = LiverFibrosisSimulation()
    >>> sim.run_demo()
"""

__version__ = "1.0.0"
__author__ = "iGEM Liver Fibrosis Team"
__email__ = "contact@liverfibrosis.simulation"

# Import main classes for easy access
try:
    from .simulation import LiverFibrosisSimulation
    from .analysis import EnhancedAnalyzer
    from .visualization import ProfessionalVisualization, AdvancedVisualization
    from .data_generation import DataGenerator
    
    __all__ = [
        'LiverFibrosisSimulation',
        'EnhancedAnalyzer', 
        'ProfessionalVisualization',
        'AdvancedVisualization',
        'DataGenerator'
    ]
    
except ImportError as e:
    # Handle case where some dependencies might be missing
    print(f"Warning: Some modules could not be imported: {e}")
    __all__ = []

# Package metadata
PACKAGE_INFO = {
    'name': 'liver-fibrosis-simulation',
    'version': __version__,
    'description': 'Advanced PhysiCell-based liver fibrosis simulation',
    'author': __author__,
    'email': __email__,
    'url': 'https://github.com/igem/liver-fibrosis-simulation'
}
