"""
Data Generation Module for Liver Fibrosis Simulation
Synthetic data generation for testing and demonstration
"""

from pathlib import Path
import sys
from typing import Optional, Union, Dict, Any

# Add parent directory to path to import data generation modules
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

class DataGenerator:
    """Wrapper for data generation functionality"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize data generator"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
        
        # Import and initialize the actual generator class
        try:
            from generate_ideal_data import IdealDataGenerator
            self.generator = IdealDataGenerator(str(self.base_dir))
        except ImportError as e:
            print(f"Warning: Could not import data generator: {e}")
            self.generator = None
    
    def generate_3d_cell_positions(self, n_cells: int, group_params: Dict):
        """Generate 3D cell positions"""
        if self.generator:
            return self.generator.generate_3d_cell_positions(n_cells, group_params)
        return None
    
    def generate_cellular_states(self, positions, group_params: Dict, time_point: float):
        """Generate cellular state data"""
        if self.generator:
            return self.generator.generate_cellular_states(positions, group_params, time_point)
        return None
    
    def generate_microenvironment_data(self, group_params: Dict, time_point: float):
        """Generate microenvironment data"""
        if self.generator:
            return self.generator.generate_microenvironment_data(group_params, time_point)
        return None
    
    def generate_time_series_analysis(self, group_name: str, group_params: Dict):
        """Generate time series analysis data"""
        if self.generator:
            return self.generator.generate_time_series_analysis(group_name, group_params)
        return None
    
    def generate_complete_dataset(self, group_name: Optional[str] = None, save_xml: bool = True):
        """Generate complete dataset"""
        if self.generator:
            return self.generator.generate_complete_dataset(group_name, save_xml)
        return False
    
    def create_summary_visualization(self):
        """Create data overview visualization"""
        if self.generator:
            return self.generator.create_summary_visualization()
        return False
    
    def run_complete_data_generation(self):
        """Run complete data generation workflow"""
        if self.generator:
            return self.generator.run_complete_data_generation()
        return False

class ExperimentalDesign:
    """Experimental design utilities"""
    
    def __init__(self):
        """Initialize experimental design"""
        self.experimental_groups = {
            'control': {
                'TGF_beta1': 2.0,
                'exosome_concentration': 0.0,
                'miR_455_3p_baseline': 0.0,
                'miR_148a_5p_baseline': 0.0,
                'expected_activation': 0.08,
                'treatment_efficacy': 0.0,
                'synergy_factor': 1.0
            },
            'positive_model': {
                'TGF_beta1': 15.0,
                'exosome_concentration': 0.0,
                'miR_455_3p_baseline': 0.0,
                'miR_148a_5p_baseline': 0.0,
                'expected_activation': 0.75,
                'treatment_efficacy': 0.0,
                'synergy_factor': 1.0
            },
            'miR_455_3p': {
                'TGF_beta1': 15.0,
                'exosome_concentration': 50.0,
                'miR_455_3p_baseline': 0.4,
                'miR_148a_5p_baseline': 0.0,
                'expected_activation': 0.45,
                'treatment_efficacy': 42.0,
                'synergy_factor': 1.2
            },
            'miR_148a_5p': {
                'TGF_beta1': 15.0,
                'exosome_concentration': 50.0,
                'miR_455_3p_baseline': 0.0,
                'miR_148a_5p_baseline': 0.4,
                'expected_activation': 0.50,
                'treatment_efficacy': 38.0,
                'synergy_factor': 1.2
            },
            'dual_miRNA_1_1': {
                'TGF_beta1': 15.0,
                'exosome_concentration': 75.0,
                'miR_455_3p_baseline': 0.25,
                'miR_148a_5p_baseline': 0.25,
                'expected_activation': 0.28,
                'treatment_efficacy': 68.0,
                'synergy_factor': 2.8
            }
        }
    
    def get_group_parameters(self, group_name: str) -> Dict:
        """Get parameters for experimental group"""
        return self.experimental_groups.get(group_name, {})
    
    def list_experimental_groups(self):
        """List all available experimental groups"""
        return list(self.experimental_groups.keys())
    
    def create_custom_group(self, name: str, parameters: Dict):
        """Create custom experimental group"""
        self.experimental_groups[name] = parameters
    
    def get_design_matrix(self):
        """Get experimental design matrix"""
        try:
            import pandas as pd
            
            design_data = []
            for group, params in self.experimental_groups.items():
                design_data.append({
                    'Group': group,
                    **params
                })
            
            return pd.DataFrame(design_data)
        except ImportError:
            print("pandas not available for design matrix")
            return None

class ParameterSpace:
    """Parameter space exploration utilities"""
    
    def __init__(self):
        """Initialize parameter space"""
        self.parameter_ranges = {
            'TGF_beta1': (0.5, 25.0),
            'exosome_concentration': (0.0, 100.0),
            'miR_455_3p_baseline': (0.0, 1.0),
            'miR_148a_5p_baseline': (0.0, 1.0),
            'synergy_factor': (1.0, 5.0)
        }
    
    def sample_parameter_space(self, n_samples: int = 100):
        """Sample parameter space for sensitivity analysis"""
        try:
            import numpy as np
            
            samples = {}
            for param, (min_val, max_val) in self.parameter_ranges.items():
                samples[param] = np.random.uniform(min_val, max_val, n_samples)
            
            return samples
        except ImportError:
            print("numpy not available for parameter sampling")
            return None
    
    def create_latin_hypercube_design(self, n_samples: int = 100):
        """Create Latin hypercube sampling design"""
        try:
            from sklearn.preprocessing import StandardScaler
            import numpy as np
            
            # Simple Latin hypercube implementation
            n_params = len(self.parameter_ranges)
            samples = np.zeros((n_samples, n_params))
            
            for i, (param, (min_val, max_val)) in enumerate(self.parameter_ranges.items()):
                # Create Latin hypercube samples
                intervals = np.linspace(0, 1, n_samples + 1)
                lower = intervals[:-1]
                upper = intervals[1:]
                points = np.random.uniform(lower, upper)
                np.random.shuffle(points)
                
                # Scale to parameter range
                samples[:, i] = points * (max_val - min_val) + min_val
            
            return samples
        except ImportError:
            print("sklearn/numpy not available for Latin hypercube sampling")
            return None

def main():
    """Command-line interface for data generation"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Liver Fibrosis Data Generation')
    parser.add_argument('--generate', action='store_true', help='Generate synthetic data')
    parser.add_argument('--group', help='Generate data for specific group')
    parser.add_argument('--all-groups', action='store_true', help='Generate data for all groups')
    parser.add_argument('--visualization', action='store_true', help='Create data visualization')
    
    args = parser.parse_args()
    
    success = True
    
    if args.generate or args.all_groups:
        print("Generating synthetic data...")
        generator = DataGenerator()
        if args.group:
            success &= generator.generate_complete_dataset(args.group)
        else:
            success &= generator.run_complete_data_generation()
    
    if args.visualization:
        print("Creating data visualization...")
        generator = DataGenerator()
        success &= generator.create_summary_visualization()
    
    if not any([args.generate, args.all_groups, args.visualization]):
        print("No operation specified. Use --help for options.")
        success = False
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
