#!/usr/bin/env python3
"""
Experiment Configuration Manager for Liver Fibrosis miRNA Treatment Study
Manages all 9 experimental configurations and their parameters
"""

from dataclasses import dataclass
from typing import Dict, List, Optional
from pathlib import Path

@dataclass
class ExperimentConfig:
    """Configuration for a single experiment"""
    name: str
    display_name: str
    description: str
    miR_455_3p_dose: float
    miR_148a_5p_dose: float
    has_exosomes: bool
    is_control: bool
    color: str
    marker_style: str

class ExperimentConfigManager:
    """Manager for all experimental configurations"""
    
    def __init__(self):
        """Initialize with all 9 experimental configurations"""
        self.configs = {
            'control': ExperimentConfig(
                name='control',
                display_name='Control',
                description='Negative control - no treatment',
                miR_455_3p_dose=0.0,
                miR_148a_5p_dose=0.0,
                has_exosomes=False,
                is_control=True,
                color='#95A5A6',  # Gray
                marker_style='o'
            ),
            
            'dual_miRNA_1_1': ExperimentConfig(
                name='dual_miRNA_1_1',
                display_name='Dual miRNA (1:1)',
                description='miR-455-3p and miR-148a-5p at 1:1 ratio',
                miR_455_3p_dose=1.0,
                miR_148a_5p_dose=1.0,
                has_exosomes=True,
                is_control=False,
                color='#9B59B6',  # Purple
                marker_style='s'
            ),
            
            'dual_miRNA_2_1': ExperimentConfig(
                name='dual_miRNA_2_1',
                display_name='Dual miRNA (2:1)',
                description='miR-455-3p and miR-148a-5p at 2:1 ratio',
                miR_455_3p_dose=2.0,
                miR_148a_5p_dose=1.0,
                has_exosomes=True,
                is_control=False,
                color='#8E44AD',  # Dark purple
                marker_style='D'
            ),
            
            'dual_miRNA_1_2': ExperimentConfig(
                name='dual_miRNA_1_2',
                display_name='Dual miRNA (1:2)',
                description='miR-455-3p and miR-148a-5p at 1:2 ratio',
                miR_455_3p_dose=1.0,
                miR_148a_5p_dose=2.0,
                has_exosomes=True,
                is_control=False,
                color='#BB8FCE',  # Light purple
                marker_style='^'
            ),
            
            'miR_455_only': ExperimentConfig(
                name='miR_455_only',
                display_name='miR-455-3p Only',
                description='Single miR-455-3p treatment',
                miR_455_3p_dose=1.0,
                miR_148a_5p_dose=0.0,
                has_exosomes=True,
                is_control=False,
                color='#E74C3C',  # Red
                marker_style='v'
            ),
            
            'miR_148_only': ExperimentConfig(
                name='miR_148_only',
                display_name='miR-148a-5p Only',
                description='Single miR-148a-5p treatment',
                miR_455_3p_dose=0.0,
                miR_148a_5p_dose=1.0,
                has_exosomes=True,
                is_control=False,
                color='#3498DB',  # Blue
                marker_style='<'
            ),
            
            'natural_exosomes': ExperimentConfig(
                name='natural_exosomes',
                display_name='Natural Exosomes',
                description='Natural exosomes without miRNA loading',
                miR_455_3p_dose=0.0,
                miR_148a_5p_dose=0.0,
                has_exosomes=True,
                is_control=False,
                color='#2ECC71',  # Green
                marker_style='>'
            ),
            
            'negative_control': ExperimentConfig(
                name='negative_control',
                display_name='NC Mimic',
                description='Negative control mimic - scrambled sequence',
                miR_455_3p_dose=0.0,
                miR_148a_5p_dose=0.0,
                has_exosomes=True,
                is_control=True,
                color='#F39C12',  # Orange
                marker_style='h'
            ),
            
            'positive_model': ExperimentConfig(
                name='positive_model',
                display_name='Positive Control',
                description='Positive control - maximum fibrosis induction',
                miR_455_3p_dose=0.0,
                miR_148a_5p_dose=0.0,
                has_exosomes=False,
                is_control=False,
                color='#C0392B',  # Dark red
                marker_style='p'
            )
        }
    
    def get_config(self, name: str) -> Optional[ExperimentConfig]:
        """Get configuration by name"""
        return self.configs.get(name)
    
    def get_all_configs(self) -> Dict[str, ExperimentConfig]:
        """Get all configurations"""
        return self.configs.copy()
    
    def get_config_names(self) -> List[str]:
        """Get list of all configuration names"""
        return list(self.configs.keys())
    
    def get_treatment_groups(self) -> List[str]:
        """Get list of treatment groups (excluding controls)"""
        return [name for name, config in self.configs.items() if not config.is_control]
    
    def get_control_groups(self) -> List[str]:
        """Get list of control groups"""
        return [name for name, config in self.configs.items() if config.is_control]
    
    def get_miRNA_groups(self) -> List[str]:
        """Get list of groups with miRNA treatment"""
        return [name for name, config in self.configs.items() 
                if config.miR_455_3p_dose > 0 or config.miR_148a_5p_dose > 0]
    
    def get_exosome_groups(self) -> List[str]:
        """Get list of groups with exosome delivery"""
        return [name for name, config in self.configs.items() if config.has_exosomes]
    
    def create_output_directories(self, base_dir: Path) -> Dict[str, Path]:
        """Create output directories for all configurations"""
        directories = {}
        
        for config_name in self.configs:
            config_dir = base_dir / f"results_{config_name}"
            config_dir.mkdir(exist_ok=True)
            
            # Create subdirectories
            (config_dir / "figures").mkdir(exist_ok=True)
            (config_dir / "pdfs").mkdir(exist_ok=True)
            (config_dir / "animations").mkdir(exist_ok=True)
            (config_dir / "interactive").mkdir(exist_ok=True)
            (config_dir / "data").mkdir(exist_ok=True)
            
            directories[config_name] = config_dir
        
        return directories
    
    def get_comparison_pairs(self) -> List[tuple]:
        """Get meaningful comparison pairs for analysis"""
        return [
            ('control', 'dual_miRNA_1_1'),
            ('control', 'miR_455_only'),
            ('control', 'miR_148_only'),
            ('miR_455_only', 'miR_148_only'),
            ('dual_miRNA_1_1', 'dual_miRNA_2_1'),
            ('dual_miRNA_1_1', 'dual_miRNA_1_2'),
            ('natural_exosomes', 'negative_control'),
            ('control', 'positive_model')
        ]
    
    def get_dose_series(self) -> Dict[str, List[str]]:
        """Get dose series for dose-response analysis"""
        return {
            'miR_455_series': ['control', 'miR_455_only', 'dual_miRNA_1_1', 'dual_miRNA_2_1'],
            'miR_148_series': ['control', 'miR_148_only', 'dual_miRNA_1_1', 'dual_miRNA_1_2'],
            'dual_ratio_series': ['dual_miRNA_2_1', 'dual_miRNA_1_1', 'dual_miRNA_1_2']
        }
    
    def validate_data_availability(self, base_dir: Path) -> Dict[str, bool]:
        """Check which configurations have available data"""
        availability = {}
        
        if (base_dir / "synthetic_results").exists():
            synthetic_dir = base_dir / "synthetic_results"
            
            for config_name in self.configs:
                output_dir = synthetic_dir / f"output_{config_name}"
                time_series_file = synthetic_dir / f"time_series_{config_name}.csv"
                
                # Check if both directory and time series exist
                availability[config_name] = (
                    output_dir.exists() and 
                    time_series_file.exists() and
                    len(list(output_dir.glob("cells_*.csv"))) > 0
                )
        else:
            # If no synthetic data, assume all configurations are available for demo data
            print("No synthetic data found, using demo data for all configurations")
            availability = {name: True for name in self.configs}
        
        return availability
    
    def get_statistical_groups(self) -> Dict[str, List[str]]:
        """Get groups for statistical analysis"""
        return {
            'all_treatments': self.get_treatment_groups(),
            'all_controls': self.get_control_groups(),
            'miRNA_treatments': self.get_miRNA_groups(),
            'exosome_treatments': self.get_exosome_groups(),
            'single_miRNA': ['miR_455_only', 'miR_148_only'],
            'dual_miRNA': ['dual_miRNA_1_1', 'dual_miRNA_2_1', 'dual_miRNA_1_2']
        }

# Global instance
experiment_manager = ExperimentConfigManager()

# Convenience functions
def get_all_experiment_names() -> List[str]:
    """Get all experiment names"""
    return experiment_manager.get_config_names()

def get_experiment_config(name: str) -> Optional[ExperimentConfig]:
    """Get specific experiment configuration"""
    return experiment_manager.get_config(name)

def create_experiment_directories(base_dir: Path) -> Dict[str, Path]:
    """Create directories for all experiments"""
    return experiment_manager.create_output_directories(base_dir)

if __name__ == "__main__":
    # Test the configuration manager
    print("Liver Fibrosis Experiment Configuration Manager")
    print("=" * 50)
    
    manager = ExperimentConfigManager()
    
    print("Available Experiments:")
    for name, config in manager.get_all_configs().items():
        print(f"  {name}: {config.display_name}")
        print(f"    {config.description}")
        print(f"    miR-455: {config.miR_455_3p_dose}, miR-148: {config.miR_148a_5p_dose}")
        print()
    
    print("Treatment Groups:", manager.get_treatment_groups())
    print("Control Groups:", manager.get_control_groups())
    print("miRNA Groups:", manager.get_miRNA_groups())
    print("Exosome Groups:", manager.get_exosome_groups())
    
    print("\nComparison Pairs:")
    for pair in manager.get_comparison_pairs():
        print(f"  {pair[0]} vs {pair[1]}")
    
    print("\nDose Series:")
    for series_name, groups in manager.get_dose_series().items():
        print(f"  {series_name}: {groups}")
