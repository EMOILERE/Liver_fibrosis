"""
Visualization Module for Liver Fibrosis Simulation
Professional visualization and analysis tools
"""

from pathlib import Path
import sys
from typing import Optional, Union

# Add parent directory to path to import visualization modules
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

class ProfessionalVisualization:
    """Wrapper for professional visualization functionality"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize professional visualization"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
        
        # Import and initialize the actual visualization class
        try:
            from professional_visualization import ProfessionalBioVisualization
            self.visualizer = ProfessionalBioVisualization(str(self.base_dir))
        except ImportError as e:
            print(f"Warning: Could not import professional visualization: {e}")
            self.visualizer = None
    
    def create_morphological_analysis(self):
        """Create morphological analysis plots"""
        if self.visualizer:
            return self.visualizer.create_morphological_analysis([])
        return False
    
    def create_molecular_heatmaps(self):
        """Create molecular distribution heatmaps"""
        if self.visualizer:
            return self.visualizer.create_molecular_heatmaps()
        return False
    
    def create_dose_response_curves(self):
        """Create dose-response analysis"""
        if self.visualizer:
            return self.visualizer.create_dose_response_curves()
        return False
    
    def create_pathway_network(self):
        """Create pathway network diagram"""
        if self.visualizer:
            return self.visualizer.create_pathway_network_diagram()
        return False
    
    def create_statistical_analysis(self):
        """Create statistical analysis plots"""
        if self.visualizer:
            return self.visualizer.create_statistical_analysis()
        return False
    
    def run_complete_suite(self):
        """Run complete professional visualization suite"""
        if self.visualizer:
            return self.visualizer.run_professional_visualization_suite()
        return False

class AdvancedVisualization:
    """Wrapper for advanced visualization functionality including 3D and animations"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize advanced visualization"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
        
        # Import and initialize the actual visualization class
        try:
            from advanced_visualization import AdvancedVisualization as AdvViz
            self.visualizer = AdvViz(str(self.base_dir))
        except ImportError as e:
            print(f"Warning: Could not import advanced visualization: {e}")
            self.visualizer = None
    
    def create_3d_heatmaps(self):
        """Create 3D molecular heatmaps"""
        if self.visualizer:
            return self.visualizer.create_3d_heatmaps()
        return False
    
    def create_2d_dynamics_gif(self):
        """Create 2D cell dynamics animation"""
        if self.visualizer:
            return self.visualizer.create_2d_cell_dynamics_gif()
        return False
    
    def create_3d_dynamics_gif(self):
        """Create 3D cell dynamics animation"""
        if self.visualizer:
            return self.visualizer.create_3d_cell_dynamics_gif()
        return False
    
    def create_treatment_analysis(self):
        """Create treatment efficacy analysis"""
        if self.visualizer:
            return self.visualizer.create_treatment_efficacy_analysis()
        return False
    
    def create_morphological_analysis(self):
        """Create morphological analysis"""
        if self.visualizer:
            return self.visualizer.create_morphological_analysis()
        return False
    
    def run_complete_advanced_visualization(self):
        """Run complete advanced visualization suite"""
        if self.visualizer:
            return self.visualizer.run_complete_advanced_visualization()
        return False

class Visualization3D:
    """Wrapper for 3D-specific visualization functionality"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize 3D visualization"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
        
        # Import and initialize the actual 3D visualization class
        try:
            from threeD_visualization import Advanced3DVisualization
            self.visualizer = Advanced3DVisualization(str(self.base_dir))
        except ImportError:
            try:
                # Alternative import name
                from three_d_visualization import Advanced3DVisualization
                self.visualizer = Advanced3DVisualization(str(self.base_dir))
            except ImportError as e:
                print(f"Warning: Could not import 3D visualization: {e}")
                self.visualizer = None
    
    def create_3d_cell_visualization(self):
        """Create 3D cell visualization"""
        if self.visualizer:
            return self.visualizer.create_3D_cell_visualization()
        return False
    
    def create_3d_molecular_gradients(self):
        """Create 3D molecular gradient visualization"""
        if self.visualizer:
            return self.visualizer.create_3D_molecular_gradients()
        return False
    
    def analyze_spatial_patterns(self):
        """Analyze 3D spatial patterns"""
        if self.visualizer:
            return self.visualizer.analyze_3D_spatial_patterns()
        return False
    
    def run_complete_3d_analysis(self):
        """Run complete 3D analysis pipeline"""
        if self.visualizer:
            return self.visualizer.run_complete_3D_analysis()
        return False

def main():
    """Command-line interface for visualization"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Liver Fibrosis Visualization')
    parser.add_argument('--professional', action='store_true', help='Run professional visualization suite')
    parser.add_argument('--advanced', action='store_true', help='Run advanced visualization suite')
    parser.add_argument('--3d', action='store_true', help='Run 3D visualization suite')
    parser.add_argument('--all', action='store_true', help='Run all visualization suites')
    
    args = parser.parse_args()
    
    success = True
    
    if args.professional or args.all:
        print("Running professional visualization suite...")
        viz = ProfessionalVisualization()
        success &= viz.run_complete_suite()
    
    if args.advanced or args.all:
        print("Running advanced visualization suite...")
        viz = AdvancedVisualization()
        success &= viz.run_complete_advanced_visualization()
    
    if args.__dict__.get('3d', False) or args.all:
        print("Running 3D visualization suite...")
        viz = Visualization3D()
        success &= viz.run_complete_3d_analysis()
    
    if not any([args.professional, args.advanced, args.__dict__.get('3d', False), args.all]):
        print("No visualization specified. Use --help for options.")
        success = False
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
