"""
Demo Module for Liver Fibrosis Simulation
Provides easy-to-use demonstration functionality
"""

from pathlib import Path
import sys
import time
from typing import Optional, Union

# Add parent directory to path to import demo modules
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

class LiverFibrosisDemo:
    """Main demo class for liver fibrosis simulation"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize the demo"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
        
        # Create output directory structure
        self.output_dir = self.base_dir / "demo_outputs"
        self.output_dir.mkdir(exist_ok=True)
        
        print("Liver Fibrosis Simulation Demo System")
        print("=" * 50)
    
    def check_dependencies(self) -> bool:
        """Check if required dependencies are available"""
        print("Checking dependencies...")
        
        required_packages = [
            'numpy', 'pandas', 'matplotlib', 'seaborn', 
            'scipy', 'plotly', 'scikit-learn'
        ]
        
        missing = []
        for package in required_packages:
            try:
                __import__(package)
                print(f"  OK: {package}")
            except ImportError:
                print(f"  MISSING: {package}")
                missing.append(package)
        
        # Special check for scikit-learn
        try:
            import sklearn
            print("  OK: scikit-learn")
        except ImportError:
            if 'scikit-learn' not in missing:
                print("  MISSING: scikit-learn")
                missing.append('scikit-learn')
        
        if missing:
            print(f"\nMissing packages: {', '.join(missing)}")
            print(f"Install with: pip install {' '.join(missing)}")
            return False
        
        print("All dependencies satisfied!")
        return True
    
    def run_data_generation_demo(self) -> bool:
        """Run data generation demonstration"""
        print("\nStep 1: Generate synthetic biological data...")
        print("-" * 40)
        
        try:
            from .data_generation import DataGenerator
            
            generator = DataGenerator(str(self.base_dir))
            success = generator.run_complete_data_generation()
            
            if success:
                print("SUCCESS: Data generation completed!")
                return True
            else:
                print("ERROR: Data generation failed")
                return False
                
        except Exception as e:
            print(f"ERROR: Data generation failed: {e}")
            return False
    
    def run_analysis_demo(self) -> bool:
        """Run analysis demonstration"""
        print("\nStep 2: Run comprehensive data analysis...")
        print("-" * 40)
        
        try:
            from .analysis import EnhancedAnalyzer
            
            analyzer = EnhancedAnalyzer(str(self.base_dir))
            success = analyzer.run_complete_analysis()
            
            if success:
                print("SUCCESS: Analysis completed!")
                return True
            else:
                print("ERROR: Analysis failed")
                return False
                
        except Exception as e:
            print(f"ERROR: Analysis failed: {e}")
            return False
    
    def run_visualization_demo(self) -> bool:
        """Run visualization demonstration"""
        print("\nStep 3: Generate professional visualizations...")
        print("-" * 40)
        
        try:
            from .visualization import ProfessionalVisualization, AdvancedVisualization
            
            # Professional visualization
            prof_viz = ProfessionalVisualization(str(self.base_dir))
            success1 = prof_viz.run_complete_suite()
            
            # Advanced visualization
            adv_viz = AdvancedVisualization(str(self.base_dir))
            success2 = adv_viz.run_complete_advanced_visualization()
            
            if success1 and success2:
                print("SUCCESS: All visualizations completed!")
                return True
            else:
                print("WARNING: Some visualizations may have failed")
                return success1 or success2
                
        except Exception as e:
            print(f"ERROR: Visualization failed: {e}")
            return False
    
    def run_complete_demo(self) -> bool:
        """Run complete demonstration workflow"""
        print("Starting Complete Liver Fibrosis Simulation Demo")
        print("=" * 60)
        
        start_time = time.time()
        
        # Check dependencies first
        if not self.check_dependencies():
            print("\nERROR: Please install missing dependencies first")
            return False
        
        print(f"\nStarting complete demo...")
        print(f"Estimated time: 3-5 minutes")
        
        # Run demo steps
        steps_completed = 0
        total_steps = 3
        
        # Step 1: Data Generation
        if self.run_data_generation_demo():
            steps_completed += 1
        
        # Step 2: Analysis
        if self.run_analysis_demo():
            steps_completed += 1
        
        # Step 3: Visualization
        if self.run_visualization_demo():
            steps_completed += 1
        
        # Summary
        duration = time.time() - start_time
        
        print("\n" + "=" * 60)
        print("Liver Fibrosis Simulation Demo Complete!")
        print("=" * 60)
        print(f"Total time: {duration:.1f} seconds")
        print(f"Completed steps: {steps_completed}/{total_steps}")
        
        if steps_completed >= 1:
            self._print_results_summary()
        
        success = steps_completed >= 2  # At least 2 out of 3 steps successful
        
        if success:
            print("\nDemo completed successfully!")
        else:
            print("\nDemo encountered issues. Check error messages above.")
        
        return success
    
    def _print_results_summary(self):
        """Print summary of generated results"""
        print("\nGenerated files:")
        
        # Check for data files
        base_dir = Path(self.base_dir)
        
        # Data files
        synthetic_dir = base_dir / "synthetic_results"
        if synthetic_dir.exists():
            print("   Data: synthetic_results/ - Complete biological dataset")
            print(f"      • Contains data for 9 experimental groups")
            print(f"      • 5 time points with detailed cell and microenvironment data per group")
        
        # Visualization files
        viz_dir = base_dir / "visualization_outputs"
        if viz_dir.exists():
            print("   Visualizations: visualization_outputs/")
            
            figures_dir = viz_dir / "figures"
            if figures_dir.exists():
                png_count = len(list(figures_dir.glob("*.png")))
                print(f"      • {png_count} high-resolution PNG figures")
            
            pdfs_dir = viz_dir / "pdfs"
            if pdfs_dir.exists():
                pdf_count = len(list(pdfs_dir.glob("*.pdf")))
                print(f"      • {pdf_count} publication-quality PDF files")
            
            animations_dir = viz_dir / "animations"
            if animations_dir.exists():
                gif_count = len(list(animations_dir.glob("*.gif")))
                print(f"      • {gif_count} animated GIF files")
            
            interactive_dir = viz_dir / "interactive"
            if interactive_dir.exists():
                html_count = len(list(interactive_dir.glob("*.html")))
                print(f"      • {html_count} interactive HTML visualizations")
        
        # Legacy visualization files
        viz_files = [
            ("morphological_analysis.png", "Morphological analysis"),
            ("molecular_heatmaps.png", "Molecular spatial distribution heatmaps"), 
            ("dose_response_analysis.png", "Dose-response analysis"),
            ("pathway_network_diagram.png", "Molecular pathway network diagram"),
            ("statistical_analysis.png", "Statistical analysis"),
            ("enhanced_time_series_analysis.png", "Enhanced time series analysis"),
            ("cellular_state_distribution.png", "Cell state distribution"),
            ("3D_spatial_analysis.png", "3D spatial pattern analysis")
        ]
        
        for filename, description in viz_files:
            if (base_dir / filename).exists():
                print(f"   Legacy: {filename} - {description}")
        
        # Interactive files
        interactive_files = [
            ("3D_cell_visualization.html", "Interactive 3D cell visualization"),
            ("3D_molecular_gradients.html", "Interactive molecular gradients")
        ]
        
        for filename, description in interactive_files:
            if (base_dir / filename).exists():
                print(f"   Interactive: {filename} - {description}")
        
        # Report files
        report_files = [
            ("enhanced_analysis_report.md", "Comprehensive analysis report"),
            ("3D_analysis_report.md", "3D analysis report")
        ]
        
        for filename, description in report_files:
            if (base_dir / filename).exists():
                print(f"   Report: {filename} - {description}")
        
        print("\nNext steps:")
        print("   1. Open *.html files to view interactive 3D visualizations")
        print("   2. View *.png and *.pdf files for publication-quality figures")
        print("   3. Read *.md analysis reports")
        print("   4. Explore synthetic_results/ directory for raw data")
        print("   5. Modify parameters and regenerate analysis")
        
        print("\nTips:")
        print("   • All data represents biologically realistic parameters")
        print("   • Visualizations showcase complete simulation capabilities")
        print("   • Methods can be applied to actual PhysiCell simulation data")

def main():
    """Command-line interface for demo"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Liver Fibrosis Simulation Demo')
    parser.add_argument('--complete', action='store_true', help='Run complete demo (default)')
    parser.add_argument('--data-only', action='store_true', help='Run data generation demo only')
    parser.add_argument('--analysis-only', action='store_true', help='Run analysis demo only')
    parser.add_argument('--visualization-only', action='store_true', help='Run visualization demo only')
    
    args = parser.parse_args()
    
    # Initialize demo
    demo = LiverFibrosisDemo()
    
    success = False
    
    if args.data_only:
        success = demo.run_data_generation_demo()
    elif args.analysis_only:
        success = demo.run_analysis_demo()
    elif args.visualization_only:
        success = demo.run_visualization_demo()
    else:
        # Default: run complete demo
        success = demo.run_complete_demo()
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
