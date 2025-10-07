#!/usr/bin/env python3
"""
Complete Demo Script - Generate ideal data and run all visualization analysis
Showcase the complete functionality of PhysiCell Enhanced Liver Fibrosis Simulation
"""

import os
import sys
import subprocess
from pathlib import Path
import time

def print_banner():
    """Print welcome banner"""
    print("=" * 80)
    print("PhysiCell Enhanced Liver Fibrosis Simulation - Complete Demo System")
    print("=" * 80)
    print("This demo will showcase:")
    print("• 33 cellular state variables for complex biological modeling")
    print("• 12 microenvironmental factors with 3D spatial distribution")  
    print("• Professional publication-quality visualization")
    print("• Comprehensive statistical and time series analysis")
    print("• Realistic 3D culture environment simulation")
    print("• miRNA synergistic therapeutic effect analysis")
    print("=" * 80)

def check_requirements():
    """Check required Python packages"""
    required_packages = [
        'numpy', 'pandas', 'matplotlib', 'seaborn', 
        'scipy', 'plotly'
    ]
    
    # Try importing scikit-learn with different module names
    sklearn_names = ['sklearn', 'scikit-learn', 'scikit_learn']
    
    print("Checking dependencies...")
    missing = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"  OK: {package}")
        except ImportError:
            print(f"  MISSING: {package}")
            missing.append(package)
    
    # Special check for scikit-learn
    sklearn_found = False
    for name in sklearn_names:
        try:
            __import__(name)
            sklearn_found = True
            print(f"  OK: scikit-learn")
            break
        except ImportError:
            continue
    
    if not sklearn_found:
        try:
            import sklearn
            sklearn_found = True
            print(f"  OK: scikit-learn")
        except ImportError:
            print(f"  MISSING: scikit-learn")
            missing.append('scikit-learn')
    
    if missing:
        print(f"\nMissing packages: {', '.join(missing)}")
        print(f"Install with: pip install {' '.join(missing)}")
        return False
    
    print("All dependencies satisfied!")
    return True

def run_step(step_name, script_name, description):
    """Run single analysis step"""
    print(f"\n{step_name}: {description}")
    print("-" * 50)
    
    try:
        # Capture output but still show important information
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=False,
            text=True,
            check=True
        )
        print(f"SUCCESS: {step_name} completed!")
        return True
    
    except subprocess.CalledProcessError as e:
        print(f"WARNING: {step_name} encountered issues: {e}")
        print("Continuing to next step...")
        return False
    
    except FileNotFoundError:
        print(f"ERROR: Script not found: {script_name}")
        return False

def main():
    """Main demo function"""
    print_banner()
    
    # Check dependencies
    if not check_requirements():
        print("\nERROR: Please install missing dependencies first")
        return False
    
    print(f"\nStarting complete demo...")
    print(f"Estimated time: 3-5 minutes")
    
    start_time = time.time()
    
    # Execute steps
    steps = [
        ("Step 1", "generate_ideal_data.py", "Generate ideal biological data"),
        ("Step 2", "professional_visualization.py", "Create professional publication-quality visualization"),
        ("Step 3", "enhanced_analysis.py", "Run enhanced data analysis"),
        ("Step 4", "3D_visualization.py", "Generate interactive 3D visualization")
    ]
    
    completed_steps = 0
    
    for step_name, script_name, description in steps:
        success = run_step(step_name, script_name, description)
        if success:
            completed_steps += 1
        time.sleep(0.5)  # Brief pause
    
    # Summary
    duration = time.time() - start_time
    
    print("\n" + "=" * 80)
    print("PhysiCell Enhanced Simulation Demo Complete!")
    print("=" * 80)
    print(f"Total time: {duration:.1f} seconds")
    print(f"Completed steps: {completed_steps}/{len(steps)}")
    
    if completed_steps >= 1:
        print("\nGenerated files:")
        
        # Check generated files
        base_dir = Path(".")
        
        # Data files
        synthetic_dir = base_dir / "synthetic_results"
        if synthetic_dir.exists():
            print("   Data: synthetic_results/ - Ideal dataset")
            print(f"      • Contains complete data for 9 experimental groups")
            print(f"      • 5 time points with detailed cell and microenvironment data per group")
        
        # Visualization files
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
                print(f"   Chart: {filename} - {description}")
        
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
            ("3D_analysis_report.md", "3D analysis report"),
            ("synthetic_results/README_SYNTHETIC_DATA.md", "Data usage instructions")
        ]
        
        for filename, description in report_files:
            if (base_dir / filename).exists():
                print(f"   Report: {filename} - {description}")
    
    print("\nNext steps:")
    print("   1. Open *.html files to view interactive 3D visualizations")
    print("   2. View *.png chart files")
    print("   3. Read *.md analysis reports")
    print("   4. Explore synthetic_results/ directory for data")
    print("   5. Modify parameters and regenerate data/analysis")
    
    print("\nTips:")
    print("   • All data is ideal data based on real biological parameters")
    print("   • Visualizations showcase complete PhysiCell simulation capabilities")
    print("   • These methods can be applied to actual PhysiCell simulation data")
    
    print("\nMore information:")
    print("   • README_3D_ENHANCED.md - Complete feature description")
    print("   • ENHANCED_FEATURES.md - Technical specifications")
    print("   • USAGE_GUIDE.md - Detailed usage guide")
    
    print("\nCongratulations! You have experienced a world-class 3D biological simulation platform!")
    
    return completed_steps > 0

if __name__ == "__main__":
    try:
        success = main()
        if not success:
            print("\nERROR: Demo could not complete successfully, please check error messages")
            sys.exit(1)
    except KeyboardInterrupt:
        print("\nINTERRUPTED: User interrupted the demo")
        sys.exit(0)
    except Exception as e:
        print(f"\nERROR: Unexpected error during demo: {e}")
        sys.exit(1)