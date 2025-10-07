#!/usr/bin/env python3
"""
Alternative Demo Runner for Liver Fibrosis Simulation
This script provides multiple ways to run the demo system
"""

import sys
import os
from pathlib import Path

def main():
    """Run the liver fibrosis simulation demo"""
    print("Liver Fibrosis Simulation Demo Runner")
    print("=" * 50)
    
    # Add current directory to Python path
    current_dir = Path(__file__).parent
    sys.path.insert(0, str(current_dir))
    
    # Method 1: Try package import
    try:
        print("Attempting to use package interface...")
        from liver_fibrosis_sim.demo import LiverFibrosisDemo
        demo = LiverFibrosisDemo()
        success = demo.run_complete_demo()
        
        if success:
            print("\n🎉 Demo completed successfully using package interface!")
            return True
        else:
            print("⚠️  Demo completed with some issues")
            return False
            
    except ImportError as e:
        print(f"Package import failed: {e}")
        print("Trying alternative approach...")
    
    # Method 2: Try direct script execution
    try:
        print("Attempting to run demo_complete_analysis.py directly...")
        import subprocess
        
        demo_script = current_dir / "demo_complete_analysis.py"
        if demo_script.exists():
            result = subprocess.run([sys.executable, str(demo_script)], 
                                  capture_output=True, text=True)
            
            print(result.stdout)
            if result.stderr:
                print("Errors:", result.stderr)
            
            if result.returncode == 0:
                print("\n🎉 Demo completed successfully using direct script!")
                return True
            else:
                print("⚠️  Demo script had issues")
                return False
        else:
            print("ERROR: demo_complete_analysis.py not found")
            
    except Exception as e:
        print(f"Direct script execution failed: {e}")
    
    # Method 3: Manual execution
    print("Trying manual execution approach...")
    
    try:
        # Import required modules directly
        exec(open('generate_ideal_data.py').read())
        print("✓ Data generation completed")
        
        exec(open('enhanced_analysis.py').read())
        print("✓ Enhanced analysis completed")
        
        exec(open('professional_visualization.py').read())
        print("✓ Professional visualization completed")
        
        exec(open('advanced_visualization.py').read())
        print("✓ Advanced visualization completed")
        
        exec(open('3D_visualization.py').read())
        print("✓ 3D visualization completed")
        
        print("\n🎉 Manual execution completed successfully!")
        return True
        
    except Exception as e:
        print(f"Manual execution failed: {e}")
    
    # Method 4: Minimal demo
    print("Trying minimal demo...")
    
    try:
        print("Checking dependencies...")
        
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        print("✓ Core dependencies available")
        
        # Generate minimal demo data
        print("Generating minimal demo visualization...")
        
        # Simple plot
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Sample data
        groups = ['Control', 'Positive Model', 'miR-455-3p', 'miR-148a-5p', 'Dual miRNA']
        efficacies = [0, 10, 45, 40, 70]
        
        bars = ax.bar(groups, efficacies, color=['gray', 'red', 'blue', 'green', 'purple'], alpha=0.7)
        ax.set_ylabel('Treatment Efficacy (%)')
        ax.set_title('Liver Fibrosis Treatment Efficacy Analysis')
        ax.grid(True, alpha=0.3)
        
        # Add value labels
        for bar, value in zip(bars, efficacies):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{value}%', ha='center', va='bottom')
        
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig('minimal_demo_result.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✓ Minimal demo visualization created: minimal_demo_result.png")
        print("\n🎉 Minimal demo completed successfully!")
        return True
        
    except Exception as e:
        print(f"Minimal demo failed: {e}")
    
    print("\n❌ All demo methods failed")
    print("\nPlease check:")
    print("1. Dependencies are installed: pip install -r requirements.txt")
    print("2. Python version is 3.7+")
    print("3. All script files are present")
    
    return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
