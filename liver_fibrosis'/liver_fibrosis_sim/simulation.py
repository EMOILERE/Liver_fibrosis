"""
Liver Fibrosis Simulation Module
Main interface for PhysiCell-based liver fibrosis simulations
"""

import os
import subprocess
import sys
from pathlib import Path
import time
from typing import Dict, List, Optional, Union

class LiverFibrosisSimulation:
    """
    Main simulation class for liver fibrosis modeling
    
    This class provides a high-level interface to the PhysiCell-based
    liver fibrosis simulation, including compilation, execution, and
    basic result management.
    """
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """
        Initialize the simulation
        
        Args:
            base_dir: Base directory for simulation files
        """
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
        
        # Define important paths
        self.config_dir = self.base_dir / "config"
        self.custom_modules_dir = self.base_dir / "custom_modules" 
        self.output_dir = self.base_dir / "output"
        self.makefile = self.base_dir / "Makefile"
        self.executable = self.base_dir / "liver_fibrosis_igem"
        
        # Experimental groups
        self.experimental_groups = [
            'control', 'positive_model', 'natural_exosomes', 'NC_mimic',
            'miR_455_3p', 'miR_148a_5p', 'dual_miRNA_1_1', 'dual_miRNA_1_2', 'dual_miRNA_2_1'
        ]
        
        # Simulation parameters
        self.default_params = {
            'duration': 2880,  # 48 hours in minutes
            'output_interval': 60,  # Every hour
            'enable_3D': False,
            'reduce_noise': True
        }
    
    def check_dependencies(self) -> bool:
        """
        Check if all required dependencies are available
        
        Returns:
            bool: True if all dependencies are satisfied
        """
        # Check for compiler
        try:
            result = subprocess.run(['g++', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                print("ERROR: g++ compiler not found")
                return False
        except FileNotFoundError:
            print("ERROR: g++ compiler not found")
            return False
        
        # Check for make
        try:
            result = subprocess.run(['make', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                print("ERROR: make not found")
                return False
        except FileNotFoundError:
            print("ERROR: make not found")
            return False
        
        # Check for Python packages
        required_packages = ['numpy', 'pandas', 'matplotlib', 'seaborn', 'scipy', 'plotly']
        missing_packages = []
        
        for package in required_packages:
            try:
                __import__(package)
            except ImportError:
                missing_packages.append(package)
        
        if missing_packages:
            print(f"ERROR: Missing Python packages: {', '.join(missing_packages)}")
            print(f"Install with: pip install {' '.join(missing_packages)}")
            return False
        
        print("All dependencies satisfied")
        return True
    
    def compile(self, clean_first: bool = True) -> bool:
        """
        Compile the PhysiCell simulation
        
        Args:
            clean_first: Whether to clean previous build first
            
        Returns:
            bool: True if compilation successful
        """
        print("Compiling PhysiCell simulation...")
        
        try:
            # Change to simulation directory
            original_dir = os.getcwd()
            os.chdir(self.base_dir)
            
            if clean_first:
                result = subprocess.run(['make', 'clean'], 
                                      capture_output=True, text=True)
                if result.returncode == 0:
                    print("Previous build cleaned")
            
            # Compile
            result = subprocess.run(['make', '-j4'], 
                                  capture_output=True, text=True)
            
            if result.returncode == 0:
                print("Compilation successful")
                return True
            else:
                print(f"Compilation failed: {result.stderr}")
                return False
                
        except Exception as e:
            print(f"Compilation error: {e}")
            return False
        finally:
            os.chdir(original_dir)
    
    def configure_experiment(self, group: str, **kwargs) -> bool:
        """
        Configure simulation for specific experimental group
        
        Args:
            group: Experimental group name
            **kwargs: Additional parameters
            
        Returns:
            bool: True if configuration successful
        """
        if group not in self.experimental_groups:
            print(f"ERROR: Unknown experimental group: {group}")
            print(f"Available groups: {', '.join(self.experimental_groups)}")
            return False
        
        print(f"Configuring experiment: {group}")
        
        # Use quick_config.py to set parameters
        config_script = self.base_dir / "config" / "quick_config.py"
        if not config_script.exists():
            print("ERROR: Configuration script not found")
            return False
        
        try:
            # Build command
            cmd = [sys.executable, str(config_script), '--experiment', group]
            
            # Add additional parameters
            if kwargs.get('enable_3D', False):
                cmd.extend(['--3d'])
            if kwargs.get('duration'):
                cmd.extend(['--duration', str(kwargs['duration'])])
            if kwargs.get('reduce_noise', True):
                cmd.extend(['--reduce-noise'])
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"Configuration successful for {group}")
                return True
            else:
                print(f"Configuration failed: {result.stderr}")
                return False
                
        except Exception as e:
            print(f"Configuration error: {e}")
            return False
    
    def run_simulation(self, timeout_minutes: int = 60) -> bool:
        """
        Run the compiled simulation
        
        Args:
            timeout_minutes: Maximum runtime in minutes
            
        Returns:
            bool: True if simulation completed successfully
        """
        if not self.executable.exists():
            print("ERROR: Executable not found. Please compile first.")
            return False
        
        print("Running simulation...")
        start_time = time.time()
        
        try:
            # Change to simulation directory
            original_dir = os.getcwd()
            os.chdir(self.base_dir)
            
            # Run simulation
            result = subprocess.run([str(self.executable)], 
                                  timeout=timeout_minutes * 60,
                                  capture_output=True, text=True)
            
            duration = time.time() - start_time
            
            if result.returncode == 0:
                print(f"Simulation completed successfully in {duration:.1f} seconds")
                return True
            else:
                print(f"Simulation failed: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            print(f"Simulation timed out after {timeout_minutes} minutes")
            return False
        except Exception as e:
            print(f"Simulation error: {e}")
            return False
        finally:
            os.chdir(original_dir)
    
    def run_single_experiment(self, group: str, **kwargs) -> bool:
        """
        Run complete simulation for single experimental group
        
        Args:
            group: Experimental group name
            **kwargs: Simulation parameters
            
        Returns:
            bool: True if successful
        """
        print(f"Running complete experiment: {group}")
        
        # Configure
        if not self.configure_experiment(group, **kwargs):
            return False
        
        # Compile if needed
        if not self.executable.exists():
            if not self.compile():
                return False
        
        # Run simulation
        return self.run_simulation(kwargs.get('timeout', 60))
    
    def run_all_experiments(self, **kwargs) -> Dict[str, bool]:
        """
        Run all experimental groups
        
        Args:
            **kwargs: Simulation parameters
            
        Returns:
            dict: Results for each group
        """
        print("Running all experimental groups...")
        results = {}
        
        # Compile once
        if not self.compile():
            print("ERROR: Compilation failed")
            return results
        
        for group in self.experimental_groups:
            print(f"\nProcessing group: {group}")
            results[group] = self.run_single_experiment(group, **kwargs)
        
        # Summary
        successful = sum(results.values())
        total = len(results)
        print(f"\nExperiment suite completed: {successful}/{total} successful")
        
        return results
    
    def run_demo(self) -> bool:
        """
        Run demonstration with synthetic data generation
        
        Returns:
            bool: True if demo completed successfully
        """
        print("Running liver fibrosis simulation demo...")
        
        try:
            # Import and run data generation
            from .data_generation import DataGenerator
            
            generator = DataGenerator(str(self.base_dir))
            if not generator.run_complete_data_generation():
                print("ERROR: Demo data generation failed")
                return False
            
            # Run visualization
            from .visualization import AdvancedVisualization
            
            visualizer = AdvancedVisualization(str(self.base_dir))
            if not visualizer.run_complete_advanced_visualization():
                print("ERROR: Demo visualization failed")
                return False
            
            print("Demo completed successfully!")
            return True
            
        except ImportError as e:
            print(f"ERROR: Missing dependencies for demo: {e}")
            return False
        except Exception as e:
            print(f"ERROR: Demo failed: {e}")
            return False
    
    def get_results_summary(self) -> Dict:
        """
        Get summary of simulation results
        
        Returns:
            dict: Summary statistics
        """
        summary = {
            'base_directory': str(self.base_dir),
            'executable_exists': self.executable.exists(),
            'config_files': list(self.config_dir.glob('*.xml')) if self.config_dir.exists() else [],
            'output_directories': list(self.base_dir.glob('output_*')),
            'experimental_groups': self.experimental_groups
        }
        
        return summary

def main():
    """Command-line interface for simulation"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Liver Fibrosis Simulation')
    parser.add_argument('--demo', action='store_true', help='Run demonstration')
    parser.add_argument('--group', help='Run specific experimental group')
    parser.add_argument('--all', action='store_true', help='Run all experimental groups')
    parser.add_argument('--3d', action='store_true', help='Enable 3D simulation')
    parser.add_argument('--duration', type=int, default=2880, help='Simulation duration in minutes')
    
    args = parser.parse_args()
    
    # Initialize simulation
    sim = LiverFibrosisSimulation()
    
    # Check dependencies
    if not sim.check_dependencies():
        sys.exit(1)
    
    # Run requested operation
    if args.demo:
        success = sim.run_demo()
    elif args.group:
        success = sim.run_single_experiment(args.group, enable_3D=args.__dict__.get('3d', False), duration=args.duration)
    elif args.all:
        results = sim.run_all_experiments(enable_3D=args.__dict__.get('3d', False), duration=args.duration)
        success = any(results.values())
    else:
        print("No operation specified. Use --help for options.")
        success = False
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
