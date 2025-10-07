#!/usr/bin/env python3
"""
Package Installation Script for Liver Fibrosis Simulation
Installs the package in development mode and fixes common import issues
"""

import sys
import os
import subprocess
from pathlib import Path

def install_package():
    """Install the liver_fibrosis_sim package"""
    print("Installing Liver Fibrosis Simulation Package...")
    print("=" * 50)
    
    # Get current directory
    current_dir = Path(__file__).parent
    print(f"Current directory: {current_dir}")
    
    # Check if setup.py exists
    setup_file = current_dir / "setup.py"
    if not setup_file.exists():
        print(f"ERROR: setup.py not found in {current_dir}")
        return False
    
    # Check if liver_fibrosis_sim directory exists
    package_dir = current_dir / "liver_fibrosis_sim"
    if not package_dir.exists():
        print(f"ERROR: Package directory {package_dir} not found")
        return False
    
    try:
        # Install in development mode
        print("Installing package in development mode...")
        result = subprocess.run([
            sys.executable, "-m", "pip", "install", "-e", "."
        ], cwd=current_dir, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("SUCCESS: Package installed successfully!")
            print("\nYou can now use:")
            print("  python -m liver_fibrosis_sim.demo")
            print("  liver-fibrosis-demo")
            print("  liver-fibrosis-analyze")
            print("  liver-fibrosis-visualize")
            return True
        else:
            print(f"ERROR: Package installation failed")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"ERROR: Installation failed with exception: {e}")
        return False

def test_import():
    """Test if package can be imported"""
    print("\nTesting package import...")
    
    try:
        import liver_fibrosis_sim
        print("SUCCESS: Package imported successfully!")
        
        # Test main classes
        from liver_fibrosis_sim import LiverFibrosisSimulation
        print("SUCCESS: LiverFibrosisSimulation imported!")
        
        from liver_fibrosis_sim.demo import LiverFibrosisDemo
        print("SUCCESS: LiverFibrosisDemo imported!")
        
        return True
        
    except ImportError as e:
        print(f"ERROR: Import failed: {e}")
        print("\nTrying alternative approach...")
        
        # Add current directory to Python path
        current_dir = Path(__file__).parent
        sys.path.insert(0, str(current_dir))
        
        try:
            import liver_fibrosis_sim
            print("SUCCESS: Package imported with sys.path modification!")
            return True
        except ImportError as e2:
            print(f"ERROR: Alternative import also failed: {e2}")
            return False

def run_quick_demo():
    """Run a quick demo to test functionality"""
    print("\nRunning quick demo test...")
    
    try:
        from liver_fibrosis_sim.demo import LiverFibrosisDemo
        demo = LiverFibrosisDemo()
        
        # Just test dependency check
        deps_ok = demo.check_dependencies()
        
        if deps_ok:
            print("SUCCESS: Demo system working correctly!")
            return True
        else:
            print("WARNING: Some dependencies missing, but demo system is accessible")
            return True
            
    except Exception as e:
        print(f"ERROR: Demo test failed: {e}")
        return False

def create_alternative_runner():
    """Create alternative runner script if package installation fails"""
    print("\nCreating alternative runner script...")
    
    runner_script = """#!/usr/bin/env python3
import sys
import os
from pathlib import Path

# Add current directory to Python path
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

# Import and run demo
try:
    from liver_fibrosis_sim.demo import LiverFibrosisDemo
    demo = LiverFibrosisDemo()
    demo.run_complete_demo()
except ImportError as e:
    print(f"Import error: {e}")
    print("Falling back to direct script execution...")
    
    # Try to run demo_complete_analysis.py directly
    demo_script = current_dir / "demo_complete_analysis.py"
    if demo_script.exists():
        import subprocess
        subprocess.run([sys.executable, str(demo_script)])
    else:
        print("ERROR: Cannot find demo script")
"""
    
    runner_file = Path(__file__).parent / "run_demo.py"
    with open(runner_file, 'w') as f:
        f.write(runner_script)
    
    # Make executable on Unix systems
    try:
        import stat
        runner_file.chmod(runner_file.stat().st_mode | stat.S_IEXEC)
    except:
        pass
    
    print(f"Created alternative runner: {runner_file}")
    print("You can run: python run_demo.py")

def main():
    """Main installation and testing workflow"""
    print("Liver Fibrosis Simulation Package Installer")
    print("=" * 60)
    
    success = True
    
    # Step 1: Install package
    if not install_package():
        print("Package installation failed, continuing with tests...")
        success = False
    
    # Step 2: Test import
    if not test_import():
        print("Import test failed, creating alternative runner...")
        success = False
    
    # Step 3: Test demo
    if success and not run_quick_demo():
        print("Demo test failed")
        success = False
    
    # Step 4: Create alternative if needed
    if not success:
        create_alternative_runner()
    
    print("\n" + "=" * 60)
    if success:
        print("🎉 Installation completed successfully!")
        print("\nRecommended usage:")
        print("  python -m liver_fibrosis_sim.demo")
        print("  # or")
        print("  liver-fibrosis-demo")
    else:
        print("⚠️  Package installation had issues")
        print("\nAlternative usage:")
        print("  python run_demo.py")
        print("  # or")
        print("  python demo_complete_analysis.py")
    
    print("\nFor advanced usage, see README.md")
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
