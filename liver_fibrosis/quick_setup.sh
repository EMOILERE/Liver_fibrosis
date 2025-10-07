#!/bin/bash

echo "================================================"
echo " Enhanced PhysiCell Liver Fibrosis Simulation"
echo " Quick Setup and Management Tool"
echo "================================================"
echo ""

if [ ! -f "config/PhysiCell_settings.xml" ]; then
    echo "WARNING: Configuration file config/PhysiCell_settings.xml not found"
    if [ -f "config/PhysiCell_settings.xml.backup" ]; then
        echo "Restoring configuration file from backup..."
        cp config/PhysiCell_settings.xml.backup config/PhysiCell_settings.xml
        echo "Configuration file restored successfully"
    else
        echo "ERROR: No configuration file or backup found"
        echo "Please ensure PhysiCell_settings.xml exists in config/ directory"
        exit 1
    fi
fi

while true; do
    echo ""
    echo "=== MAIN MENU ==="
    echo "1. Compile project"
    echo "2. Run single experimental group"
    echo "3. Run all experimental groups"
    echo "4. Generate ideal data demo"
    echo "5. Run professional analysis"
    echo "6. Modify simulation parameters"
    echo "7. Clean output files"
    echo "8. View analysis results"
    echo "9. Exit"
    echo ""
    
    read -p "Please select option (1-9): " choice
    
    case "$choice" in
        1)
            echo ""
            echo "=== COMPILING PROJECT ==="
            make clean
            make -j4
            if [ $? -eq 0 ]; then
                echo "SUCCESS: Compilation completed successfully!"
            else
                echo "ERROR: Compilation failed. Please check error messages above."
            fi
            read -p "Press Enter to continue..."
            ;;
        2)
            echo ""
            echo "=== AVAILABLE EXPERIMENTAL GROUPS ==="
            echo "1. control - Control group (low TGF-beta1)"
            echo "2. positive_model - Positive model (high TGF-beta1)"
            echo "3. natural_exosomes - Natural exosomes treatment"
            echo "4. NC_mimic - Negative control mimic"
            echo "5. miR_455_3p - miR-455-3p treatment"
            echo "6. miR_148a_5p - miR-148a-5p treatment"
            echo "7. dual_miRNA_1_1 - Dual miRNA ratio 1:1"
            echo "8. dual_miRNA_1_2 - Dual miRNA ratio 1:2"
            echo "9. dual_miRNA_2_1 - Dual miRNA ratio 2:1"
            echo ""
            
            read -p "Select experimental group (1-9): " exp_choice
            
            case "$exp_choice" in
                1) exp_name="control" ;;
                2) exp_name="positive_model" ;;
                3) exp_name="natural_exosomes" ;;
                4) exp_name="NC_mimic" ;;
                5) exp_name="miR_455_3p" ;;
                6) exp_name="miR_148a_5p" ;;
                7) exp_name="dual_miRNA_1_1" ;;
                8) exp_name="dual_miRNA_1_2" ;;
                9) exp_name="dual_miRNA_2_1" ;;
                *) echo "Invalid choice. Please try again."; continue ;;
            esac
            
            echo ""
            echo "=== SIMULATION DIMENSIONS ==="
            echo "1. 2D simulation (faster, classic)"
            echo "2. 3D simulation (comprehensive, realistic)"
            echo ""
            read -p "Select simulation type (1-2): " sim_type
            
            dimension_flag=""
            if [ "$sim_type" = "2" ]; then
                dimension_flag="--3d"
                echo "Selected: 3D simulation"
            else
                echo "Selected: 2D simulation"
            fi
            
            echo "Setting up experimental group: $exp_name"
            python3 config/quick_config.py --experiment "$exp_name" $dimension_flag
            if [ $? -eq 0 ]; then
                echo "Running simulation..."
                if [ -f "./project" ]; then
                    ./project
                elif [ -f "./project.exe" ]; then
                    ./project.exe
                else
                    echo "ERROR: Executable not found. Please compile first."
                    read -p "Press Enter to continue..."
                    continue
                fi
                echo "SUCCESS: Simulation completed!"
            else
                echo "ERROR: Configuration failed. Please check Python script."
            fi
            read -p "Press Enter to continue..."
            ;;
        3)
            echo ""
            echo "=== RUNNING ALL EXPERIMENTAL GROUPS ==="
            echo "This will run all 9 experimental groups sequentially."
            echo "Estimated time: 2-6 hours depending on settings."
            echo ""
            
            echo "=== SIMULATION DIMENSIONS ==="
            echo "1. 2D simulation (faster)"
            echo "2. 3D simulation (comprehensive)"
            echo ""
            read -p "Select simulation type (1-2): " sim_type
            
            dimension_flag=""
            if [ "$sim_type" = "2" ]; then
                dimension_flag="--3d"
                echo "Running all groups in 3D mode..."
            else
                echo "Running all groups in 2D mode..."
            fi
            
            if [ -f "run_complete_simulation.py" ]; then
                python3 run_complete_simulation.py $dimension_flag
            else
                echo "ERROR: run_complete_simulation.py not found"
            fi
            echo ""
            echo "SUCCESS: All experimental groups completed!"
            read -p "Press Enter to continue..."
            ;;
        4)
            echo ""
            echo "=== GENERATING IDEAL DATA DEMO ==="
            echo "This generates ideal synthetic data for testing and visualization."
            echo "Processing time: 2-5 minutes"
            echo ""
            
            if [ -f "demo_complete_analysis.py" ]; then
                python3 demo_complete_analysis.py
                echo "SUCCESS: Ideal data demo completed!"
            else
                echo "ERROR: demo_complete_analysis.py not found"
            fi
            read -p "Press Enter to continue..."
            ;;
        5)
            echo ""
            echo "=== PROFESSIONAL ANALYSIS SUITE ==="
            echo "Running comprehensive analysis and visualization..."
            echo ""
            
            if [ -f "enhanced_analysis.py" ]; then
                echo "Running enhanced_analysis.py..."
                python3 enhanced_analysis.py
            else
                echo "WARNING: enhanced_analysis.py not found, skipping..."
            fi
            
            if [ -f "professional_visualization.py" ]; then
                echo "Running professional_visualization.py..."
                python3 professional_visualization.py
            else
                echo "WARNING: professional_visualization.py not found, skipping..."
            fi
            
            if [ -f "3D_visualization.py" ]; then
                echo "Running 3D_visualization.py..."
                python3 3D_visualization.py
            else
                echo "WARNING: 3D_visualization.py not found, skipping..."
            fi
            
            echo "SUCCESS: Professional analysis completed!"
            read -p "Press Enter to continue..."
            ;;
        6)
            echo ""
            echo "=== PARAMETER MODIFICATION ==="
            echo "1. Modify synergy factor (default: 2.5)"
            echo "2. Modify TGF-beta1 concentration"
            echo "3. Modify simulation duration"
            echo "4. Reset to default parameters"
            echo ""
            read -p "Select option (1-4): " param_choice
            
            case "$param_choice" in
                1)
                    echo "Current synergy factor range: 1.0-4.0 (optimal: 2.5)"
                    read -p "Enter new synergy factor: " synergy_value
                    if [ ! -z "$synergy_value" ]; then
                        echo "Setting synergy factor to: $synergy_value"
                        python3 config/quick_config.py --synergy "$synergy_value"
                    fi
                    ;;
                2)
                    echo "TGF-beta1 concentration settings:"
                    echo "Control: 2.0 ng/mL"
                    echo "Treatment: 15.0 ng/mL"
                    echo "Range: 5-25 ng/mL"
                    read -p "Enter TGF-beta1 concentration (ng/mL): " tgf_value
                    if [ ! -z "$tgf_value" ]; then
                        echo "Setting TGF-beta1 concentration to: $tgf_value ng/mL"
                        python3 config/quick_config.py --tgf "$tgf_value"
                    fi
                    ;;
                3)
                    echo "Simulation duration settings:"
                    echo "Short: 1440 minutes (24 hours)"
                    echo "Medium: 2880 minutes (48 hours)"
                    echo "Long: 4320 minutes (72 hours)"
                    echo "Extended: 7200 minutes (120 hours)"
                    read -p "Enter duration in minutes: " duration_value
                    if [ ! -z "$duration_value" ]; then
                        echo "Setting simulation duration to: $duration_value minutes"
                        python3 config/quick_config.py --duration "$duration_value"
                    fi
                    ;;
                4)
                    echo "Resetting to default parameters..."
                    python3 config/quick_config.py --reset-defaults
                    echo "Parameters reset to defaults"
                    ;;
            esac
            read -p "Press Enter to continue..."
            ;;
        7)
            echo ""
            echo "=== CLEANING OUTPUT FILES ==="
            echo "This will remove all simulation output and analysis files."
            read -p "Are you sure? (y/N): " confirm
            if [ "$confirm" = "y" ] || [ "$confirm" = "Y" ]; then
                rm -rf output output_* synthetic_results *.svg *.png *.html *.txt *.md
                mkdir -p output
                echo "SUCCESS: All output files cleaned!"
            else
                echo "Operation cancelled."
            fi
            read -p "Press Enter to continue..."
            ;;
        8)
            echo ""
            echo "=== ANALYSIS RESULTS ==="
            
            found_results=false
            
            if ls *.png 1> /dev/null 2>&1; then
                echo "Found: PNG chart files"
                found_results=true
            fi
            
            if ls *.html 1> /dev/null 2>&1; then
                echo "Found: HTML interactive visualizations"
                found_results=true
            fi
            
            if ls *_report.md 1> /dev/null 2>&1; then
                echo "Found: Analysis reports"
                found_results=true
            fi
            
            if [ "$found_results" = false ]; then
                echo "No analysis results found."
                echo "Please run simulations or generate demo data first."
            else
                echo ""
                echo "Analysis results are available in the current directory."
                echo "Open .html files in a web browser for interactive visualization."
                echo "View .png files for static charts."
                echo "Read .md files for detailed reports."
            fi
            
            read -p "Press Enter to continue..."
            ;;
        9)
            echo ""
            echo "================================================"
            echo " Thank you for using Enhanced PhysiCell!"
            echo " Your liver fibrosis simulation is ready!"
            echo "================================================"
            echo ""
            exit 0
            ;;
        *)
            echo "Invalid choice. Please select 1-9."
            ;;
    esac
done