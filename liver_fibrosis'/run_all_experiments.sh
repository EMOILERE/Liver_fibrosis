#!/bin/bash
# Linuxversion - runliver fibrosisexperimentgroup
# iGEM PhysiCell Simulation Project

echo "Starting liver fibrosis simulation suite on Linux..."
echo "=================================================="

# checkfile
if [ ! -f "config/simulation_parameters.xml" ]; then
 echo "Warning: config/simulation_parameters.xml not found, creating from template..."
 if [ -f "example_configs/basic_test.xml" ]; then
 cp example_configs/basic_test.xml config/simulation_parameters.xml
 echo "Created config/simulation_parameters.xml from basic_test template"
 else
 echo "Error: No template configuration found!"
 exit 1
 fi
fi

# verifyquick_config.pywork
echo "Testing quick_config.py..."
if ! python3 config/quick_config.py --config config/PhysiCell_settings.xml --list > /dev/null 2>&1; then
 echo "Warning: quick_config.py not working properly, using direct execution method"
 USE_DIRECT_METHOD=true
else
 USE_DIRECT_METHOD=false
fi

# checkcompile
if [ ! -f "liver_fibrosis_igem" ]; then
 echo "Executable not found, compiling project..."
 make clean
 make
 if [ $? -ne 0 ]; then
 echo "Compilation failed! Please check error messages above."
 exit 1
 fi
fi

# createoutputdirectory
mkdir -p output_control output_positive_model output_natural_exosomes output_NC_mimic
mkdir -p output_miR_455_3p output_miR_148a_5p output_dual_miRNA_1_1 output_dual_miRNA_1_2 output_dual_miRNA_2_1

# function：runexperimentgroup
run_experiment() {
 local exp_name=$1
 local exp_label=$2
 local output_dir=$3
 
 echo ""
 echo "${exp_label}..."
 
 if [ "$USE_DIRECT_METHOD" = true ]; then
 # method：modifyPhysiCell_settings.xmlexperimental_groupparameters
 cp config/PhysiCell_settings.xml config/temp_settings.xml
 sed -i "s/<string.*experimental_group.*>.*<\/string>/<string name=\"experimental_group\" units=\"dimensionless\">${exp_name}<\/string>/" config/temp_settings.xml
 ./liver_fibrosis_igem config/temp_settings.xml
 else
 # quick_config.pymethod
 if python3 config/quick_config.py --config config/PhysiCell_settings.xml --experiment ${exp_name} --output config/temp_config.xml; then
 ./liver_fibrosis_igem config/temp_config.xml
 else
 echo "Warning: quick_config.py failed, trying direct method for ${exp_name}"
 cp config/PhysiCell_settings.xml config/temp_settings.xml
 sed -i "s/<string.*experimental_group.*>.*<\/string>/<string name=\"experimental_group\" units=\"dimensionless\">${exp_name}<\/string>/" config/temp_settings.xml
 ./liver_fibrosis_igem config/temp_settings.xml
 fi
 fi
 
 # outputfile
 if [ -d "output" ]; then
 mv output/* ${output_dir}/ 2>/dev/null || true
 fi
 echo "${exp_label} completed."
}

# runexperimentgroup
run_experiment "control" "1. Running Control Group (No TGF-beta1, No Exosomes)" "output_control"
run_experiment "positive_model" "2. Running Positive Model Group (TGF-beta1 only)" "output_positive_model"
run_experiment "natural_exosomes" "3. Running Natural Exosomes Group" "output_natural_exosomes"
run_experiment "NC_mimic" "4. Running NC mimic Group" "output_NC_mimic"
run_experiment "miR_455_3p" "5. Running miR-455-3p Single Loading Group" "output_miR_455_3p"
run_experiment "miR_148a_5p" "6. Running miR-148a-5p Single Loading Group" "output_miR_148a_5p"
run_experiment "dual_miRNA_1_1" "7. Running Dual miRNA 1:1 Ratio Group" "output_dual_miRNA_1_1"
run_experiment "dual_miRNA_1_2" "8. Running Dual miRNA 1:2 Ratio Group" "output_dual_miRNA_1_2"
run_experiment "dual_miRNA_2_1" "9. Running Dual miRNA 2:1 Ratio Group" "output_dual_miRNA_2_1"

# file
rm -f config/temp_config.xml config/temp_settings.xml

echo ""
echo "=================================================="
echo "All experimental groups completed!"
echo ""
echo "Output directories created:"
echo "- output_control"
echo "- output_positive_model" 
echo "- output_natural_exosomes"
echo "- output_NC_mimic"
echo "- output_miR_455_3p"
echo "- output_miR_148a_5p"
echo "- output_dual_miRNA_1_1"
echo "- output_dual_miRNA_1_2"
echo "- output_dual_miRNA_2_1"
echo ""
echo "Next steps:"
echo "1. Review SVG files for visual results"
echo "2. Analyze MultiCellDS data files"
echo "3. Compare activation rates and fibrosis markers"
echo "4. Look for synergistic effects in dual miRNA groups"
echo ""
echo "Run: python3 analyze_results.py"
