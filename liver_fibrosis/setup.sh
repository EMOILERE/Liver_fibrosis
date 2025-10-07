#!/bin/bash
echo "Enhanced PhysiCell Liver Fibrosis Simulation"
echo "============================================="
echo ""

echo "Quick Options:"
echo "1. Generate ideal data demo (recommended)"
echo "2. Compile PhysiCell project"  
echo "3. Run analysis suite"
echo "4. Exit"
echo ""

read -p "Choose option (1-4): " opt

if [ "$opt" = "1" ]; then
    echo "Generating ideal data demo..."
    python3 demo_complete_analysis.py
elif [ "$opt" = "2" ]; then
    echo "Compiling PhysiCell..."
    make clean
    make
    if [ $? -eq 0 ]; then
        echo "SUCCESS: Compilation completed!"
        echo "Executable: liver_fibrosis_igem"
        echo "You can now run: ./liver_fibrosis_igem"
    else
        echo "ERROR: Compilation failed"
    fi
elif [ "$opt" = "3" ]; then
    echo "Running analysis..."
    python3 enhanced_analysis.py
    python3 professional_visualization.py
    python3 3D_visualization.py
elif [ "$opt" = "4" ]; then
    echo "Goodbye!"
    exit 0
else
    echo "Invalid option"
fi

echo "Done!"
