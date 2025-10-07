/*
###############################################################################
# Parameter Parser Implementation #
# Utility for loading simulation parameters from XML files #
###############################################################################
*/

#include "parameter_parser.h"
#include <iostream>
#include <fstream>

bool load_simulation_parameters(const std::string& filename) {
    // Simple file existence check
    std::ifstream file(filename);
    if (file.good()) {
        file.close();
        return true;
    }
    return false;
}