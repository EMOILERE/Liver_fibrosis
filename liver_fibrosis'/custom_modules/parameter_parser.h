/*
###############################################################################
# Parameter Parser for iGEM Liver Fibrosis Simulation #
# Utility for loading simulation parameters from XML files #
###############################################################################
*/

#ifndef PARAMETER_PARSER_H
#define PARAMETER_PARSER_H

#include <string>
#include <map>

bool load_simulation_parameters(const std::string& filename);

#endif // PARAMETER_PARSER_H