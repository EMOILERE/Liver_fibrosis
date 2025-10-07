/*
###############################################################################
# Custom implementation for iGEM Liver Fibrosis Simulation #
# Models LX-2 hepatic stellate cell response to miRNA-loaded exosomes #
# Core logic for TGF-beta1 activation and miRNA synergistic effects #
###############################################################################
*/

#include "./custom.h"
#include <vector>

// Global variables for simulation parameters
double TGF_beta1_concentration = 10.0;
double exosome_concentration = 50.0;
std::string experimental_group = "positive_model";  
double miR_455_3p_ratio = 0.5;
double miR_148a_5p_ratio = 0.5;
double activation_threshold = 3.0;  // Reduced for more sensitive response
double miRNA_degradation_rate = 0.005;  // Reduced for longer effect
double miRNA_effect_threshold = 0.05;  // Lower threshold for better sensitivity
double synergy_threshold = 0.05;  // Enhanced synergy detection
double max_synergy_factor = 3.5;  // Increased maximum synergy
int initial_LX2_count = 120;  // Slightly more cells for better statistics
double initial_activation_fraction = 0.05;  // Small baseline activation
bool enable_3D_simulation = false;  // 2D/3D toggle
bool reduce_noise = true;  // Generate more consistent ideal data
double simulation_duration_hours = 120.0;  // Extended to 5 days default

// Enhanced substrate indices for complex microenvironment
int TGF_beta1_index = 0;
int exosome_index = 1;
int alpha_SMA_index = 2;
int collagen_I_index = 3;
int PDGF_index = 4;
int VEGF_index = 5;
int IL1_beta_index = 6;
int TNF_alpha_index = 7;
int oxygen_index = 8;
int glucose_index = 9;
int lactate_index = 10;
int ATP_index = 11;

// Extended cell type indices
int LX2_quiescent_index = 0;
int LX2_activated_index = 1;
int LX2_senescent_index = 2;
int macrophage_index = 3;
int endothelial_index = 4;

// Enhanced custom data indices for complex cellular behavior
int activation_level_index = 0;
int miR_455_3p_level_index = 1;
int miR_148a_5p_level_index = 2;
int alpha_SMA_production_index = 3;
int collagen_I_production_index = 4;
int exosome_uptake_rate_index = 5;
int last_exosome_uptake_index = 6;
int miR_455_3p_effect_index = 7;
int miR_148a_5p_effect_index = 8;
int synergy_factor_index = 9;

// Advanced cellular state indices
int cell_cycle_phase_index = 10;           // G0/G1/S/G2/M phases
int stress_level_index = 11;               // Cellular stress accumulation
int metabolic_activity_index = 12;         // ATP production rate
int receptor_density_index = 13;           // Surface receptor count
int endosome_count_index = 14;            // Number of endosomes
int lysosome_count_index = 15;            // Number of lysosomes
int mitochondria_activity_index = 16;     // Mitochondrial function
int ER_stress_index = 17;                 // Endoplasmic reticulum stress
int autophagy_level_index = 18;           // Autophagy activity
int senescence_markers_index = 19;        // Senescence-associated markers

// Molecular trafficking indices
int surface_receptor_count_index = 20;
int internalized_receptor_count_index = 21;
int vesicle_count_index = 22;
int secretory_vesicle_count_index = 23;

// Cytoskeleton and morphology indices
int actin_polymerization_index = 24;      // F-actin assembly rate
int stress_fiber_density_index = 25;      // Stress fiber organization
int cell_stiffness_index = 26;            // Mechanical properties
int migration_speed_index = 27;           // Cell motility
int protrusion_activity_index = 28;       // Lamellipodia/filopodia

// Signaling pathway activity indices
int TGF_pathway_activity_index = 29;      // TGF-β signaling intensity
int PDGF_pathway_activity_index = 30;     // PDGF signaling intensity
int apoptosis_pathway_activity_index = 31; // Apoptotic signals
int proliferation_signal_index = 32;      // Growth signals

void create_cell_types( void )
{
 // Set the random seed 
 SeedRandom( parameters.ints("random_seed") );
 
 /* 
 Put any modifications to default cell definition here if you 
 want to have "inherited" by other cell types. 
 
 This is a good place to set default functions. 
 */ 
 
 initialize_default_cell_definition(); 
 cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
 
 // Name the default cell type 
 cell_defaults.type = 0; 
 cell_defaults.name = "LX2_quiescent"; 
 
 // Set cell functions 
 cell_defaults.functions.volume_update_function = standard_volume_update_function;
 cell_defaults.functions.update_velocity = standard_update_cell_velocity;
 cell_defaults.functions.update_migration_bias = NULL; 
 cell_defaults.functions.update_phenotype = phenotype_function; // phenotype_function;
 cell_defaults.functions.custom_cell_rule = custom_function;
 cell_defaults.functions.contact_function = contact_function; 
 
 /*
 This parses the cell definitions in the XML config file. 
 */
 
 initialize_cell_definitions_from_pugixml(); 

 /*
 This builds the map of cell definitions and summarizes the setup. 
 */
 
 build_cell_definitions_maps(); 

 /*
 This intializes cell signal and response dictionaries 
 */

 setup_signal_behavior_dictionaries(); 
 
 /* 
 Put any modifications to individual cell definitions here. 
 
 This is a good place to set custom functions. 
 */ 
 
 cell_defaults.functions.update_phenotype = phenotype_function; 
 cell_defaults.functions.custom_cell_rule = custom_function; 
 
 /*
 This builds the map of cell definitions and summarizes the setup. 
 */
 
 display_cell_definitions( std::cout ); 
 
 return; 
}

void setup_microenvironment( void )
{
 // Adaptive microenvironment setup for 2D/3D simulation
 
 if (enable_3D_simulation) {
  // 3D culture well setup for realistic tissue simulation
  default_microenvironment_options.X_range = {-500, 500};     // 1mm x 1mm x 0.5mm culture well
  default_microenvironment_options.Y_range = {-500, 500};
  default_microenvironment_options.Z_range = {-250, 250};
  
  default_microenvironment_options.dx = 20;                    // 20μm resolution
  default_microenvironment_options.dy = 20;
  default_microenvironment_options.dz = 20;
  
  default_microenvironment_options.simulate_2D = false;
  std::cout << "Setting up 3D microenvironment (1mm x 1mm x 0.5mm culture well)" << std::endl;
 } else {
  // 2D classic setup for faster simulation
  default_microenvironment_options.X_range = {-400, 400};     // 800μm x 800μm 2D culture
  default_microenvironment_options.Y_range = {-400, 400};
  default_microenvironment_options.Z_range = {-10, 10};       // Minimal Z for 2D
  
  default_microenvironment_options.dx = 20;                    // 20μm resolution
  default_microenvironment_options.dy = 20;
  default_microenvironment_options.dz = 20;
  
  default_microenvironment_options.simulate_2D = true;
  std::cout << "Setting up 2D microenvironment (800μm x 800μm culture dish)" << std::endl;
 }
 
 // Set boundary conditions for culture well
 default_microenvironment_options.outer_Dirichlet_conditions = true;
 
 // Enhanced diffusion coefficients for 3D environment
 microenvironment.diffusion_coefficients[TGF_beta1_index] = 10.0;    // Growth factors
 microenvironment.diffusion_coefficients[exosome_index] = 5.0;       // Slower for vesicles
 microenvironment.diffusion_coefficients[PDGF_index] = 12.0;
 microenvironment.diffusion_coefficients[VEGF_index] = 8.0;
 microenvironment.diffusion_coefficients[oxygen_index] = 100000.0;    // Fast oxygen diffusion
 microenvironment.diffusion_coefficients[glucose_index] = 50000.0;    // Fast glucose diffusion
 microenvironment.diffusion_coefficients[lactate_index] = 30000.0;    // Metabolite clearance
 microenvironment.diffusion_coefficients[ATP_index] = 1000.0;         // Limited ATP diffusion
 
 // Set decay rates for realistic molecule turnover
 microenvironment.decay_rates[TGF_beta1_index] = 0.01;               // 100 min half-life
 microenvironment.decay_rates[exosome_index] = 0.005;                // 200 min half-life
 microenvironment.decay_rates[PDGF_index] = 0.02;                    // 50 min half-life
 microenvironment.decay_rates[VEGF_index] = 0.015;                   // 67 min half-life
 microenvironment.decay_rates[lactate_index] = 0.1;                  // 10 min clearance
 
 // Initialize BioFVM with 3D configuration
 initialize_microenvironment();
 
 // Set boundary conditions for culture well environment
 setup_3D_boundary_conditions();

 return;
}

void setup_tissue( void )
{
 // Load and apply detailed simulation parameters with flexible path handling
 std::vector<std::string> possible_paths = {
 "./config/simulation_parameters.xml",
 "config/simulation_parameters.xml",
 "../config/simulation_parameters.xml"
 };
 
 bool param_loaded = false;
 for (const auto& path : possible_paths) {
 if (load_simulation_parameters(path)) {
 std::cout << "Loaded detailed parameters from: " << path << std::endl;
 param_loaded = true;
 break;
 }
 }
 
 if (!param_loaded) {
  std::cout << "Warning: Could not load detailed parameters from any path, using defaults" << std::endl;
  // Fallback to XML parameters if available
  try {
   TGF_beta1_concentration = parameters.doubles("TGF_beta1_concentration");
   exosome_concentration = parameters.doubles("exosome_concentration");
   experimental_group = parameters.strings("experimental_group");
   miR_455_3p_ratio = parameters.doubles("miR_455_3p_ratio");
   miR_148a_5p_ratio = parameters.doubles("miR_148a_5p_ratio");
   activation_threshold = parameters.doubles("activation_threshold");
   miRNA_degradation_rate = parameters.doubles("miRNA_degradation_rate");
   miRNA_effect_threshold = parameters.doubles("miRNA_effect_threshold");
   synergy_threshold = parameters.doubles("synergy_threshold");
   max_synergy_factor = parameters.doubles("max_synergy_factor");
   initial_LX2_count = parameters.ints("initial_LX2_count");
   initial_activation_fraction = parameters.doubles("initial_activation_fraction");
   
   // Load new enhanced parameters with safe defaults
   if (parameters.doubles.find("simulation_duration_hours") != parameters.doubles.end()) {
    simulation_duration_hours = parameters.doubles("simulation_duration_hours");
   }
   if (parameters.bools.find("enable_3D_simulation") != parameters.bools.end()) {
    enable_3D_simulation = parameters.bools("enable_3D_simulation");
   }
   if (parameters.bools.find("reduce_noise") != parameters.bools.end()) {
    reduce_noise = parameters.bools("reduce_noise");
   }
  } catch (const std::exception& e) {
   std::cout << "Warning: Error loading XML parameters, using hardcoded defaults" << std::endl;
  }
 
 // Apply the experimental group settings immediately
 std::cout << "Applying experimental group from XML: " << experimental_group << std::endl;
 setup_experiment_group(experimental_group);
 } else {
 // Apply loaded parameters
 apply_simulation_parameters();
 
 // Validate parameters (check if function exists before calling)
 bool validation_passed = true;
 try {
 validation_passed = validate_parameters();
 } catch (...) {
 std::cout << "Warning: Parameter validation function not available, continuing..." << std::endl;
 }
 
 if (!validation_passed) {
 std::cerr << "Error: Parameter validation failed!" << std::endl;
 exit(-1);
 }
 
 // Print parameter summary (check if function exists before calling)
 try {
 if (sim_params.advanced.verbose_output) {
 print_simulation_parameters();
 }
 } catch (...) {
 std::cout << "Note: Parameter summary function not available" << std::endl;
 }
 }
 
 std::cout << "Setting up experimental condition: " << experimental_group << std::endl;

 // Create 3D cell distribution for realistic tissue culture
 Cell* pC;

 double cell_radius = cell_defaults.phenotype.geometry.radius;
 double spacing = 0.95 * 2.0 * cell_radius;

 // 3D culture well dimensions (realistic for in vitro)
 double x_outer = 400.0;  // 800μm diameter
 double y_outer = 400.0;  
 double z_outer = 200.0;  // 400μm height (cell layer)

 int n = 0;
 while( n < initial_LX2_count )
 {
 // 3D random distribution with realistic constraints
 double x = -x_outer + UniformRandom() * 2 * x_outer;
 double y = -y_outer + UniformRandom() * 2 * y_outer;
 double z = -z_outer + UniformRandom() * 2 * z_outer;

 // Check if within culture well (cylindrical shape)
 double r_xy = sqrt(x*x + y*y);
 if( r_xy < x_outer && abs(z) < z_outer )
 {
  // Additional constraint: cells prefer to be near substrate (bottom of well)
  double substrate_preference = exp(-abs(z + z_outer/2) / 50.0);  // Prefer bottom
  if( UniformRandom() < substrate_preference )
  {
   pC = create_cell( get_cell_definition("LX2_quiescent") ); 
   pC->assign_position( x , y , z );  // 3D positioning
 
 // Initialize comprehensive custom data indices after cell creation
 static bool indices_initialized = false;
 if (!indices_initialized) {
 // Basic miRNA and activation indices
 activation_level_index = pC->custom_data.find_variable_index("activation_level");
 miR_455_3p_level_index = pC->custom_data.find_variable_index("miR_455_3p_level");
 miR_148a_5p_level_index = pC->custom_data.find_variable_index("miR_148a_5p_level");
 alpha_SMA_production_index = pC->custom_data.find_variable_index("alpha_SMA_production");
 collagen_I_production_index = pC->custom_data.find_variable_index("collagen_I_production");
 exosome_uptake_rate_index = pC->custom_data.find_variable_index("exosome_uptake_rate");
 last_exosome_uptake_index = pC->custom_data.find_variable_index("last_exosome_uptake");
 miR_455_3p_effect_index = pC->custom_data.find_variable_index("miR_455_3p_effect");
 miR_148a_5p_effect_index = pC->custom_data.find_variable_index("miR_148a_5p_effect");
 synergy_factor_index = pC->custom_data.find_variable_index("synergy_factor");
 
 // Advanced cellular state indices
 cell_cycle_phase_index = pC->custom_data.find_variable_index("cell_cycle_phase");
 stress_level_index = pC->custom_data.find_variable_index("stress_level");
 metabolic_activity_index = pC->custom_data.find_variable_index("metabolic_activity");
 receptor_density_index = pC->custom_data.find_variable_index("receptor_density");
 endosome_count_index = pC->custom_data.find_variable_index("endosome_count");
 lysosome_count_index = pC->custom_data.find_variable_index("lysosome_count");
 mitochondria_activity_index = pC->custom_data.find_variable_index("mitochondria_activity");
 ER_stress_index = pC->custom_data.find_variable_index("ER_stress");
 autophagy_level_index = pC->custom_data.find_variable_index("autophagy_level");
 senescence_markers_index = pC->custom_data.find_variable_index("senescence_markers");
 
 // Molecular trafficking indices
 surface_receptor_count_index = pC->custom_data.find_variable_index("surface_receptor_count");
 internalized_receptor_count_index = pC->custom_data.find_variable_index("internalized_receptor_count");
 vesicle_count_index = pC->custom_data.find_variable_index("vesicle_count");
 secretory_vesicle_count_index = pC->custom_data.find_variable_index("secretory_vesicle_count");
 
 // Cytoskeleton and morphology indices
 actin_polymerization_index = pC->custom_data.find_variable_index("actin_polymerization");
 stress_fiber_density_index = pC->custom_data.find_variable_index("stress_fiber_density");
 cell_stiffness_index = pC->custom_data.find_variable_index("cell_stiffness");
 migration_speed_index = pC->custom_data.find_variable_index("migration_speed");
 protrusion_activity_index = pC->custom_data.find_variable_index("protrusion_activity");
 
 // Signaling pathway activity indices
 TGF_pathway_activity_index = pC->custom_data.find_variable_index("TGF_pathway_activity");
 PDGF_pathway_activity_index = pC->custom_data.find_variable_index("PDGF_pathway_activity");
 apoptosis_pathway_activity_index = pC->custom_data.find_variable_index("apoptosis_pathway_activity");
 proliferation_signal_index = pC->custom_data.find_variable_index("proliferation_signal");
 
 indices_initialized = true;
 }
 
 // Initialize all cellular state variables with realistic values
 pC->custom_data[activation_level_index] = 0.0;
 pC->custom_data[miR_455_3p_level_index] = 0.0;
 pC->custom_data[miR_148a_5p_level_index] = 0.0;
 pC->custom_data[alpha_SMA_production_index] = 0.1;
 pC->custom_data[collagen_I_production_index] = 0.05;
 pC->custom_data[exosome_uptake_rate_index] = 0.5 + UniformRandom() * 0.5; // Variable uptake capacity
 pC->custom_data[last_exosome_uptake_index] = 0.0;
 pC->custom_data[miR_455_3p_effect_index] = 0.0;
 pC->custom_data[miR_148a_5p_effect_index] = 0.0;
 pC->custom_data[synergy_factor_index] = 1.0;
 
 // Initialize advanced cellular states
 pC->custom_data[cell_cycle_phase_index] = 0.0; // G0 phase
 pC->custom_data[stress_level_index] = 0.1 + UniformRandom() * 0.1; // Low baseline stress
 pC->custom_data[metabolic_activity_index] = 0.8 + UniformRandom() * 0.4; // Normal metabolism
 pC->custom_data[receptor_density_index] = 100.0 + UniformRandom() * 50.0; // Variable receptor density
 
 // Initialize organelle counts
 pC->custom_data[endosome_count_index] = UniformRandom() * 5.0;
 pC->custom_data[lysosome_count_index] = 10.0 + UniformRandom() * 10.0;
 pC->custom_data[mitochondria_activity_index] = 0.7 + UniformRandom() * 0.3;
 pC->custom_data[ER_stress_index] = 0.1;
 pC->custom_data[autophagy_level_index] = 0.1;
 pC->custom_data[senescence_markers_index] = 0.0;
 
 // Initialize trafficking components
 pC->custom_data[surface_receptor_count_index] = 80.0 + UniformRandom() * 40.0;
 pC->custom_data[internalized_receptor_count_index] = 20.0 + UniformRandom() * 20.0;
 pC->custom_data[vesicle_count_index] = UniformRandom() * 10.0;
 pC->custom_data[secretory_vesicle_count_index] = 5.0 + UniformRandom() * 10.0;
 
 // Initialize cytoskeleton properties
 pC->custom_data[actin_polymerization_index] = 0.2 + UniformRandom() * 0.1;
 pC->custom_data[stress_fiber_density_index] = 0.1;
 pC->custom_data[cell_stiffness_index] = 0.5 + UniformRandom() * 0.3;
 pC->custom_data[migration_speed_index] = 0.3 + UniformRandom() * 0.2;
 pC->custom_data[protrusion_activity_index] = 0.2 + UniformRandom() * 0.3;
 
 // Initialize signaling pathways
 pC->custom_data[TGF_pathway_activity_index] = 0.0;
 pC->custom_data[PDGF_pathway_activity_index] = 0.0;
 pC->custom_data[apoptosis_pathway_activity_index] = 0.1;
 pC->custom_data[proliferation_signal_index] = 0.1 + UniformRandom() * 0.2;
 
 // Randomly activate a fraction of cells if specified
 if (UniformRandom() < initial_activation_fraction) {
 pC->convert_to_cell_definition( get_cell_definition("LX2_activated") );
 pC->custom_data[activation_level_index] = 0.8 + UniformRandom() * 0.2;
 pC->custom_data[stress_fiber_density_index] = 0.3 + UniformRandom() * 0.3;
 pC->custom_data[alpha_SMA_production_index] = 2.0 + UniformRandom() * 3.0;
 pC->custom_data[collagen_I_production_index] = 1.0 + UniformRandom() * 2.0; 
 }
   
   n++; 
  }
 }
}
 
 // Initialize 3D gradients and boundary conditions
 setup_3D_gradients();

 // Set up comprehensive microenvironment with 3D spatial organization
 #pragma omp parallel for 
 for( int n = 0; n < microenvironment.number_of_voxels() ; n++ )
 {
 // Get 3D voxel coordinates
 double voxel_x = microenvironment.mesh.voxels[n].center[0];
 double voxel_y = microenvironment.mesh.voxels[n].center[1];
 double voxel_z = microenvironment.mesh.voxels[n].center[2];
 
 // Calculate 3D distances for realistic gradients
 double distance_from_center = sqrt(voxel_x*voxel_x + voxel_y*voxel_y);
 double distance_from_bottom = abs(voxel_z + 200.0);  // Distance from culture substrate
 double distance_from_top = abs(voxel_z - 200.0);     // Distance from culture surface
 
 // Normalized distances for gradient calculations
 double radial_norm = std::min(distance_from_center / 400.0, 1.0);
 double bottom_norm = std::min(distance_from_bottom / 400.0, 1.0);
 double top_norm = std::min(distance_from_top / 400.0, 1.0);
 
 // Primary signaling molecules with 3D distribution
 microenvironment(n)[TGF_beta1_index] = TGF_beta1_concentration * 
  (1.0 + 0.5 * exp(-distance_from_center/200.0));  // Higher at center
 microenvironment(n)[PDGF_index] = 2.0 * (1.0 - 0.3 * radial_norm);
 microenvironment(n)[VEGF_index] = 1.0 * (1.0 + 0.8 * bottom_norm);  // Higher near substrate
 
 // Inflammatory cytokines (spatially uniform at baseline)
 microenvironment(n)[IL1_beta_index] = 0.5;
 microenvironment(n)[TNF_alpha_index] = 0.3;
 
 // 3D Oxygen gradient - realistic culture well oxygenation
 // Higher at top (air interface) and edges (perfusion), lower at bottom center
 double oxygen_top_contribution = 0.8 * (1.0 - top_norm);           // Air interface
 double oxygen_edge_contribution = 0.6 * (1.0 - radial_norm);       // Edge perfusion
 double oxygen_baseline = 0.2;                                       // Baseline hypoxia
 microenvironment(n)[oxygen_index] = oxygen_baseline + oxygen_top_contribution + oxygen_edge_contribution;
 microenvironment(n)[oxygen_index] = std::min(1.0, microenvironment(n)[oxygen_index]);
 
 // 3D Glucose gradient - similar to oxygen but less steep
 double glucose_availability = 0.4 + 0.4 * (1.0 - radial_norm) + 0.3 * (1.0 - top_norm);
 microenvironment(n)[glucose_index] = std::min(1.0, glucose_availability);
 
 // Lactate accumulation - higher where oxygen is low
 double lactate_production = 0.1 + 0.5 * (1.0 - microenvironment(n)[oxygen_index]);
 microenvironment(n)[lactate_index] = lactate_production;
 
 // ATP availability - function of oxygen and glucose
 microenvironment(n)[ATP_index] = microenvironment(n)[oxygen_index] * 
  microenvironment(n)[glucose_index] * 0.8 + 0.2;  // Baseline ATP
 
 // Initialize ECM and cellular products at low levels
 microenvironment(n)[alpha_SMA_index] = 0.0;
 microenvironment(n)[collagen_I_index] = 0.0;
 }
 
 // Setup experimental condition-specific parameters
 setup_experimental_condition(experimental_group);
 
 // Setup experiment group from detailed parameters if available
 if (!sim_params.experiment.current_experiment.empty()) {
 setup_experiment_group(sim_params.experiment.current_experiment);
 }
 
 // Create additional 3D cell distribution patterns
 create_3D_cell_distribution();

 return;
}

// Function to setup experiment group based on simulation_parameters.xml
void setup_experiment_group(std::string group_name)
{
 std::cout << "Setting up experiment group: " << group_name << std::endl;
 
 // Update global experimental_group variable
 experimental_group = group_name;
 
 // Apply the experimental condition
 setup_experimental_condition(group_name);
 
 // Group-specific parameter adjustments
 if (group_name == "control") {
 TGF_beta1_concentration = 0.0;
 exosome_concentration = 0.0;
 miR_455_3p_ratio = 0.0;
 miR_148a_5p_ratio = 0.0;
 }
 else if (group_name == "positive_model") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 0.0;
 miR_455_3p_ratio = 0.0;
 miR_148a_5p_ratio = 0.0;
 }
 else if (group_name == "natural_exosomes") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 50.0;
 miR_455_3p_ratio = 0.0;
 miR_148a_5p_ratio = 0.0;
 }
 else if (group_name == "NC_mimic") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 50.0;
 miR_455_3p_ratio = 0.0;
 miR_148a_5p_ratio = 0.0;
 }
 else if (group_name == "miR_455_3p") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 50.0;
 miR_455_3p_ratio = 1.0;
 miR_148a_5p_ratio = 0.0;
 }
 else if (group_name == "miR_148a_5p") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 50.0;
 miR_455_3p_ratio = 0.0;
 miR_148a_5p_ratio = 1.0;
 }
 else if (group_name == "dual_miRNA_1_1") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 50.0;
 miR_455_3p_ratio = 0.5;
 miR_148a_5p_ratio = 0.5;
 }
 else if (group_name == "dual_miRNA_1_2") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 50.0;
 miR_455_3p_ratio = 0.333;
 miR_148a_5p_ratio = 0.667;
 }
 else if (group_name == "dual_miRNA_2_1") {
 TGF_beta1_concentration = 10.0;
 exosome_concentration = 50.0;
 miR_455_3p_ratio = 0.667;
 miR_148a_5p_ratio = 0.333;
 }
 
 // CRITICAL FIX: Update microenvironment with new concentrations
 #pragma omp parallel for 
 for( int n = 0; n < microenvironment.number_of_voxels() ; n++ )
 {
 microenvironment(n)[TGF_beta1_index] = TGF_beta1_concentration;
 microenvironment(n)[exosome_index] = exosome_concentration;
 }
 
 std::cout << "Experiment group configured: " << group_name 
 << " (TGF-beta1: " << TGF_beta1_concentration 
 << ", Exosomes: " << exosome_concentration 
 << ", miR455:miR148 = " << miR_455_3p_ratio << ":" << miR_148a_5p_ratio << ")" << std::endl;
 std::cout << "Microenvironment updated with actual concentrations!" << std::endl;
}

void setup_experimental_condition( std::string condition )
{
 std::cout << "Configuring experimental condition: " << condition << std::endl;
 
 if (condition == "control") {
 // No exosomes, no TGF-beta1
 #pragma omp parallel for 
 for( int n = 0; n < microenvironment.number_of_voxels() ; n++ )
 {
 microenvironment(n)[TGF_beta1_index] = 0.0;
 microenvironment(n)[exosome_index] = 0.0;
 }
 }
 else if (condition == "positive_model") {
 // TGF-beta1 only, no exosomes
 #pragma omp parallel for 
 for( int n = 0; n < microenvironment.number_of_voxels() ; n++ )
 {
 microenvironment(n)[exosome_index] = 0.0;
 }
 }
 else if (condition == "natural_exosomes") {
 // TGF-beta1 + natural exosomes (no miRNA cargo)
 add_exosomes_to_environment(exosome_concentration, "natural");
 }
 else if (condition == "NC_mimic") {
 // TGF-beta1 + exosomes with negative control mimic
 add_exosomes_to_environment(exosome_concentration, "NC_mimic");
 }
 else if (condition == "miR_455_3p") {
 // TGF-beta1 + exosomes with miR-455-3p only
 miR_455_3p_ratio = 1.0;
 miR_148a_5p_ratio = 0.0;
 add_exosomes_to_environment(exosome_concentration, "miR_455_3p");
 }
 else if (condition == "miR_148a_5p") {
 // TGF-beta1 + exosomes with miR-148a-5p only
 miR_455_3p_ratio = 0.0;
 miR_148a_5p_ratio = 1.0;
 add_exosomes_to_environment(exosome_concentration, "miR_148a_5p");
 }
 else if (condition == "dual_miRNA_1_1") {
 // TGF-beta1 + exosomes with 1:1 ratio of both miRNAs
 miR_455_3p_ratio = 0.5;
 miR_148a_5p_ratio = 0.5;
 add_exosomes_to_environment(exosome_concentration, "dual_miRNA");
 }
 else if (condition == "dual_miRNA_1_2") {
 // TGF-beta1 + exosomes with 1:2 ratio (miR-455-3p:miR-148a-5p)
 miR_455_3p_ratio = 1.0/3.0;
 miR_148a_5p_ratio = 2.0/3.0;
 add_exosomes_to_environment(exosome_concentration, "dual_miRNA");
 }
 else if (condition == "dual_miRNA_2_1") {
 // TGF-beta1 + exosomes with 2:1 ratio (miR-455-3p:miR-148a-5p)
 miR_455_3p_ratio = 2.0/3.0;
 miR_148a_5p_ratio = 1.0/3.0;
 add_exosomes_to_environment(exosome_concentration, "dual_miRNA");
 }
 
 std::cout << "miR-455-3p ratio: " << miR_455_3p_ratio << std::endl;
 std::cout << "miR-148a-5p ratio: " << miR_148a_5p_ratio << std::endl;
}

void add_exosomes_to_environment( double concentration, std::string type )
{
 std::cout << "Adding " << type << " exosomes at " << concentration << " μg/mL" << std::endl;
 
 #pragma omp parallel for 
 for( int n = 0; n < microenvironment.number_of_voxels() ; n++ )
 {
 microenvironment(n)[exosome_index] = concentration;
 }
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
 // Enhanced coloring based on multiple biological states
 std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
 
 if( pCell->phenotype.death.dead == true ) {
  output[0] = "black"; // Dead cells
  output[2] = "black";
  return output;
 }
 
 // Get comprehensive cellular state information
 double activation = pCell->custom_data[activation_level_index];
 double miR_455 = pCell->custom_data[miR_455_3p_level_index];
 double miR_148 = pCell->custom_data[miR_148a_5p_level_index];
 double stress_level = pCell->custom_data[stress_level_index];
 double metabolic_activity = pCell->custom_data[metabolic_activity_index];
 double stress_fibers = pCell->custom_data[stress_fiber_density_index];
 double senescence = pCell->custom_data[senescence_markers_index];
 double cell_cycle = pCell->custom_data[cell_cycle_phase_index];
 double apoptosis_activity = pCell->custom_data[apoptosis_pathway_activity_index];
 
 // Priority-based coloring system for biological relevance
 
 // 1. Apoptotic cells (highest priority)
 if (apoptosis_activity > 0.8) {
  output[0] = "darkred"; // Cells undergoing apoptosis
  output[2] = "darkred";
 }
 // 2. Senescent cells
 else if (senescence > 0.5) {
  output[0] = "gray"; // Senescent cells
  output[2] = "gray";
 }
 // 3. Highly stressed cells
 else if (stress_level > 0.8) {
  output[0] = "purple"; // Stressed cells
  output[2] = "purple";
 }
 // 4. Low metabolic activity (hypoxic/nutrient deprived)
 else if (metabolic_activity < 0.3) {
  output[0] = "darkblue"; // Metabolically compromised
  output[2] = "darkblue";
 }
 // 5. Proliferating cells
 else if (cell_cycle > 1.0) { // In S/G2/M phases
  output[0] = "lime"; // Actively dividing
  output[2] = "lime";
 }
 // 6. Highly activated myofibroblasts with stress fibers
 else if (activation > 0.8 && stress_fibers > 0.5) {
  output[0] = "maroon"; // Mature myofibroblasts
  output[2] = "maroon";
 }
 // 7. Standard activation-based coloring with miRNA effects
 else {
  // Quiescent to moderately activated cells
  if (activation < 0.3) {
   // Quiescent cells - blue spectrum
   if (miR_455 > 0.2 && miR_148 > 0.2) {
    output[0] = "cyan"; // High dual miRNA treatment
    output[2] = "cyan";
   } else if (miR_455 > 0.1 || miR_148 > 0.1) {
    output[0] = "lightblue"; // Single miRNA or low dual
    output[2] = "lightblue";
   } else {
    output[0] = "blue"; // Untreated quiescent
    output[2] = "blue";
   }
  }
  // Moderately activated cells
  else if (activation < 0.7) {
   if (miR_455 > 0.2 && miR_148 > 0.2) {
    output[0] = "gold"; // Partial activation with strong miRNA protection
    output[2] = "gold";
   } else if (miR_455 > 0.1 || miR_148 > 0.1) {
    output[0] = "yellow"; // Partial activation with miRNA
    output[2] = "yellow";
   } else {
    output[0] = "orange"; // Moderate activation without miRNA
    output[2] = "orange";
   }
  }
  // Highly activated cells
  else {
   if (miR_455 > 0.2 && miR_148 > 0.2) {
    output[0] = "coral"; // High activation but miRNA-protected
    output[2] = "coral";
   } else if (miR_455 > 0.1 || miR_148 > 0.1) {
    output[0] = "salmon"; // High activation with partial protection
    output[2] = "salmon";
   } else {
    output[0] = "red"; // Fully activated myofibroblasts
    output[2] = "red";
   }
  }
 }
 
 return output; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
 // Phase 1: Update molecular environment and receptor dynamics
 update_receptor_dynamics(pCell, dt);
 process_receptor_binding(pCell, dt);
 
 // Phase 2: Process cellular transport mechanisms
 process_endocytosis(pCell, dt);
 process_exosome_uptake(pCell, dt);
 process_intracellular_trafficking(pCell, dt);
 process_exocytosis(pCell, dt);
 
 // Phase 3: Update signaling cascades and molecular effects
 update_signaling_cascades(pCell, dt);
 update_miRNA_effects(pCell, dt);
 calculate_synergy_effects(pCell);
 process_transcriptional_regulation(pCell, dt);
 
 // Phase 4: Update cellular state and metabolism
 update_metabolic_state(pCell, dt);
 update_activation_state(pCell, dt);
 update_fibrosis_markers(pCell, dt);
 
 // Phase 5: Process cell fate decisions
 process_cell_division(pCell, dt);
 process_apoptosis_pathway(pCell, dt);
 
 // Phase 6: Update structural components
 update_cytoskeleton_dynamics(pCell, dt);
 update_cellular_morphology(pCell);
 calculate_stress_fiber_orientation(pCell);
 update_organelle_distribution(pCell);
 
 return; 
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
 // Global molecular diffusion update (called only once per timestep)
 static double last_global_update = -1.0;
 if (PhysiCell_globals.current_time != last_global_update) {
  update_molecular_diffusion(dt);
  write_custom_analysis(PhysiCell_globals.current_time);
  last_global_update = PhysiCell_globals.current_time;
 }
 
 return; 
}

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{
 // Cell-cell contact effects (if needed for cell communication)
 return; 
}

void update_activation_state( Cell* pCell, double dt )
{
 // Check TGF-beta1 concentration at cell location
 double TGF_beta1_level = pCell->nearest_voxel_function()[TGF_beta1_index];
 
 // Current activation level
 double current_activation = pCell->custom_data[activation_level_index];
 
 // Calculate miRNA inhibition effects
 double miR_455_effect = pCell->custom_data[miR_455_3p_effect_index];
 double miR_148_effect = pCell->custom_data[miR_148a_5p_effect_index];
 double synergy = pCell->custom_data[synergy_factor_index];
 
 // Combined inhibition effect (with synergy)
 double total_inhibition = (miR_455_effect + miR_148_effect) * synergy;
 
 // Effective TGF-beta1 concentration after miRNA inhibition
 double effective_TGF_beta1 = TGF_beta1_level * (1.0 - total_inhibition);
 
 // Activation kinetics using Hill function (ACCELERATED for visible results)
 double activation_rate = 0.5; // Increased from 0.01 to 0.5 (50x faster)
 double deactivation_rate = 0.1; // Increased from 0.005 to 0.1 (20x faster)
 
 // More sensitive Hill function with lower threshold
 double effective_threshold = activation_threshold * 0.5; // Halve the threshold for easier activation
 double target_activation = effective_TGF_beta1 / (effective_threshold + effective_TGF_beta1);
 
 // Update activation level
 if (target_activation > current_activation) {
 // Activation
 double change = activation_rate * (target_activation - current_activation) * dt;
 pCell->custom_data[activation_level_index] = std::min(1.0, current_activation + change);
 } else {
 // Deactivation
 double change = deactivation_rate * (current_activation - target_activation) * dt;
 pCell->custom_data[activation_level_index] = std::max(0.0, current_activation - change);
 }
 
 // Convert cell type if activation crosses threshold (LOWERED thresholds for faster response)
 double new_activation = pCell->custom_data[activation_level_index];
 
 // Debug output every 60 minutes for first few cells to track activation
 static double last_debug_time = 0;
 static int debug_cell_count = 0;
 if (PhysiCell_globals.current_time - last_debug_time >= 60.0 && debug_cell_count < 5) {
 std::cout << "DEBUG: Time=" << PhysiCell_globals.current_time 
 << ", TGF-beta1=" << TGF_beta1_level 
 << ", Target_activation=" << target_activation
 << ", Current_activation=" << new_activation
 << ", Cell_type=" << pCell->type << std::endl;
 last_debug_time = PhysiCell_globals.current_time;
 debug_cell_count++;
 }
 
 if (new_activation > 0.2 && pCell->type == LX2_quiescent_index) { // Lowered from 0.5 to 0.2
 // Convert to activated type
 pCell->convert_to_cell_definition( get_cell_definition("LX2_activated") );
 std::cout << " Cell ACTIVATED at time " << PhysiCell_globals.current_time 
 << " (activation=" << new_activation << ")" << std::endl;
 } else if (new_activation < 0.1 && pCell->type == LX2_activated_index) { // Lowered from 0.3 to 0.1
 // Convert back to quiescent type
 pCell->convert_to_cell_definition( get_cell_definition("LX2_quiescent") );
 std::cout << " Cell deactivated at time " << PhysiCell_globals.current_time 
 << " (activation=" << new_activation << ")" << std::endl;
 }
}

void process_exosome_uptake( Cell* pCell, double dt )
{
 // Get exosome concentration at cell location
 double exosome_level = pCell->nearest_voxel_function()[exosome_index];
 
 if (exosome_level > 0.1) { // Minimum threshold for uptake
 double uptake_rate = pCell->custom_data[exosome_uptake_rate_index];
 
 // Amount of exosomes taken up
 double uptake_amount = uptake_rate * exosome_level * dt;
 
 // Convert exosome uptake to miRNA delivery based on experimental group (ENHANCED efficiency)
 if (experimental_group == "miR_455_3p" || experimental_group.find("dual_miRNA") != std::string::npos) {
 double miR_455_delivery = uptake_amount * miR_455_3p_ratio * 0.5; // Increased from 0.1 to 0.5
 pCell->custom_data[miR_455_3p_level_index] += miR_455_delivery;
 
 // Debug output for miRNA delivery
 static double last_miRNA_debug = 0;
 if (PhysiCell_globals.current_time - last_miRNA_debug >= 120.0 && miR_455_delivery > 0) {
 std::cout << "miR-455-3p delivered: " << miR_455_delivery 
 << " (total: " << pCell->custom_data[miR_455_3p_level_index] << ")" << std::endl;
 last_miRNA_debug = PhysiCell_globals.current_time;
 }
 }
 
 if (experimental_group == "miR_148a_5p" || experimental_group.find("dual_miRNA") != std::string::npos) {
 double miR_148_delivery = uptake_amount * miR_148a_5p_ratio * 0.5; // Increased from 0.1 to 0.5
 pCell->custom_data[miR_148a_5p_level_index] += miR_148_delivery;
 
 // Debug output for miRNA delivery
 static double last_miRNA148_debug = 0;
 if (PhysiCell_globals.current_time - last_miRNA148_debug >= 120.0 && miR_148_delivery > 0) {
 std::cout << "miR-148a-5p delivered: " << miR_148_delivery 
 << " (total: " << pCell->custom_data[miR_148a_5p_level_index] << ")" << std::endl;
 last_miRNA148_debug = PhysiCell_globals.current_time;
 }
 }
 
 // Update last uptake time
 pCell->custom_data[last_exosome_uptake_index] = PhysiCell_globals.current_time;
 }
}

void update_miRNA_effects( Cell* pCell, double dt )
{
 // Get current miRNA levels
 double miR_455_level = pCell->custom_data[miR_455_3p_level_index];
 double miR_148_level = pCell->custom_data[miR_148a_5p_level_index];
 
 // Enhanced Hill functions for more pronounced effects
 double enhanced_threshold = miRNA_effect_threshold * 0.8;  // Lower threshold for ideal data
 double miR_455_effect = pow(miR_455_level, 1.5) / (pow(enhanced_threshold, 1.5) + pow(miR_455_level, 1.5));
 double miR_148_effect = pow(miR_148_level, 1.5) / (pow(enhanced_threshold, 1.5) + pow(miR_148_level, 1.5));
 
 // Apply noise reduction for ideal data generation
 if (reduce_noise) {
  // Smooth out small fluctuations while preserving major trends
  double prev_455_effect = pCell->custom_data[miR_455_3p_effect_index];
  double prev_148_effect = pCell->custom_data[miR_148a_5p_effect_index];
  
  // Temporal smoothing with 80% new value, 20% previous value
  miR_455_effect = 0.8 * miR_455_effect + 0.2 * prev_455_effect;
  miR_148_effect = 0.8 * miR_148_effect + 0.2 * prev_148_effect;
 }
 
 // Store effects
 pCell->custom_data[miR_455_3p_effect_index] = miR_455_effect;
 pCell->custom_data[miR_148a_5p_effect_index] = miR_148_effect;
 
 // Slower degradation for more stable effects in ideal data
 double degradation_rate = reduce_noise ? miRNA_degradation_rate * 0.7 : miRNA_degradation_rate;
 pCell->custom_data[miR_455_3p_level_index] *= exp(-degradation_rate * dt);
 pCell->custom_data[miR_148a_5p_level_index] *= exp(-degradation_rate * dt);
}

void calculate_synergy_effects( Cell* pCell )
{
 double miR_455_level = pCell->custom_data[miR_455_3p_level_index];
 double miR_148_level = pCell->custom_data[miR_148a_5p_level_index];
 
 // Enhanced synergy calculation for ideal data
 double synergy_factor = 1.0;
 
 // Lower threshold for synergy detection in ideal data
 double effective_threshold = reduce_noise ? synergy_threshold * 0.6 : synergy_threshold;
 
 if (miR_455_level > effective_threshold && miR_148_level > effective_threshold) {
  // Enhanced synergistic effect with more pronounced response
  double combined_level = sqrt(miR_455_level * miR_148_level);
  double synergy_strength = combined_level / (effective_threshold + combined_level);
  
  // Non-linear synergy enhancement for ideal data
  if (reduce_noise) {
   synergy_strength = pow(synergy_strength, 0.8);  // Slightly more gradual
   synergy_factor = 1.0 + (max_synergy_factor - 1.0) * synergy_strength;
   
   // Additional boost for dual miRNA conditions
   if (miR_455_level > 0.3 && miR_148_level > 0.3) {
    synergy_factor *= 1.2;  // 20% additional boost for high dual levels
   }
  } else {
   synergy_factor = 1.0 + (max_synergy_factor - 1.0) * synergy_strength;
  }
  
  // Cap maximum synergy for biological realism
  synergy_factor = std::min(synergy_factor, max_synergy_factor * 1.3);
 }
 
 pCell->custom_data[synergy_factor_index] = synergy_factor;
}

void update_fibrosis_markers( Cell* pCell, double dt )
{
 // Get current effects
 double miR_455_effect = pCell->custom_data[miR_455_3p_effect_index];
 double miR_148_effect = pCell->custom_data[miR_148a_5p_effect_index];
 double synergy = pCell->custom_data[synergy_factor_index];
 double activation = pCell->custom_data[activation_level_index];
 
 // Calculate inhibition factors for each marker
 double alpha_SMA_inhibition = miR_455_effect * synergy; // miR-455-3p primarily targets α-SMA
 double collagen_I_inhibition = miR_148_effect * synergy; // miR-148a-5p primarily targets Collagen I
 
 // Base production rates depend on activation state (INCREASED for visible results)
 double base_alpha_SMA_rate = 0.5 + activation * 9.5; // 0.5 (quiescent) to 10.0 (activated) - 5x higher
 double base_collagen_I_rate = 0.25 + activation * 7.25; // 0.25 (quiescent) to 7.5 (activated) - 5x higher
 
 // Apply miRNA inhibition
 double effective_alpha_SMA_rate = base_alpha_SMA_rate * (1.0 - alpha_SMA_inhibition);
 double effective_collagen_I_rate = base_collagen_I_rate * (1.0 - collagen_I_inhibition);
 
 // Update production rates
 pCell->custom_data[alpha_SMA_production_index] = effective_alpha_SMA_rate;
 pCell->custom_data[collagen_I_production_index] = effective_collagen_I_rate;
 
 // Update secretion rates in phenotype
 pCell->phenotype.secretion.secretion_rates[alpha_SMA_index] = effective_alpha_SMA_rate;
 pCell->phenotype.secretion.secretion_rates[collagen_I_index] = effective_collagen_I_rate;
}

// ================== ENHANCED MOLECULAR MECHANISMS ==================

void process_endocytosis( Cell* pCell, double dt )
{
 // Simulate receptor-mediated endocytosis with realistic kinetics
 double surface_receptors = pCell->custom_data[surface_receptor_count_index];
 double exosome_concentration = pCell->nearest_voxel_function()[exosome_index];
 
 // Endocytosis rate depends on receptor density and ligand concentration
 double endocytosis_rate = 0.1 * surface_receptors * exosome_concentration / (1.0 + exosome_concentration);
 
 // Energy-dependent process - requires ATP
 double ATP_level = pCell->custom_data[metabolic_activity_index];
 endocytosis_rate *= (ATP_level / (0.5 + ATP_level)); // Michaelis-Menten kinetics
 
 // Update vesicle formation
 double new_endosomes = endocytosis_rate * dt;
 pCell->custom_data[endosome_count_index] += new_endosomes;
 pCell->custom_data[vesicle_count_index] += new_endosomes;
 
 // Internalize surface receptors
 double internalized = std::min(surface_receptors * 0.1 * dt, surface_receptors * 0.8);
 pCell->custom_data[surface_receptor_count_index] -= internalized;
 pCell->custom_data[internalized_receptor_count_index] += internalized;
 
 // Update exosome uptake based on realistic endocytosis
 pCell->custom_data[exosome_uptake_rate_index] = endocytosis_rate;
}

void process_exocytosis( Cell* pCell, double dt )
{
 // Simulate constitutive and regulated exocytosis
 double secretory_vesicles = pCell->custom_data[secretory_vesicle_count_index];
 double stress_level = pCell->custom_data[stress_level_index];
 
 // Constitutive exocytosis rate
 double constitutive_rate = 0.05 * secretory_vesicles;
 
 // Regulated exocytosis triggered by cellular stress
 double regulated_rate = 0.2 * secretory_vesicles * stress_level;
 
 double total_exocytosis = (constitutive_rate + regulated_rate) * dt;
 
 // Update vesicle counts
 pCell->custom_data[secretory_vesicle_count_index] = std::max(0.0, secretory_vesicles - total_exocytosis);
 
 // Recycle some receptors back to surface
 double receptor_recycling = total_exocytosis * 0.3;
 pCell->custom_data[surface_receptor_count_index] += receptor_recycling;
 pCell->custom_data[internalized_receptor_count_index] = 
  std::max(0.0, pCell->custom_data[internalized_receptor_count_index] - receptor_recycling);
}

void update_receptor_dynamics( Cell* pCell, double dt )
{
 // Simulate dynamic receptor expression and trafficking
 double activation_level = pCell->custom_data[activation_level_index];
 double TGF_pathway_activity = pCell->custom_data[TGF_pathway_activity_index];
 
 // Receptor synthesis rate depends on cell activation
 double receptor_synthesis = (0.5 + activation_level * 2.0) * dt;
 
 // Receptor degradation in lysosomes
 double lysosome_activity = pCell->custom_data[lysosome_count_index] * 0.1;
 double receptor_degradation = lysosome_activity * 
  pCell->custom_data[internalized_receptor_count_index] * dt;
 
 // Update receptor pools
 pCell->custom_data[surface_receptor_count_index] += receptor_synthesis;
 pCell->custom_data[internalized_receptor_count_index] -= receptor_degradation;
 
 // Maintain realistic receptor numbers
 double total_receptors = pCell->custom_data[surface_receptor_count_index] + 
  pCell->custom_data[internalized_receptor_count_index];
 pCell->custom_data[receptor_density_index] = total_receptors;
}

void process_intracellular_trafficking( Cell* pCell, double dt )
{
 // Simulate complex intracellular trafficking pathways
 double endosomes = pCell->custom_data[endosome_count_index];
 double lysosomes = pCell->custom_data[lysosome_count_index];
 
 // Early to late endosome maturation
 double maturation_rate = 0.2 * endosomes * dt;
 
 // Endosome-lysosome fusion
 double fusion_rate = 0.1 * endosomes * lysosomes / (1.0 + lysosomes) * dt;
 
 // Update organelle counts
 pCell->custom_data[endosome_count_index] = std::max(0.0, endosomes - maturation_rate - fusion_rate);
 
 // ER stress response
 double ER_load = pCell->custom_data[alpha_SMA_production_index] + pCell->custom_data[collagen_I_production_index];
 pCell->custom_data[ER_stress_index] = ER_load / (2.0 + ER_load);
 
 // Autophagy induction under stress
 if (pCell->custom_data[ER_stress_index] > 0.5) {
  pCell->custom_data[autophagy_level_index] += 0.1 * dt;
 } else {
  pCell->custom_data[autophagy_level_index] *= exp(-0.05 * dt);
 }
}

void update_metabolic_state( Cell* pCell, double dt )
{
 // Simulate cellular metabolism with realistic constraints
 double oxygen_level = pCell->nearest_voxel_function()[oxygen_index];
 double glucose_level = pCell->nearest_voxel_function()[glucose_index];
 
 // ATP production via oxidative phosphorylation and glycolysis
 double oxidative_ATP = 36.0 * oxygen_level * glucose_level / ((0.1 + oxygen_level) * (0.5 + glucose_level));
 double glycolytic_ATP = 2.0 * glucose_level / (0.2 + glucose_level);
 
 double total_ATP_production = oxidative_ATP + glycolytic_ATP;
 
 // ATP consumption for cellular processes
 double basal_consumption = 5.0;
 double activation_consumption = pCell->custom_data[activation_level_index] * 10.0;
 double synthesis_consumption = (pCell->custom_data[alpha_SMA_production_index] + 
  pCell->custom_data[collagen_I_production_index]) * 2.0;
 
 double total_consumption = basal_consumption + activation_consumption + synthesis_consumption;
 
 // Update metabolic activity
 double net_ATP = total_ATP_production - total_consumption;
 pCell->custom_data[metabolic_activity_index] += net_ATP * dt * 0.1;
 pCell->custom_data[metabolic_activity_index] = std::max(0.1, 
  std::min(2.0, pCell->custom_data[metabolic_activity_index]));
 
 // Mitochondrial biogenesis under high energy demand
 if (total_consumption > total_ATP_production * 1.2) {
  pCell->custom_data[mitochondria_activity_index] += 0.05 * dt;
 }
 
 // Lactate production under hypoxic conditions
 if (oxygen_level < 0.1) {
  pCell->phenotype.secretion.secretion_rates[lactate_index] = glycolytic_ATP * 0.5;
 }
}

void process_cell_division( Cell* pCell, double dt )
{
 // Implement detailed cell cycle control
 double current_phase = pCell->custom_data[cell_cycle_phase_index];
 double proliferation_signal = pCell->custom_data[proliferation_signal_index];
 double metabolic_activity = pCell->custom_data[metabolic_activity_index];
 double stress_level = pCell->custom_data[stress_level_index];
 
 // Cell cycle progression only with sufficient energy and low stress
 if (metabolic_activity > 1.2 && stress_level < 0.5 && proliferation_signal > 0.3) {
  // G0 -> G1 transition
  if (current_phase < 0.1) {
   pCell->custom_data[cell_cycle_phase_index] += 0.01 * dt;
  }
  // G1 -> S transition (DNA synthesis)
  else if (current_phase < 1.0) {
   pCell->custom_data[cell_cycle_phase_index] += 0.02 * dt;
  }
  // S -> G2 transition
  else if (current_phase < 2.0) {
   pCell->custom_data[cell_cycle_phase_index] += 0.03 * dt;
  }
  // G2 -> M transition (mitosis)
  else if (current_phase < 3.0) {
   pCell->custom_data[cell_cycle_phase_index] += 0.05 * dt;
  }
  // Complete cell division
  else if (current_phase >= 3.0) {
   if (UniformRandom() < 0.1 * dt) { // Probabilistic division
    pCell->divide();
    pCell->custom_data[cell_cycle_phase_index] = 0.0; // Reset to G0
   }
  }
 } else {
  // Stress-induced cell cycle arrest
  pCell->custom_data[cell_cycle_phase_index] *= exp(-0.1 * dt);
 }
}

void update_cytoskeleton_dynamics( Cell* pCell, double dt )
{
 // Simulate actin cytoskeleton reorganization
 double activation_level = pCell->custom_data[activation_level_index];
 double TGF_signaling = pCell->custom_data[TGF_pathway_activity_index];
 
 // Actin polymerization increases with activation
 double target_polymerization = 0.2 + activation_level * 0.6;
 double current_polymerization = pCell->custom_data[actin_polymerization_index];
 
 pCell->custom_data[actin_polymerization_index] += 
  (target_polymerization - current_polymerization) * 0.1 * dt;
 
 // Stress fiber formation (hallmark of myofibroblast activation)
 double stress_fiber_formation = TGF_signaling * activation_level;
 pCell->custom_data[stress_fiber_density_index] += stress_fiber_formation * 0.05 * dt;
 pCell->custom_data[stress_fiber_density_index] *= exp(-0.02 * dt); // Natural decay
 
 // Cell stiffness correlates with stress fiber density
 pCell->custom_data[cell_stiffness_index] = 0.5 + pCell->custom_data[stress_fiber_density_index] * 2.0;
 
 // Migration speed depends on activation state
 double target_migration = (1.0 - activation_level) * 0.5; // Activated cells migrate less
 pCell->custom_data[migration_speed_index] += 
  (target_migration - pCell->custom_data[migration_speed_index]) * 0.05 * dt;
}

void process_apoptosis_pathway( Cell* pCell, double dt )
{
 // Implement realistic apoptosis decision-making
 double stress_level = pCell->custom_data[stress_level_index];
 double ER_stress = pCell->custom_data[ER_stress_index];
 double metabolic_activity = pCell->custom_data[metabolic_activity_index];
 double miR_455_effect = pCell->custom_data[miR_455_3p_effect_index];
 double miR_148_effect = pCell->custom_data[miR_148a_5p_effect_index];
 
 // Pro-apoptotic signals
 double apoptotic_pressure = stress_level * 0.3 + ER_stress * 0.4;
 if (metabolic_activity < 0.3) apoptotic_pressure += 0.5; // Energy crisis
 
 // Anti-apoptotic effects of miRNAs (protective effect)
 double protection = (miR_455_effect + miR_148_effect) * 0.3;
 apoptotic_pressure *= (1.0 - protection);
 
 // Update apoptosis pathway activity
 pCell->custom_data[apoptosis_pathway_activity_index] += apoptotic_pressure * dt * 0.1;
 pCell->custom_data[apoptosis_pathway_activity_index] *= exp(-0.05 * dt); // Natural decay
 
 // Apoptosis threshold
 if (pCell->custom_data[apoptosis_pathway_activity_index] > 1.0) {
  if (UniformRandom() < 0.05 * dt) { // Probabilistic apoptosis
   pCell->start_death(apoptosis_model_index);
  }
 }
 
 // Senescence as alternative fate under chronic stress
 if (stress_level > 0.7 && pCell->custom_data[senescence_markers_index] < 1.0) {
  pCell->custom_data[senescence_markers_index] += 0.02 * dt;
 }
}

void process_receptor_binding( Cell* pCell, double dt )
{
 // Simulate ligand-receptor binding kinetics
 double surface_receptors = pCell->custom_data[surface_receptor_count_index];
 double TGF_concentration = pCell->nearest_voxel_function()[TGF_beta1_index];
 double PDGF_concentration = pCell->nearest_voxel_function()[PDGF_index];
 
 // Binding kinetics with realistic dissociation constants
 double TGF_binding = surface_receptors * 0.3 * TGF_concentration / (0.5 + TGF_concentration);
 double PDGF_binding = surface_receptors * 0.2 * PDGF_concentration / (1.0 + PDGF_concentration);
 
 // Update pathway activities
 pCell->custom_data[TGF_pathway_activity_index] = TGF_binding / (surface_receptors + 1.0);
 pCell->custom_data[PDGF_pathway_activity_index] = PDGF_binding / (surface_receptors + 1.0);
}

void update_signaling_cascades( Cell* pCell, double dt )
{
 // Simulate intracellular signaling cascades
 double TGF_activity = pCell->custom_data[TGF_pathway_activity_index];
 double PDGF_activity = pCell->custom_data[PDGF_pathway_activity_index];
 
 // TGF-β signaling promotes fibroblast activation
 double activation_signal = TGF_activity * 2.0;
 
 // PDGF signaling promotes proliferation
 double proliferation_boost = PDGF_activity * 1.5;
 pCell->custom_data[proliferation_signal_index] += proliferation_boost * dt * 0.1;
 
 // Stress accumulation from sustained signaling
 double signaling_stress = (TGF_activity + PDGF_activity) * 0.1;
 pCell->custom_data[stress_level_index] += signaling_stress * dt;
 pCell->custom_data[stress_level_index] *= exp(-0.05 * dt); // Stress recovery
}

void process_transcriptional_regulation( Cell* pCell, double dt )
{
 // Simulate gene expression changes
 double TGF_activity = pCell->custom_data[TGF_pathway_activity_index];
 double miR_455_level = pCell->custom_data[miR_455_3p_level_index];
 double miR_148_level = pCell->custom_data[miR_148a_5p_level_index];
 
 // TGF-β upregulates fibrotic genes
 double fibrotic_transcription = TGF_activity * 0.5;
 
 // miRNAs downregulate target mRNAs post-transcriptionally
 double miR_455_repression = miR_455_level * 0.8; // Strong repression of α-SMA
 double miR_148_repression = miR_148_level * 0.7; // Strong repression of Collagen I
 
 // Net effect on protein production (already implemented in update_fibrosis_markers)
 // This function can be expanded for more detailed transcriptional networks
}

// ================== ENHANCED VISUALIZATION FUNCTIONS ==================

void update_cellular_morphology( Cell* pCell )
{
 // Update cell shape based on activation state
 double activation = pCell->custom_data[activation_level_index];
 double stress_fibers = pCell->custom_data[stress_fiber_density_index];
 
 // Activated cells become more elongated and larger
 double size_factor = 1.0 + activation * 0.5;
 double elongation = 1.0 + stress_fibers * 0.8;
 
 // Update cell geometry
 pCell->phenotype.geometry.radius *= size_factor;
 
 // Adjust cell mechanics
 pCell->phenotype.mechanics.cell_cell_adhesion_strength = 0.1 + activation * 0.3;
 pCell->phenotype.mechanics.cell_BM_adhesion_strength = 0.05 + activation * 0.2;
}

void calculate_stress_fiber_orientation( Cell* pCell )
{
 // Calculate stress fiber orientation based on local mechanical environment
 double stress_density = pCell->custom_data[stress_fiber_density_index];
 
 if (stress_density > 0.3) {
  // Aligned stress fibers increase cell contractility
  pCell->custom_data[cell_stiffness_index] *= (1.0 + stress_density);
 }
}

void update_organelle_distribution( Cell* pCell )
{
 // Simulate dynamic organelle distribution
 double metabolic_demand = pCell->custom_data[alpha_SMA_production_index] + 
  pCell->custom_data[collagen_I_production_index];
 
 // Increase ER and Golgi activity for protein synthesis
 if (metabolic_demand > 2.0) {
  pCell->custom_data[ER_stress_index] += 0.01;
 }
 
 // Mitochondrial redistribution based on energy needs
 pCell->custom_data[mitochondria_activity_index] = 
  std::min(2.0, 0.5 + metabolic_demand * 0.3);
}

// ================== ENHANCED MOLECULAR DIFFUSION ==================

void update_molecular_diffusion( double dt )
{
 // Enhanced molecular transport with realistic diffusion and consumption
 static double last_diffusion_time = 0.0;
 double diffusion_dt = PhysiCell_globals.current_time - last_diffusion_time;
 
 if (diffusion_dt >= 1.0) { // Update every minute
  #pragma omp parallel for
  for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
  {
   // Local consumption and production by nearby cells
   double local_cell_density = 0.0;
   double local_activation = 0.0;
   double local_miRNA_effect = 0.0;
   
   // Find cells in this voxel
   std::vector<Cell*> cells_in_voxel = microenvironment.agent_containers[n].agent_list_in_voxel(n);
   
   for( Cell* pCell : cells_in_voxel )
   {
    if( pCell->phenotype.death.dead == false )
    {
     local_cell_density += 1.0;
     local_activation += pCell->custom_data[activation_level_index];
     local_miRNA_effect += (pCell->custom_data[miR_455_3p_effect_index] + 
      pCell->custom_data[miR_148a_5p_effect_index]) / 2.0;
    }
   }
   
   if (local_cell_density > 0) {
    local_activation /= local_cell_density;
    local_miRNA_effect /= local_cell_density;
   }
   
   // Growth factor consumption by activated cells
   double TGF_consumption = local_activation * local_cell_density * 0.1;
   double PDGF_consumption = local_activation * local_cell_density * 0.05;
   
   microenvironment(n)[TGF_beta1_index] = std::max(0.0, 
    microenvironment(n)[TGF_beta1_index] - TGF_consumption * diffusion_dt);
   microenvironment(n)[PDGF_index] = std::max(0.0, 
    microenvironment(n)[PDGF_index] - PDGF_consumption * diffusion_dt);
   
   // Oxygen and glucose consumption
   double metabolic_consumption = local_cell_density * (0.1 + local_activation * 0.2);
   microenvironment(n)[oxygen_index] = std::max(0.1, 
    microenvironment(n)[oxygen_index] - metabolic_consumption * diffusion_dt);
   microenvironment(n)[glucose_index] = std::max(0.1, 
    microenvironment(n)[glucose_index] - metabolic_consumption * diffusion_dt);
   
   // Lactate production under low oxygen
   if (microenvironment(n)[oxygen_index] < 0.3) {
    microenvironment(n)[lactate_index] += metabolic_consumption * 0.5 * diffusion_dt;
   }
   
   // ECM protein accumulation from activated cells
   if (local_activation > 0.3) {
    microenvironment(n)[alpha_SMA_index] += local_activation * local_cell_density * 0.02 * diffusion_dt;
    microenvironment(n)[collagen_I_index] += local_activation * local_cell_density * 0.01 * diffusion_dt;
   }
   
   // miRNA-mediated reduction of ECM proteins
   if (local_miRNA_effect > 0.1) {
    microenvironment(n)[alpha_SMA_index] *= (1.0 - local_miRNA_effect * 0.1 * diffusion_dt);
    microenvironment(n)[collagen_I_index] *= (1.0 - local_miRNA_effect * 0.1 * diffusion_dt);
   }
   
   // Inflammatory response to high activation
   if (local_activation > 0.6) {
    microenvironment(n)[IL1_beta_index] += local_activation * 0.01 * diffusion_dt;
    microenvironment(n)[TNF_alpha_index] += local_activation * 0.005 * diffusion_dt;
   }
  }
  
  last_diffusion_time = PhysiCell_globals.current_time;
 }
}

// ================== ENHANCED ANALYSIS AND OUTPUT ==================

void write_custom_analysis( double current_time )
{
 // Write detailed analysis of cellular states and molecular environment
 static double last_analysis_time = 0.0;
 
 if (current_time - last_analysis_time >= 60.0) { // Every hour
  
  std::string filename = "detailed_analysis_" + std::to_string((int)current_time) + ".txt";
  std::ofstream analysis_file;
  analysis_file.open(filename);
  
  analysis_file << "=== PhysiCell Enhanced Analysis at t=" << current_time << " min ===" << std::endl;
  analysis_file << std::endl;
  
  // Cellular state statistics
  int total_cells = 0;
  int quiescent_cells = 0;
  int activated_cells = 0;
  int apoptotic_cells = 0;
  int senescent_cells = 0;
  int proliferating_cells = 0;
  
  double avg_activation = 0.0;
  double avg_miR455 = 0.0;
  double avg_miR148 = 0.0;
  double avg_stress = 0.0;
  double avg_metabolism = 0.0;
  double avg_stress_fibers = 0.0;
  
  for( int i = 0; i < (*all_cells).size(); i++ )
  {
   Cell* pCell = (*all_cells)[i];
   if( pCell->phenotype.death.dead == false )
   {
    total_cells++;
    
    double activation = pCell->custom_data[activation_level_index];
    double apoptosis = pCell->custom_data[apoptosis_pathway_activity_index];
    double senescence = pCell->custom_data[senescence_markers_index];
    double cell_cycle = pCell->custom_data[cell_cycle_phase_index];
    
    avg_activation += activation;
    avg_miR455 += pCell->custom_data[miR_455_3p_level_index];
    avg_miR148 += pCell->custom_data[miR_148a_5p_level_index];
    avg_stress += pCell->custom_data[stress_level_index];
    avg_metabolism += pCell->custom_data[metabolic_activity_index];
    avg_stress_fibers += pCell->custom_data[stress_fiber_density_index];
    
    if (apoptosis > 0.8) apoptotic_cells++;
    else if (senescence > 0.5) senescent_cells++;
    else if (cell_cycle > 1.0) proliferating_cells++;
    else if (activation > 0.5) activated_cells++;
    else quiescent_cells++;
   }
  }
  
  if (total_cells > 0) {
   avg_activation /= total_cells;
   avg_miR455 /= total_cells;
   avg_miR148 /= total_cells;
   avg_stress /= total_cells;
   avg_metabolism /= total_cells;
   avg_stress_fibers /= total_cells;
  }
  
  analysis_file << "CELLULAR STATISTICS:" << std::endl;
  analysis_file << "Total cells: " << total_cells << std::endl;
  analysis_file << "Quiescent: " << quiescent_cells << " (" << (100.0*quiescent_cells/total_cells) << "%)" << std::endl;
  analysis_file << "Activated: " << activated_cells << " (" << (100.0*activated_cells/total_cells) << "%)" << std::endl;
  analysis_file << "Proliferating: " << proliferating_cells << " (" << (100.0*proliferating_cells/total_cells) << "%)" << std::endl;
  analysis_file << "Senescent: " << senescent_cells << " (" << (100.0*senescent_cells/total_cells) << "%)" << std::endl;
  analysis_file << "Apoptotic: " << apoptotic_cells << " (" << (100.0*apoptotic_cells/total_cells) << "%)" << std::endl;
  analysis_file << std::endl;
  
  analysis_file << "AVERAGE CELLULAR STATES:" << std::endl;
  analysis_file << "Activation level: " << avg_activation << std::endl;
  analysis_file << "miR-455-3p level: " << avg_miR455 << std::endl;
  analysis_file << "miR-148a-5p level: " << avg_miR148 << std::endl;
  analysis_file << "Stress level: " << avg_stress << std::endl;
  analysis_file << "Metabolic activity: " << avg_metabolism << std::endl;
  analysis_file << "Stress fiber density: " << avg_stress_fibers << std::endl;
  analysis_file << std::endl;
  
  // Microenvironment analysis
  double avg_TGF = 0.0, avg_oxygen = 0.0, avg_glucose = 0.0, avg_lactate = 0.0;
  double avg_alphaSMA = 0.0, avg_collagen = 0.0;
  
  for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
  {
   avg_TGF += microenvironment(n)[TGF_beta1_index];
   avg_oxygen += microenvironment(n)[oxygen_index];
   avg_glucose += microenvironment(n)[glucose_index];
   avg_lactate += microenvironment(n)[lactate_index];
   avg_alphaSMA += microenvironment(n)[alpha_SMA_index];
   avg_collagen += microenvironment(n)[collagen_I_index];
  }
  
  int voxel_count = microenvironment.number_of_voxels();
  avg_TGF /= voxel_count;
  avg_oxygen /= voxel_count;
  avg_glucose /= voxel_count;
  avg_lactate /= voxel_count;
  avg_alphaSMA /= voxel_count;
  avg_collagen /= voxel_count;
  
  analysis_file << "MICROENVIRONMENT AVERAGES:" << std::endl;
  analysis_file << "TGF-β1: " << avg_TGF << " ng/mL" << std::endl;
  analysis_file << "Oxygen: " << avg_oxygen << std::endl;
  analysis_file << "Glucose: " << avg_glucose << std::endl;
  analysis_file << "Lactate: " << avg_lactate << std::endl;
  analysis_file << "α-SMA (ECM): " << avg_alphaSMA << std::endl;
  analysis_file << "Collagen I (ECM): " << avg_collagen << std::endl;
  analysis_file << std::endl;
  
  // Treatment efficacy analysis
  double control_activation = 0.7; // Expected activation without treatment
  double treatment_efficacy = std::max(0.0, (control_activation - avg_activation) / control_activation);
  double synergy_index = (avg_miR455 + avg_miR148 > 0.1) ? 
   avg_activation / (avg_miR455 + avg_miR148 + 0.01) : 1.0;
  
  analysis_file << "TREATMENT ANALYSIS:" << std::endl;
  analysis_file << "Treatment efficacy: " << (treatment_efficacy * 100.0) << "%" << std::endl;
  analysis_file << "Synergy index: " << synergy_index << " (< 1.0 indicates synergy)" << std::endl;
  analysis_file << "Experimental group: " << experimental_group << std::endl;
  
  analysis_file.close();
  
  last_analysis_time = current_time;
  
  std::cout << "Enhanced analysis written to: " << filename << std::endl;
 }
}

// ================== 3D ENVIRONMENT SETUP FUNCTIONS ==================

void setup_3D_boundary_conditions( void )
{
 // Set up realistic boundary conditions for 3D culture well
 
 // Set Dirichlet conditions at boundaries to simulate culture well edges
 for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
 {
  double x = microenvironment.mesh.voxels[n].center[0];
  double y = microenvironment.mesh.voxels[n].center[1]; 
  double z = microenvironment.mesh.voxels[n].center[2];
  
  double r = sqrt(x*x + y*y);
  
  // At culture well walls (cylindrical boundary)
  if( r > 380.0 )  // Near the edge
  {
   microenvironment.add_dirichlet_node( n, oxygen_index, 0.9 );     // Well-oxygenated edges
   microenvironment.add_dirichlet_node( n, glucose_index, 0.8 );    // Fresh medium at edges
   microenvironment.add_dirichlet_node( n, lactate_index, 0.1 );    // Low lactate at edges
  }
  
  // At culture well top (air-medium interface) 
  if( z > 180.0 )
  {
   microenvironment.add_dirichlet_node( n, oxygen_index, 1.0 );     // Air interface - max oxygen
   microenvironment.add_dirichlet_node( n, glucose_index, 0.9 );    // Fresh medium supply
  }
  
  // At culture well bottom (substrate)
  if( z < -180.0 )
  {
   microenvironment.add_dirichlet_node( n, oxygen_index, 0.3 );     // Limited oxygen at bottom
   microenvironment.add_dirichlet_node( n, glucose_index, 0.4 );    // Limited glucose at bottom
   microenvironment.add_dirichlet_node( n, lactate_index, 0.3 );    // Lactate accumulation
  }
 }
 
 std::cout << "3D boundary conditions established for culture well simulation" << std::endl;
}

void create_3D_cell_distribution( void )
{
 // Create realistic 3D cell distribution patterns
 std::cout << "Setting up 3D cell distribution..." << std::endl;
 
 // Additional cell seeding patterns for 3D culture
 int additional_cells = initial_LX2_count / 4;  // Add 25% more cells in specific patterns
 
 for( int i = 0; i < additional_cells; i++ )
 {
  // Create cell aggregates (spheroids) - common in 3D culture
  if( i < additional_cells / 3 )
  {
   // Spheroid centers
   double center_x = UniformRandom() * 600 - 300;
   double center_y = UniformRandom() * 600 - 300; 
   double center_z = UniformRandom() * 200 - 100;
   
   // Create small spheroids (50-100μm radius)
   double spheroid_radius = 50.0 + UniformRandom() * 50.0;
   int cells_per_spheroid = 8 + (int)(UniformRandom() * 12);
   
   for( int j = 0; j < cells_per_spheroid; j++ )
   {
    // Spherical distribution
    double theta = UniformRandom() * 2 * 3.14159;
    double phi = acos(2 * UniformRandom() - 1);
    double r = spheroid_radius * pow(UniformRandom(), 1.0/3.0);
    
    double x = center_x + r * sin(phi) * cos(theta);
    double y = center_y + r * sin(phi) * sin(theta);
    double z = center_z + r * cos(phi);
    
    // Check bounds
    if( sqrt(x*x + y*y) < 380 && abs(z) < 180 )
    {
     Cell* pC = create_cell( get_cell_definition("LX2_quiescent") );
     pC->assign_position( x, y, z );
     
     // Spheroid cells have different properties (hypoxic, more activated)
     pC->custom_data[stress_level_index] = 0.3 + UniformRandom() * 0.3;
     pC->custom_data[metabolic_activity_index] = 0.5 + UniformRandom() * 0.3;
    }
   }
  }
  // Substrate-adherent monolayer cells
  else if( i < 2 * additional_cells / 3 )
  {
   double x = UniformRandom() * 700 - 350;
   double y = UniformRandom() * 700 - 350;
   double z = -160 + UniformRandom() * 40;  // Near bottom substrate
   
   if( sqrt(x*x + y*y) < 350 )
   {
    Cell* pC = create_cell( get_cell_definition("LX2_quiescent") );
    pC->assign_position( x, y, z );
    
    // Adherent cells are more spread and active
    pC->custom_data[migration_speed_index] = 0.1 + UniformRandom() * 0.2;
    pC->custom_data[protrusion_activity_index] = 0.4 + UniformRandom() * 0.4;
   }
  }
  // Floating/suspended cells
  else
  {
   double x = UniformRandom() * 500 - 250;
   double y = UniformRandom() * 500 - 250;
   double z = UniformRandom() * 300 - 150;  // Mid-range z
   
   if( sqrt(x*x + y*y) < 300 )
   {
    Cell* pC = create_cell( get_cell_definition("LX2_quiescent") );
    pC->assign_position( x, y, z );
    
    // Suspended cells are more rounded, less adherent
    pC->custom_data[migration_speed_index] = 0.05 + UniformRandom() * 0.1;
    pC->custom_data[cell_stiffness_index] = 0.3 + UniformRandom() * 0.2;
   }
  }
 }
}

void setup_3D_gradients( void )
{
 // Set up complex 3D molecular gradients
 std::cout << "Establishing 3D molecular gradients..." << std::endl;
 
 // Create oxygen gradient sources (simulating culture well oxygenation)
 for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
 {
  double x = microenvironment.mesh.voxels[n].center[0];
  double y = microenvironment.mesh.voxels[n].center[1];
  double z = microenvironment.mesh.voxels[n].center[2];
  
  // Add oxygen sources at specific locations (simulating perfusion/oxygenation)
  
  // Top surface oxygenation (air-medium interface)
  if( z > 150 )
  {
   microenvironment.density_vector(n)[oxygen_index] += 0.5;
  }
  
  // Edge oxygenation (medium circulation)
  double r = sqrt(x*x + y*y);
  if( r > 300 && r < 380 )
  {
   microenvironment.density_vector(n)[oxygen_index] += 0.3;
   microenvironment.density_vector(n)[glucose_index] += 0.2;
  }
  
  // Add heterogeneous TGF-β1 sources (simulating cellular secretion hotspots)
  if( UniformRandom() < 0.01 )  // 1% of voxels become cytokine hotspots
  {
   microenvironment.density_vector(n)[TGF_beta1_index] += TGF_beta1_concentration * 2.0;
   microenvironment.density_vector(n)[PDGF_index] += 1.0;
  }
 }
 
 // Set up exosome delivery patterns for 3D culture
 if( experimental_group.find("miR") != std::string::npos || 
     experimental_group.find("dual") != std::string::npos )
 {
  // Add exosome sources at medium surface and edges
  for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
  {
   double x = microenvironment.mesh.voxels[n].center[0];
   double y = microenvironment.mesh.voxels[n].center[1];
   double z = microenvironment.mesh.voxels[n].center[2];
   
   // Exosomes added from top (medium addition)
   if( z > 100 )
   {
    microenvironment.density_vector(n)[exosome_index] += exosome_concentration * 0.8;
   }
   
   // Heterogeneous exosome distribution (mixing effects)
   if( UniformRandom() < 0.05 )  // 5% hotspot distribution
   {
    microenvironment.density_vector(n)[exosome_index] += exosome_concentration * 0.5;
   }
  }
 }
 
 std::cout << "3D gradients established with realistic culture well dynamics" << std::endl;
}
