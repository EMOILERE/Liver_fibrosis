/*
###############################################################################
# Custom module for iGEM Liver Fibrosis Simulation #
# Models LX-2 hepatic stellate cell response to miRNA-loaded exosomes #
# Includes TGF-beta1 activation and miRNA synergistic effects #
###############################################################################
*/

#include "../../../core/PhysiCell.h"
#include "../../../modules/PhysiCell_standard_modules.h" 
#include "parameter_parser.h"

using namespace BioFVM; 
using namespace PhysiCell;

// Setup functions
void create_cell_types( void );
void setup_tissue( void ); 
void setup_microenvironment( void ); 

// Custom cell behavior functions
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );
void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt );

// Enhanced molecular mechanism functions
void update_miRNA_effects( Cell* pCell, double dt );
void process_exosome_uptake( Cell* pCell, double dt );
void calculate_synergy_effects( Cell* pCell );
void update_activation_state( Cell* pCell, double dt );
void update_fibrosis_markers( Cell* pCell, double dt );

// Advanced cellular behavior functions
void process_endocytosis( Cell* pCell, double dt );
void process_exocytosis( Cell* pCell, double dt );
void update_receptor_dynamics( Cell* pCell, double dt );
void process_intracellular_trafficking( Cell* pCell, double dt );
void update_metabolic_state( Cell* pCell, double dt );
void process_cell_division( Cell* pCell, double dt );
void update_cytoskeleton_dynamics( Cell* pCell, double dt );
void process_apoptosis_pathway( Cell* pCell, double dt );

// Enhanced molecular transport functions
void update_molecular_diffusion( double dt );
void process_receptor_binding( Cell* pCell, double dt );
void update_signaling_cascades( Cell* pCell, double dt );
void process_transcriptional_regulation( Cell* pCell, double dt );

// Advanced visualization functions
void update_cellular_morphology( Cell* pCell );
void calculate_stress_fiber_orientation( Cell* pCell );
void update_organelle_distribution( Cell* pCell );

// Experimental condition setup
void setup_experimental_condition( std::string condition );
void setup_experiment_group( std::string group_name );
void add_exosomes_to_environment( double concentration, std::string type );

// 3D environment setup functions
void setup_3D_boundary_conditions( void );
void create_3D_cell_distribution( void );
void setup_3D_gradients( void );

// Analysis and output functions
std::vector<std::string> my_coloring_function( Cell* pCell );
void write_custom_analysis( double current_time );

// Global variables for simulation parameters
extern double TGF_beta1_concentration;
extern double exosome_concentration;
extern std::string experimental_group;
extern double miR_455_3p_ratio;
extern double miR_148a_5p_ratio;
extern double activation_threshold;
extern double miRNA_degradation_rate;
extern double miRNA_effect_threshold;
extern double synergy_threshold;
extern double max_synergy_factor;
extern int initial_LX2_count;
extern double initial_activation_fraction;

// Enhanced simulation control parameters
extern bool enable_3D_simulation;
extern bool reduce_noise;
extern double simulation_duration_hours;

// Substrate and cell type indices
extern int TGF_beta1_index;
extern int exosome_index;
extern int alpha_SMA_index;
extern int collagen_I_index;
extern int PDGF_index;
extern int VEGF_index;
extern int IL1_beta_index;
extern int TNF_alpha_index;
extern int oxygen_index;
extern int glucose_index;
extern int lactate_index;
extern int ATP_index;

extern int LX2_quiescent_index;
extern int LX2_activated_index;
extern int LX2_senescent_index;
extern int macrophage_index;
extern int endothelial_index;

// Enhanced custom data indices
extern int activation_level_index;
extern int miR_455_3p_level_index;
extern int miR_148a_5p_level_index;
extern int alpha_SMA_production_index;
extern int collagen_I_production_index;
extern int exosome_uptake_rate_index;
extern int last_exosome_uptake_index;
extern int miR_455_3p_effect_index;
extern int miR_148a_5p_effect_index;
extern int synergy_factor_index;

// Advanced cellular state indices
extern int cell_cycle_phase_index;
extern int stress_level_index;
extern int metabolic_activity_index;
extern int receptor_density_index;
extern int endosome_count_index;
extern int lysosome_count_index;
extern int mitochondria_activity_index;
extern int ER_stress_index;
extern int autophagy_level_index;
extern int senescence_markers_index;

// Molecular trafficking indices
extern int surface_receptor_count_index;
extern int internalized_receptor_count_index;
extern int vesicle_count_index;
extern int secretory_vesicle_count_index;

// Cytoskeleton and morphology indices
extern int actin_polymerization_index;
extern int stress_fiber_density_index;
extern int cell_stiffness_index;
extern int migration_speed_index;
extern int protrusion_activity_index;

// Signaling pathway activity indices
extern int TGF_pathway_activity_index;
extern int PDGF_pathway_activity_index;
extern int apoptosis_pathway_activity_index;
extern int proliferation_signal_index;
