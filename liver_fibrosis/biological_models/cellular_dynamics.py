"""
Cellular Dynamics Module

Implements comprehensive cellular behavior models including cell cycle,
metabolism, stress response, and other key cellular processes for
hepatic stellate cell simulation.

Key Features:
- Cell cycle progression with checkpoint control
- Metabolic activity modeling with ATP dynamics
- Cellular stress response and adaptation
- Receptor dynamics and trafficking
- Cytoskeleton remodeling and mechanics
- Cell migration and motility
- Apoptosis and senescence pathways
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable
from dataclasses import dataclass
from enum import Enum
import warnings


class CellCyclePhase(Enum):
    """Cell cycle phases enumeration."""
    G0 = 0  # Quiescent
    G1 = 1  # Gap 1
    S = 2   # Synthesis
    G2 = 3  # Gap 2
    M = 4   # Mitosis


class StressType(Enum):
    """Types of cellular stress."""
    OXIDATIVE = "oxidative"
    ER_STRESS = "er_stress"
    MECHANICAL = "mechanical"
    METABOLIC = "metabolic"
    DNA_DAMAGE = "dna_damage"


@dataclass
class CellularState:
    """Container for cellular state variables."""
    activation_level: float = 0.0
    miR_455_3p_level: float = 0.0
    miR_148a_5p_level: float = 0.0
    alpha_SMA_production: float = 0.0
    collagen_I_production: float = 0.0
    cell_cycle_phase: CellCyclePhase = CellCyclePhase.G0
    stress_level: float = 0.0
    metabolic_activity: float = 1.0
    receptor_density: float = 1.0
    endosome_count: float = 0.0
    lysosome_count: float = 10.0
    mitochondria_activity: float = 1.0
    ER_stress: float = 0.0
    autophagy_level: float = 0.1
    senescence_markers: float = 0.0
    surface_receptor_count: float = 100.0
    internalized_receptor_count: float = 0.0
    vesicle_count: float = 0.0
    secretory_vesicle_count: float = 5.0
    actin_polymerization: float = 0.5
    stress_fiber_density: float = 0.3
    cell_stiffness: float = 1.0
    migration_speed: float = 0.0
    protrusion_activity: float = 0.0
    TGF_pathway_activity: float = 0.0
    PDGF_pathway_activity: float = 0.0
    apoptosis_pathway_activity: float = 0.0
    proliferation_signal: float = 0.0


class CellCycleModel:
    """
    Cell cycle progression model with checkpoint control.
    
    Models G0/G1/S/G2/M transitions with growth factor dependence
    and checkpoint mechanisms.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """
        Initialize cell cycle model.
        
        Args:
            parameters: Model parameters dictionary
        """
        default_params = {
            'G1_duration': 8.0,      # hours
            'S_duration': 6.0,       # hours  
            'G2_duration': 4.0,      # hours
            'M_duration': 1.0,       # hours
            'G1_S_checkpoint_threshold': 0.7,
            'G2_M_checkpoint_threshold': 0.8,
            'growth_factor_sensitivity': 1.0,
            'DNA_damage_sensitivity': 2.0,
            'p53_activation_threshold': 0.5
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
            
        self.phase_durations = {
            CellCyclePhase.G1: self.params['G1_duration'],
            CellCyclePhase.S: self.params['S_duration'],
            CellCyclePhase.G2: self.params['G2_duration'],
            CellCyclePhase.M: self.params['M_duration']
        }
        
    def update_cell_cycle(self, state: CellularState, dt: float,
                         growth_factors: float, DNA_damage: float = 0.0) -> CellularState:
        """
        Update cell cycle progression.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            growth_factors: Growth factor concentration
            DNA_damage: DNA damage level
            
        Returns:
            Updated cellular state
        """
        new_state = state
        
        # G0 -> G1 transition (growth factor dependent)
        if state.cell_cycle_phase == CellCyclePhase.G0:
            G0_exit_probability = self._calculate_G0_exit_rate(growth_factors) * dt
            if np.random.rand() < G0_exit_probability:
                new_state.cell_cycle_phase = CellCyclePhase.G1
                
        # G1 -> S checkpoint
        elif state.cell_cycle_phase == CellCyclePhase.G1:
            # Check G1/S checkpoint
            checkpoint_signal = self._evaluate_G1_S_checkpoint(growth_factors, DNA_damage)
            if checkpoint_signal > self.params['G1_S_checkpoint_threshold']:
                new_state.cell_cycle_phase = CellCyclePhase.S
                
        # S -> G2 transition (automatic after duration)
        elif state.cell_cycle_phase == CellCyclePhase.S:
            # Simplified: automatic progression after duration
            new_state.cell_cycle_phase = CellCyclePhase.G2
            
        # G2 -> M checkpoint  
        elif state.cell_cycle_phase == CellCyclePhase.G2:
            checkpoint_signal = self._evaluate_G2_M_checkpoint(DNA_damage)
            if checkpoint_signal > self.params['G2_M_checkpoint_threshold']:
                new_state.cell_cycle_phase = CellCyclePhase.M
                
        # M -> G1 transition (division)
        elif state.cell_cycle_phase == CellCyclePhase.M:
            # Return to G1 after mitosis
            new_state.cell_cycle_phase = CellCyclePhase.G1
            
        return new_state
    
    def _calculate_G0_exit_rate(self, growth_factors: float) -> float:
        """Calculate rate of G0 exit based on growth factors."""
        # Hill function for growth factor response
        K_gf = 0.5  # Half-saturation constant
        n_gf = 2.0  # Hill coefficient
        
        return self.params['growth_factor_sensitivity'] * (
            growth_factors**n_gf / (K_gf**n_gf + growth_factors**n_gf)
        )
    
    def _evaluate_G1_S_checkpoint(self, growth_factors: float, DNA_damage: float) -> float:
        """Evaluate G1/S checkpoint signal."""
        # Growth promoting signal
        growth_signal = growth_factors / (0.5 + growth_factors)
        
        # DNA damage inhibitory signal
        damage_signal = 1.0 / (1.0 + self.params['DNA_damage_sensitivity'] * DNA_damage)
        
        return growth_signal * damage_signal
    
    def _evaluate_G2_M_checkpoint(self, DNA_damage: float) -> float:
        """Evaluate G2/M checkpoint signal."""
        # Primarily DNA damage checkpoint
        return 1.0 / (1.0 + 5.0 * DNA_damage)


class MetabolismModel:
    """
    Cellular metabolism model with ATP dynamics and metabolic switching.
    
    Models glycolysis, oxidative phosphorylation, and metabolic adaptation
    to environmental conditions.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize metabolism model."""
        default_params = {
            'basal_ATP_consumption': 0.1,     # per hour
            'max_glycolysis_rate': 2.0,       # ATP/hour
            'max_oxphos_rate': 10.0,          # ATP/hour
            'glucose_Km': 1.0,                # mM
            'oxygen_Km': 0.1,                 # mM
            'lactate_production_rate': 0.5,   # per glucose
            'mitochondrial_efficiency': 0.8,
            'stress_metabolic_cost': 0.2,     # additional ATP cost under stress
            'activation_metabolic_cost': 0.3   # additional cost when activated
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_metabolism(self, state: CellularState, dt: float,
                         glucose: float, oxygen: float, lactate: float) -> Tuple[CellularState, Dict[str, float]]:
        """
        Update cellular metabolism.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            glucose: Glucose concentration (mM)
            oxygen: Oxygen concentration (mM)
            lactate: Lactate concentration (mM)
            
        Returns:
            Updated state and metabolic fluxes
        """
        new_state = state
        
        # Calculate ATP production rates
        glycolysis_rate = self._calculate_glycolysis_rate(glucose, lactate)
        oxphos_rate = self._calculate_oxphos_rate(oxygen, state.mitochondria_activity)
        
        total_ATP_production = glycolysis_rate + oxphos_rate
        
        # Calculate ATP consumption
        basal_consumption = self.params['basal_ATP_consumption']
        stress_consumption = self.params['stress_metabolic_cost'] * state.stress_level
        activation_consumption = self.params['activation_metabolic_cost'] * state.activation_level
        
        total_ATP_consumption = basal_consumption + stress_consumption + activation_consumption
        
        # Update metabolic activity (ATP balance)
        net_ATP_change = (total_ATP_production - total_ATP_consumption) * dt
        new_state.metabolic_activity = max(0.1, min(2.0, 
            state.metabolic_activity + 0.1 * net_ATP_change))
        
        # Update mitochondrial activity based on oxygen availability
        new_state.mitochondria_activity = self._update_mitochondrial_activity(
            state.mitochondria_activity, oxygen, dt)
        
        # Calculate metabolic fluxes
        fluxes = {
            'glucose_consumption': glycolysis_rate / 2.0,  # 2 ATP per glucose
            'oxygen_consumption': oxphos_rate / 30.0,      # ~30 ATP per O2
            'lactate_production': glycolysis_rate * self.params['lactate_production_rate'],
            'ATP_net': net_ATP_change / dt
        }
        
        return new_state, fluxes
    
    def _calculate_glycolysis_rate(self, glucose: float, lactate: float) -> float:
        """Calculate glycolytic ATP production rate."""
        # Michaelis-Menten kinetics with product inhibition
        glucose_term = glucose / (self.params['glucose_Km'] + glucose)
        lactate_inhibition = 1.0 / (1.0 + lactate / 5.0)  # Product inhibition
        
        return self.params['max_glycolysis_rate'] * glucose_term * lactate_inhibition
    
    def _calculate_oxphos_rate(self, oxygen: float, mitochondrial_activity: float) -> float:
        """Calculate oxidative phosphorylation ATP production rate."""
        oxygen_term = oxygen / (self.params['oxygen_Km'] + oxygen)
        
        return (self.params['max_oxphos_rate'] * oxygen_term * 
                mitochondrial_activity * self.params['mitochondrial_efficiency'])
    
    def _update_mitochondrial_activity(self, current_activity: float, 
                                     oxygen: float, dt: float) -> float:
        """Update mitochondrial activity based on oxygen availability."""
        target_activity = oxygen / (0.2 + oxygen)  # Oxygen-dependent target
        adaptation_rate = 0.1  # per hour
        
        activity_change = adaptation_rate * (target_activity - current_activity) * dt
        return max(0.1, min(2.0, current_activity + activity_change))


class StressResponseModel:
    """
    Cellular stress response model including oxidative stress,
    ER stress, and adaptive mechanisms.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize stress response model."""
        default_params = {
            'stress_accumulation_rate': 0.1,
            'stress_recovery_rate': 0.05,
            'oxidative_stress_threshold': 0.5,
            'ER_stress_threshold': 0.6,
            'stress_adaptation_rate': 0.02,
            'antioxidant_capacity': 1.0,
            'chaperone_capacity': 1.0,
            'autophagy_induction_threshold': 0.7
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_stress_response(self, state: CellularState, dt: float,
                              stress_stimuli: Dict[StressType, float]) -> CellularState:
        """
        Update cellular stress response.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            stress_stimuli: Dictionary of stress stimuli levels
            
        Returns:
            Updated cellular state
        """
        new_state = state
        
        # Calculate total stress accumulation
        total_stress_input = sum(stress_stimuli.values())
        
        # Update overall stress level
        stress_accumulation = self.params['stress_accumulation_rate'] * total_stress_input
        stress_recovery = self.params['stress_recovery_rate'] * state.stress_level
        
        new_state.stress_level = max(0.0, min(2.0, 
            state.stress_level + (stress_accumulation - stress_recovery) * dt))
        
        # Update ER stress specifically
        er_stress_input = stress_stimuli.get(StressType.ER_STRESS, 0.0)
        protein_load = state.alpha_SMA_production + state.collagen_I_production
        
        er_stress_change = (er_stress_input + 0.1 * protein_load - 
                           self.params['chaperone_capacity'] * state.ER_stress) * dt
        new_state.ER_stress = max(0.0, min(2.0, state.ER_stress + er_stress_change))
        
        # Autophagy induction under stress
        if new_state.stress_level > self.params['autophagy_induction_threshold']:
            autophagy_induction = self.params['stress_adaptation_rate'] * new_state.stress_level * dt
            new_state.autophagy_level = min(1.0, state.autophagy_level + autophagy_induction)
        
        # Stress adaptation (hormesis)
        if state.stress_level > 0.3 and state.stress_level < 0.8:
            adaptation_factor = 1.0 + 0.1 * state.stress_level
            new_state.metabolic_activity *= adaptation_factor
        
        return new_state


class ReceptorDynamicsModel:
    """
    Model for receptor expression, trafficking, and signaling dynamics.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize receptor dynamics model."""
        default_params = {
            'receptor_synthesis_rate': 10.0,    # receptors/hour
            'receptor_degradation_rate': 0.1,   # per hour
            'endocytosis_rate': 0.2,           # per hour
            'recycling_rate': 0.15,            # per hour
            'lysosomal_degradation_rate': 0.05, # per hour
            'activation_upregulation_factor': 2.0,
            'ligand_induced_endocytosis_factor': 3.0
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_receptor_dynamics(self, state: CellularState, dt: float,
                                ligand_concentration: float) -> CellularState:
        """
        Update receptor dynamics.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            ligand_concentration: Extracellular ligand concentration
            
        Returns:
            Updated cellular state
        """
        new_state = state
        
        # Receptor synthesis (activation-dependent)
        synthesis_factor = 1.0 + self.params['activation_upregulation_factor'] * state.activation_level
        receptor_synthesis = self.params['receptor_synthesis_rate'] * synthesis_factor * dt
        
        # Receptor endocytosis (ligand-dependent)
        basal_endocytosis = self.params['endocytosis_rate'] * state.surface_receptor_count
        ligand_induced_endocytosis = (self.params['ligand_induced_endocytosis_factor'] * 
                                    ligand_concentration * state.surface_receptor_count)
        total_endocytosis = (basal_endocytosis + ligand_induced_endocytosis) * dt
        
        # Receptor recycling
        receptor_recycling = self.params['recycling_rate'] * state.internalized_receptor_count * dt
        
        # Lysosomal degradation
        receptor_degradation = (self.params['lysosomal_degradation_rate'] * 
                              state.internalized_receptor_count * dt)
        
        # Update receptor pools
        new_state.surface_receptor_count = max(0.0, 
            state.surface_receptor_count + receptor_synthesis + receptor_recycling - total_endocytosis)
        
        new_state.internalized_receptor_count = max(0.0,
            state.internalized_receptor_count + total_endocytosis - receptor_recycling - receptor_degradation)
        
        # Update total receptor density
        new_state.receptor_density = (new_state.surface_receptor_count + 
                                    new_state.internalized_receptor_count) / 100.0
        
        return new_state


class CytoskeletonModel:
    """
    Model for cytoskeletal dynamics including actin polymerization
    and stress fiber formation.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize cytoskeleton model."""
        default_params = {
            'actin_polymerization_rate': 0.5,   # per hour
            'actin_depolymerization_rate': 0.3, # per hour
            'stress_fiber_assembly_rate': 0.2,  # per hour
            'stress_fiber_disassembly_rate': 0.1, # per hour
            'mechanical_stress_threshold': 0.5,
            'TGF_beta_effect_strength': 2.0,
            'stiffness_actin_coupling': 0.8,
            'migration_speed_coupling': 1.5
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_cytoskeleton(self, state: CellularState, dt: float,
                           mechanical_stress: float) -> CellularState:
        """
        Update cytoskeletal dynamics.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            mechanical_stress: External mechanical stress
            
        Returns:
            Updated cellular state
        """
        new_state = state
        
        # Actin polymerization dynamics
        polymerization_rate = self.params['actin_polymerization_rate']
        depolymerization_rate = self.params['actin_depolymerization_rate']
        
        # TGF-beta enhances actin polymerization
        tgf_enhancement = 1.0 + self.params['TGF_beta_effect_strength'] * state.TGF_pathway_activity
        
        actin_change = ((polymerization_rate * tgf_enhancement - 
                        depolymerization_rate * state.actin_polymerization) * dt)
        
        new_state.actin_polymerization = max(0.0, min(2.0, 
            state.actin_polymerization + actin_change))
        
        # Stress fiber dynamics
        assembly_stimulus = (mechanical_stress > self.params['mechanical_stress_threshold'])
        
        if assembly_stimulus or state.activation_level > 0.3:
            stress_fiber_change = self.params['stress_fiber_assembly_rate'] * dt
        else:
            stress_fiber_change = -self.params['stress_fiber_disassembly_rate'] * dt
        
        new_state.stress_fiber_density = max(0.0, min(2.0,
            state.stress_fiber_density + stress_fiber_change))
        
        # Update cell stiffness based on cytoskeleton
        new_state.cell_stiffness = (1.0 + self.params['stiffness_actin_coupling'] * 
                                  (state.actin_polymerization + state.stress_fiber_density))
        
        # Update migration capacity
        migration_factor = (state.actin_polymerization - 0.5 * state.stress_fiber_density)
        new_state.migration_speed = max(0.0, 
            self.params['migration_speed_coupling'] * migration_factor)
        
        return new_state


class MigrationModel:
    """
    Cell migration model including chemotaxis and random motility.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize migration model."""
        default_params = {
            'random_motility_coefficient': 0.1,  # μm²/min
            'chemotaxis_sensitivity': 1.0,
            'gradient_threshold': 0.01,
            'persistence_time': 30.0,            # minutes
            'speed_activation_coupling': 0.5,
            'contact_inhibition_strength': 2.0
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def calculate_migration_velocity(self, state: CellularState,
                                   gradient: np.ndarray,
                                   cell_density: float) -> np.ndarray:
        """
        Calculate cell migration velocity.
        
        Args:
            state: Current cellular state
            gradient: Chemical gradient vector
            cell_density: Local cell density
            
        Returns:
            Migration velocity vector (μm/min)
        """
        # Base migration speed from cytoskeleton
        base_speed = state.migration_speed * (1.0 + self.params['speed_activation_coupling'] * 
                                            state.activation_level)
        
        # Contact inhibition
        contact_inhibition = 1.0 / (1.0 + self.params['contact_inhibition_strength'] * cell_density)
        effective_speed = base_speed * contact_inhibition
        
        # Random motility component
        random_direction = np.random.normal(0, 1, 2)
        random_direction /= np.linalg.norm(random_direction)
        random_velocity = (np.sqrt(2 * self.params['random_motility_coefficient']) * 
                          random_direction)
        
        # Chemotactic component
        gradient_magnitude = np.linalg.norm(gradient)
        if gradient_magnitude > self.params['gradient_threshold']:
            chemotactic_direction = gradient / gradient_magnitude
            chemotactic_velocity = (self.params['chemotaxis_sensitivity'] * 
                                  gradient_magnitude * chemotactic_direction)
        else:
            chemotactic_velocity = np.zeros(2)
        
        # Combine components
        total_velocity = effective_speed * (random_velocity + chemotactic_velocity)
        
        return total_velocity


class ApoptosisModel:
    """
    Apoptosis pathway model with multiple triggers and execution phases.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize apoptosis model."""
        default_params = {
            'stress_apoptosis_threshold': 1.5,
            'DNA_damage_apoptosis_threshold': 0.8,
            'p53_activation_rate': 0.5,        # per hour
            'caspase_activation_rate': 2.0,    # per hour
            'apoptosis_execution_time': 2.0,   # hours
            'survival_signal_strength': 1.0,
            'mitochondrial_threshold': 0.3
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_apoptosis_pathway(self, state: CellularState, dt: float,
                                DNA_damage: float, survival_signals: float) -> CellularState:
        """
        Update apoptosis pathway activity.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            DNA_damage: DNA damage level
            survival_signals: Growth factor/survival signals
            
        Returns:
            Updated cellular state
        """
        new_state = state
        
        # Calculate apoptotic stimuli
        stress_stimulus = max(0.0, state.stress_level - self.params['stress_apoptosis_threshold'])
        damage_stimulus = max(0.0, DNA_damage - self.params['DNA_damage_apoptosis_threshold'])
        mitochondrial_stimulus = max(0.0, self.params['mitochondrial_threshold'] - 
                                   state.mitochondria_activity)
        
        # Calculate survival signals
        survival_factor = survival_signals * self.params['survival_signal_strength']
        
        # Net apoptotic signal
        pro_apoptotic = stress_stimulus + damage_stimulus + mitochondrial_stimulus
        anti_apoptotic = survival_factor
        
        net_apoptotic_signal = max(0.0, pro_apoptotic - anti_apoptotic)
        
        # Update apoptosis pathway activity
        if net_apoptotic_signal > 0.1:
            apoptosis_increase = self.params['caspase_activation_rate'] * net_apoptotic_signal * dt
            new_state.apoptosis_pathway_activity = min(1.0, 
                state.apoptosis_pathway_activity + apoptosis_increase)
        else:
            # Gradual decay if no stimulus
            decay_rate = 0.1  # per hour
            new_state.apoptosis_pathway_activity = max(0.0,
                state.apoptosis_pathway_activity - decay_rate * dt)
        
        return new_state
    
    def is_cell_apoptotic(self, state: CellularState) -> bool:
        """Check if cell should undergo apoptosis."""
        return state.apoptosis_pathway_activity > 0.8


class SenescenceModel:
    """
    Cellular senescence model including replicative and stress-induced senescence.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize senescence model."""
        default_params = {
            'replicative_senescence_rate': 0.001,  # per division
            'stress_senescence_rate': 0.01,        # per hour under high stress
            'senescence_threshold': 0.7,
            'SASP_production_rate': 0.1,           # per hour
            'telomere_shortening_rate': 0.02,      # per division
            'DNA_damage_accumulation_rate': 0.005   # per hour
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_senescence(self, state: CellularState, dt: float,
                         division_occurred: bool) -> CellularState:
        """
        Update senescence markers.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            division_occurred: Whether cell divided this timestep
            
        Returns:
            Updated cellular state
        """
        new_state = state
        
        # Replicative senescence
        if division_occurred:
            replicative_increment = self.params['replicative_senescence_rate']
            new_state.senescence_markers += replicative_increment
        
        # Stress-induced senescence
        if state.stress_level > 1.0:
            stress_increment = (self.params['stress_senescence_rate'] * 
                              (state.stress_level - 1.0) * dt)
            new_state.senescence_markers += stress_increment
        
        # DNA damage accumulation
        damage_increment = self.params['DNA_damage_accumulation_rate'] * state.stress_level * dt
        new_state.senescence_markers += damage_increment
        
        # Cap senescence markers
        new_state.senescence_markers = min(1.0, new_state.senescence_markers)
        
        return new_state
    
    def is_cell_senescent(self, state: CellularState) -> bool:
        """Check if cell is senescent."""
        return state.senescence_markers > self.params['senescence_threshold']


class CellularBehaviorManager:
    """
    Manager class that coordinates all cellular behavior models.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize cellular behavior manager."""
        self.cell_cycle = CellCycleModel(parameters.get('cell_cycle', {}) if parameters else None)
        self.metabolism = MetabolismModel(parameters.get('metabolism', {}) if parameters else None)
        self.stress_response = StressResponseModel(parameters.get('stress_response', {}) if parameters else None)
        self.receptor_dynamics = ReceptorDynamicsModel(parameters.get('receptor_dynamics', {}) if parameters else None)
        self.cytoskeleton = CytoskeletonModel(parameters.get('cytoskeleton', {}) if parameters else None)
        self.migration = MigrationModel(parameters.get('migration', {}) if parameters else None)
        self.apoptosis = ApoptosisModel(parameters.get('apoptosis', {}) if parameters else None)
        self.senescence = SenescenceModel(parameters.get('senescence', {}) if parameters else None)
    
    def update_all_behaviors(self, state: CellularState, dt: float,
                           environment: Dict[str, float]) -> Tuple[CellularState, Dict[str, float]]:
        """
        Update all cellular behaviors in coordinated manner.
        
        Args:
            state: Current cellular state
            dt: Time step (hours)
            environment: Environmental conditions
            
        Returns:
            Updated state and metabolic fluxes
        """
        new_state = state
        
        # Extract environmental variables
        glucose = environment.get('glucose', 5.0)
        oxygen = environment.get('oxygen', 0.2)
        lactate = environment.get('lactate', 1.0)
        growth_factors = environment.get('growth_factors', 0.5)
        TGF_beta = environment.get('TGF_beta', 0.0)
        mechanical_stress = environment.get('mechanical_stress', 0.0)
        
        # Stress stimuli
        stress_stimuli = {
            StressType.OXIDATIVE: environment.get('oxidative_stress', 0.0),
            StressType.ER_STRESS: environment.get('ER_stress', 0.0),
            StressType.MECHANICAL: mechanical_stress,
            StressType.METABOLIC: max(0.0, 1.0 - new_state.metabolic_activity)
        }
        
        # Update metabolism first (affects other processes)
        new_state, metabolic_fluxes = self.metabolism.update_metabolism(
            new_state, dt, glucose, oxygen, lactate)
        
        # Update stress response
        new_state = self.stress_response.update_stress_response(
            new_state, dt, stress_stimuli)
        
        # Update receptor dynamics
        new_state = self.receptor_dynamics.update_receptor_dynamics(
            new_state, dt, TGF_beta)
        
        # Update cytoskeleton
        new_state = self.cytoskeleton.update_cytoskeleton(
            new_state, dt, mechanical_stress)
        
        # Update cell cycle
        DNA_damage = environment.get('DNA_damage', 0.0)
        new_state = self.cell_cycle.update_cell_cycle(
            new_state, dt, growth_factors, DNA_damage)
        
        # Update apoptosis
        survival_signals = growth_factors
        new_state = self.apoptosis.update_apoptosis_pathway(
            new_state, dt, DNA_damage, survival_signals)
        
        # Update senescence
        division_occurred = (state.cell_cycle_phase == CellCyclePhase.M and 
                           new_state.cell_cycle_phase == CellCyclePhase.G1)
        new_state = self.senescence.update_senescence(
            new_state, dt, division_occurred)
        
        return new_state, metabolic_fluxes
    
    def get_cell_fate(self, state: CellularState) -> str:
        """
        Determine cell fate based on current state.
        
        Args:
            state: Current cellular state
            
        Returns:
            Cell fate string
        """
        if self.apoptosis.is_cell_apoptotic(state):
            return "apoptotic"
        elif self.senescence.is_cell_senescent(state):
            return "senescent"
        elif state.cell_cycle_phase == CellCyclePhase.G0:
            return "quiescent"
        elif state.activation_level > 0.5:
            return "activated"
        else:
            return "proliferative"
