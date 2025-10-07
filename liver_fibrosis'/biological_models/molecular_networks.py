"""
Molecular Networks Module

Implements advanced molecular interaction networks for biological system modeling
including signal transduction pathways, gene regulatory networks, and
miRNA regulation systems.

Key Features:
- TGF-beta signaling pathway modeling
- miRNA regulatory network dynamics
- Boolean-continuous hybrid models
- Synergy effect calculations
- Master equation solvers
- Signal transduction cascades
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable, Union
from dataclasses import dataclass
from enum import Enum
import warnings
from scipy.integrate import odeint
from scipy.optimize import fsolve


class SignalState(Enum):
    """Signal activation states."""
    INACTIVE = 0
    ACTIVE = 1
    HYPERACTIVE = 2


@dataclass
class MolecularSpecies:
    """Container for molecular species properties."""
    name: str
    concentration: float = 0.0
    degradation_rate: float = 0.1
    synthesis_rate: float = 0.0
    diffusion_coefficient: float = 1.0
    molecular_weight: float = 50.0  # kDa
    binding_affinity: Dict[str, float] = None
    
    def __post_init__(self):
        if self.binding_affinity is None:
            self.binding_affinity = {}


@dataclass
class ReactionParameters:
    """Parameters for biochemical reactions."""
    forward_rate: float
    reverse_rate: float = 0.0
    michaelis_constant: float = 1.0
    hill_coefficient: float = 1.0
    cooperativity: float = 1.0
    allosteric_factor: float = 1.0


class TGFBetaSignalingNetwork:
    """
    Comprehensive TGF-beta signaling pathway model including
    Smad-dependent and Smad-independent pathways.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize TGF-beta signaling network."""
        default_params = {
            # Receptor dynamics
            'TGF_beta_binding_rate': 1.0,      # per nM per hour
            'receptor_dissociation_rate': 0.1,  # per hour
            'receptor_internalization_rate': 0.5, # per hour
            
            # Smad pathway
            'Smad2_phosphorylation_rate': 2.0,  # per hour
            'Smad3_phosphorylation_rate': 1.8,  # per hour
            'Smad_dephosphorylation_rate': 0.3, # per hour
            'Smad4_binding_rate': 5.0,          # per hour
            'Smad_complex_dissociation_rate': 0.2, # per hour
            'Smad_nuclear_import_rate': 1.0,    # per hour
            'Smad_nuclear_export_rate': 0.5,    # per hour
            
            # Target gene regulation
            'alpha_SMA_transcription_rate': 0.8, # per hour
            'collagen_transcription_rate': 0.6,  # per hour
            'PAI1_transcription_rate': 1.2,      # per hour
            
            # Non-Smad pathways
            'TAK1_activation_rate': 1.5,        # per hour
            'p38_activation_rate': 2.0,         # per hour
            'JNK_activation_rate': 1.8,         # per hour
            
            # Feedback mechanisms
            'Smad7_induction_rate': 0.4,        # per hour
            'Smad7_inhibition_strength': 2.0,
            'miRNA_inhibition_strength': 3.0
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
            
        # Initialize molecular species
        self.species = self._initialize_species()
        
    def _initialize_species(self) -> Dict[str, MolecularSpecies]:
        """Initialize all molecular species in the pathway."""
        species = {}
        
        # Ligands and receptors
        species['TGF_beta'] = MolecularSpecies('TGF_beta', 0.0, 0.1, 0.0, 100.0, 25.0)
        species['TGF_beta_receptor'] = MolecularSpecies('TGF_beta_receptor', 1.0, 0.05, 0.1, 0.0, 65.0)
        species['TGF_beta_complex'] = MolecularSpecies('TGF_beta_complex', 0.0, 0.2, 0.0, 0.0, 90.0)
        
        # Smad proteins
        species['Smad2'] = MolecularSpecies('Smad2', 1.0, 0.1, 0.05, 10.0, 52.0)
        species['Smad3'] = MolecularSpecies('Smad3', 1.0, 0.1, 0.05, 10.0, 48.0)
        species['Smad4'] = MolecularSpecies('Smad4', 0.8, 0.08, 0.04, 8.0, 60.0)
        species['pSmad2'] = MolecularSpecies('pSmad2', 0.0, 0.3, 0.0, 8.0, 52.0)
        species['pSmad3'] = MolecularSpecies('pSmad3', 0.0, 0.3, 0.0, 8.0, 48.0)
        species['Smad_complex'] = MolecularSpecies('Smad_complex', 0.0, 0.2, 0.0, 5.0, 160.0)
        species['nuclear_Smad_complex'] = MolecularSpecies('nuclear_Smad_complex', 0.0, 0.15, 0.0, 2.0, 160.0)
        
        # Target genes and proteins
        species['alpha_SMA_mRNA'] = MolecularSpecies('alpha_SMA_mRNA', 0.1, 0.5, 0.0, 50.0, 0.5)
        species['alpha_SMA_protein'] = MolecularSpecies('alpha_SMA_protein', 0.2, 0.1, 0.0, 1.0, 42.0)
        species['collagen_mRNA'] = MolecularSpecies('collagen_mRNA', 0.1, 0.4, 0.0, 45.0, 0.6)
        species['collagen_protein'] = MolecularSpecies('collagen_protein', 0.1, 0.05, 0.0, 0.5, 300.0)
        
        # Inhibitory factors
        species['Smad7'] = MolecularSpecies('Smad7', 0.2, 0.2, 0.02, 12.0, 46.0)
        species['miR_455_3p'] = MolecularSpecies('miR_455_3p', 0.0, 0.1, 0.0, 80.0, 0.007)
        species['miR_148a_5p'] = MolecularSpecies('miR_148a_5p', 0.0, 0.1, 0.0, 80.0, 0.007)
        
        # Non-Smad pathway components
        species['TAK1'] = MolecularSpecies('TAK1', 1.0, 0.1, 0.05, 5.0, 67.0)
        species['p38'] = MolecularSpecies('p38', 1.0, 0.1, 0.05, 8.0, 38.0)
        species['JNK'] = MolecularSpecies('JNK', 1.0, 0.1, 0.05, 8.0, 46.0)
        species['pTAK1'] = MolecularSpecies('pTAK1', 0.0, 0.3, 0.0, 5.0, 67.0)
        species['pp38'] = MolecularSpecies('pp38', 0.0, 0.3, 0.0, 8.0, 38.0)
        species['pJNK'] = MolecularSpecies('pJNK', 0.0, 0.3, 0.0, 8.0, 46.0)
        
        return species
    
    def update_network(self, dt: float, external_TGF_beta: float = 0.0) -> Dict[str, float]:
        """
        Update the entire TGF-beta signaling network.
        
        Args:
            dt: Time step (hours)
            external_TGF_beta: External TGF-beta concentration
            
        Returns:
            Updated species concentrations
        """
        # Set external TGF-beta
        self.species['TGF_beta'].concentration = external_TGF_beta
        
        # Calculate reaction rates
        rates = self._calculate_reaction_rates()
        
        # Update species concentrations
        for species_name, species in self.species.items():
            if species_name != 'TGF_beta':  # External species
                concentration_change = rates.get(species_name, 0.0) * dt
                species.concentration = max(0.0, species.concentration + concentration_change)
        
        # Return current concentrations
        return {name: species.concentration for name, species in self.species.items()}
    
    def _calculate_reaction_rates(self) -> Dict[str, float]:
        """Calculate reaction rates for all species."""
        rates = {}
        
        # Get current concentrations
        c = {name: species.concentration for name, species in self.species.items()}
        
        # Receptor binding and activation
        receptor_binding_rate = (self.params['TGF_beta_binding_rate'] * 
                               c['TGF_beta'] * c['TGF_beta_receptor'])
        receptor_dissociation_rate = (self.params['receptor_dissociation_rate'] * 
                                    c['TGF_beta_complex'])
        
        rates['TGF_beta_receptor'] = (-receptor_binding_rate + receptor_dissociation_rate + 
                                    self.species['TGF_beta_receptor'].synthesis_rate - 
                                    self.species['TGF_beta_receptor'].degradation_rate * c['TGF_beta_receptor'])
        
        rates['TGF_beta_complex'] = (receptor_binding_rate - receptor_dissociation_rate - 
                                   self.params['receptor_internalization_rate'] * c['TGF_beta_complex'])
        
        # Smad phosphorylation
        smad2_phosphorylation = (self.params['Smad2_phosphorylation_rate'] * 
                               c['TGF_beta_complex'] * c['Smad2'])
        smad3_phosphorylation = (self.params['Smad3_phosphorylation_rate'] * 
                               c['TGF_beta_complex'] * c['Smad3'])
        
        # Smad dephosphorylation
        smad2_dephosphorylation = self.params['Smad_dephosphorylation_rate'] * c['pSmad2']
        smad3_dephosphorylation = self.params['Smad_dephosphorylation_rate'] * c['pSmad3']
        
        rates['Smad2'] = (-smad2_phosphorylation + smad2_dephosphorylation + 
                         self.species['Smad2'].synthesis_rate - 
                         self.species['Smad2'].degradation_rate * c['Smad2'])
        
        rates['Smad3'] = (-smad3_phosphorylation + smad3_dephosphorylation + 
                         self.species['Smad3'].synthesis_rate - 
                         self.species['Smad3'].degradation_rate * c['Smad3'])
        
        rates['pSmad2'] = (smad2_phosphorylation - smad2_dephosphorylation - 
                          self.species['pSmad2'].degradation_rate * c['pSmad2'])
        
        rates['pSmad3'] = (smad3_phosphorylation - smad3_dephosphorylation - 
                          self.species['pSmad3'].degradation_rate * c['pSmad3'])
        
        # Smad complex formation
        complex_formation = (self.params['Smad4_binding_rate'] * 
                           (c['pSmad2'] + c['pSmad3']) * c['Smad4'])
        complex_dissociation = (self.params['Smad_complex_dissociation_rate'] * 
                              c['Smad_complex'])
        
        rates['Smad4'] = (-complex_formation + complex_dissociation + 
                         self.species['Smad4'].synthesis_rate - 
                         self.species['Smad4'].degradation_rate * c['Smad4'])
        
        rates['Smad_complex'] = (complex_formation - complex_dissociation - 
                               self.params['Smad_nuclear_import_rate'] * c['Smad_complex'] + 
                               self.params['Smad_nuclear_export_rate'] * c['nuclear_Smad_complex'])
        
        # Nuclear translocation
        nuclear_import = self.params['Smad_nuclear_import_rate'] * c['Smad_complex']
        nuclear_export = self.params['Smad_nuclear_export_rate'] * c['nuclear_Smad_complex']
        
        rates['nuclear_Smad_complex'] = (nuclear_import - nuclear_export - 
                                       self.species['nuclear_Smad_complex'].degradation_rate * 
                                       c['nuclear_Smad_complex'])
        
        # Target gene transcription with miRNA inhibition
        miRNA_inhibition_factor = 1.0 / (1.0 + self.params['miRNA_inhibition_strength'] * 
                                        (c['miR_455_3p'] + c['miR_148a_5p']))
        
        alpha_SMA_transcription = (self.params['alpha_SMA_transcription_rate'] * 
                                 c['nuclear_Smad_complex'] * miRNA_inhibition_factor)
        collagen_transcription = (self.params['collagen_transcription_rate'] * 
                                c['nuclear_Smad_complex'] * miRNA_inhibition_factor)
        
        rates['alpha_SMA_mRNA'] = (alpha_SMA_transcription - 
                                 self.species['alpha_SMA_mRNA'].degradation_rate * c['alpha_SMA_mRNA'])
        
        rates['collagen_mRNA'] = (collagen_transcription - 
                                self.species['collagen_mRNA'].degradation_rate * c['collagen_mRNA'])
        
        # Protein translation
        alpha_SMA_translation = 0.5 * c['alpha_SMA_mRNA']  # Translation rate
        collagen_translation = 0.3 * c['collagen_mRNA']
        
        rates['alpha_SMA_protein'] = (alpha_SMA_translation - 
                                    self.species['alpha_SMA_protein'].degradation_rate * 
                                    c['alpha_SMA_protein'])
        
        rates['collagen_protein'] = (collagen_translation - 
                                   self.species['collagen_protein'].degradation_rate * 
                                   c['collagen_protein'])
        
        # Smad7 feedback inhibition
        smad7_induction = (self.params['Smad7_induction_rate'] * c['nuclear_Smad_complex'])
        
        rates['Smad7'] = (smad7_induction - 
                         self.species['Smad7'].degradation_rate * c['Smad7'])
        
        # Non-Smad pathways
        tak1_activation = (self.params['TAK1_activation_rate'] * c['TGF_beta_complex'] * 
                          c['TAK1'] / (1.0 + self.params['Smad7_inhibition_strength'] * c['Smad7']))
        
        rates['TAK1'] = (-tak1_activation + 0.3 * c['pTAK1'] + 
                        self.species['TAK1'].synthesis_rate - 
                        self.species['TAK1'].degradation_rate * c['TAK1'])
        
        rates['pTAK1'] = (tak1_activation - 0.3 * c['pTAK1'] - 
                         self.species['pTAK1'].degradation_rate * c['pTAK1'])
        
        # p38 and JNK activation downstream of TAK1
        p38_activation = self.params['p38_activation_rate'] * c['pTAK1'] * c['p38']
        jnk_activation = self.params['JNK_activation_rate'] * c['pTAK1'] * c['JNK']
        
        rates['p38'] = (-p38_activation + 0.4 * c['pp38'] + 
                       self.species['p38'].synthesis_rate - 
                       self.species['p38'].degradation_rate * c['p38'])
        
        rates['pp38'] = (p38_activation - 0.4 * c['pp38'] - 
                        self.species['pp38'].degradation_rate * c['pp38'])
        
        rates['JNK'] = (-jnk_activation + 0.4 * c['pJNK'] + 
                       self.species['JNK'].synthesis_rate - 
                       self.species['JNK'].degradation_rate * c['JNK'])
        
        rates['pJNK'] = (jnk_activation - 0.4 * c['pJNK'] - 
                        self.species['pJNK'].degradation_rate * c['pJNK'])
        
        return rates
    
    def get_activation_level(self) -> float:
        """Calculate overall pathway activation level."""
        smad_activation = (self.species['nuclear_Smad_complex'].concentration / 
                          (0.1 + self.species['nuclear_Smad_complex'].concentration))
        
        nonsmad_activation = ((self.species['pp38'].concentration + 
                             self.species['pJNK'].concentration) / 
                            (0.2 + self.species['pp38'].concentration + 
                             self.species['pJNK'].concentration))
        
        return 0.7 * smad_activation + 0.3 * nonsmad_activation


class MiRNARegulationNetwork:
    """
    miRNA regulatory network model including miRNA processing,
    target binding, and regulatory effects.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize miRNA regulation network."""
        default_params = {
            # miRNA processing
            'miRNA_processing_rate': 0.8,      # per hour
            'miRNA_degradation_rate': 0.1,     # per hour
            'RISC_loading_rate': 2.0,          # per hour
            'RISC_dissociation_rate': 0.05,    # per hour
            
            # Target binding
            'miR_455_alpha_SMA_binding': 5.0,  # binding strength
            'miR_148_collagen_binding': 4.5,   # binding strength
            'target_repression_efficiency': 0.8,
            
            # Synergy parameters
            'synergy_threshold': 0.1,
            'max_synergy_factor': 3.0,
            'cooperativity_factor': 2.0
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
            
        # Target binding affinities
        self.binding_affinities = {
            'miR_455_3p': {
                'alpha_SMA_mRNA': self.params['miR_455_alpha_SMA_binding'],
                'TGF_beta_receptor_mRNA': 2.0,
                'Smad2_mRNA': 1.5
            },
            'miR_148a_5p': {
                'collagen_mRNA': self.params['miR_148_collagen_binding'],
                'alpha_SMA_mRNA': 2.5,
                'TGF_beta_mRNA': 3.0
            }
        }
    
    def calculate_miRNA_effects(self, miR_455_concentration: float,
                              miR_148_concentration: float,
                              target_mRNAs: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate miRNA regulatory effects on target mRNAs.
        
        Args:
            miR_455_concentration: miR-455-3p concentration
            miR_148_concentration: miR-148a-5p concentration  
            target_mRNAs: Dictionary of target mRNA concentrations
            
        Returns:
            Repression factors for each target
        """
        repression_factors = {}
        
        for target, mRNA_level in target_mRNAs.items():
            # Individual miRNA effects
            miR_455_effect = self._calculate_single_miRNA_effect(
                miR_455_concentration, target, 'miR_455_3p', mRNA_level)
            
            miR_148_effect = self._calculate_single_miRNA_effect(
                miR_148_concentration, target, 'miR_148a_5p', mRNA_level)
            
            # Combined effect with synergy
            combined_effect = self._calculate_synergistic_effect(
                miR_455_effect, miR_148_effect, miR_455_concentration, miR_148_concentration)
            
            repression_factors[target] = combined_effect
            
        return repression_factors
    
    def _calculate_single_miRNA_effect(self, miRNA_concentration: float,
                                     target: str, miRNA_type: str,
                                     mRNA_level: float) -> float:
        """Calculate single miRNA repression effect."""
        if target not in self.binding_affinities.get(miRNA_type, {}):
            return 1.0  # No effect
            
        binding_strength = self.binding_affinities[miRNA_type][target]
        
        # Hill function for miRNA-target interaction
        hill_coefficient = 2.0
        K_d = 0.1  # Dissociation constant
        
        occupancy = (miRNA_concentration**hill_coefficient / 
                    (K_d**hill_coefficient + miRNA_concentration**hill_coefficient))
        
        # Repression factor (1 = no repression, 0 = complete repression)
        repression = 1.0 - (self.params['target_repression_efficiency'] * 
                          binding_strength * occupancy / (1.0 + binding_strength * occupancy))
        
        return max(0.1, repression)  # Minimum 10% activity
    
    def _calculate_synergistic_effect(self, effect1: float, effect2: float,
                                    miR_455_conc: float, miR_148_conc: float) -> float:
        """Calculate synergistic miRNA effects."""
        # Individual effects (1 - repression_factor = repression)
        repression1 = 1.0 - effect1
        repression2 = 1.0 - effect2
        
        # Check if both miRNAs are present above threshold
        if (miR_455_conc > self.params['synergy_threshold'] and 
            miR_148_conc > self.params['synergy_threshold']):
            
            # Synergy factor based on concentrations
            synergy_factor = min(self.params['max_synergy_factor'],
                               1.0 + self.params['cooperativity_factor'] * 
                               min(miR_455_conc, miR_148_conc))
            
            # Enhanced combined repression
            combined_repression = 1.0 - (1.0 - repression1) * (1.0 - repression2)
            enhanced_repression = min(0.95, combined_repression * synergy_factor)
            
            return 1.0 - enhanced_repression
        else:
            # Additive effect without synergy
            combined_repression = repression1 + repression2 - repression1 * repression2
            return 1.0 - combined_repression


class BooleanContinuousModel:
    """
    Hybrid Boolean-continuous model for gene regulatory networks.
    """
    
    def __init__(self, network_topology: Dict[str, Dict[str, str]]):
        """
        Initialize Boolean-continuous model.
        
        Args:
            network_topology: Dictionary defining network structure
                             {target: {regulator: 'activation'/'inhibition'}}
        """
        self.topology = network_topology
        self.genes = list(network_topology.keys())
        self.boolean_states = {gene: False for gene in self.genes}
        self.continuous_levels = {gene: 0.0 for gene in self.genes}
        
        # Thresholds for Boolean switching
        self.activation_thresholds = {gene: 0.5 for gene in self.genes}
        self.deactivation_thresholds = {gene: 0.3 for gene in self.genes}
    
    def update_network(self, dt: float, external_signals: Dict[str, float]) -> Dict[str, float]:
        """
        Update Boolean-continuous network state.
        
        Args:
            dt: Time step
            external_signals: External regulatory signals
            
        Returns:
            Updated continuous levels
        """
        new_boolean_states = {}
        new_continuous_levels = {}
        
        for gene in self.genes:
            # Calculate regulatory input
            regulatory_input = self._calculate_regulatory_input(gene, external_signals)
            
            # Update Boolean state based on thresholds
            current_boolean = self.boolean_states[gene]
            current_continuous = self.continuous_levels[gene]
            
            if not current_boolean and regulatory_input > self.activation_thresholds[gene]:
                new_boolean_states[gene] = True
            elif current_boolean and regulatory_input < self.deactivation_thresholds[gene]:
                new_boolean_states[gene] = False
            else:
                new_boolean_states[gene] = current_boolean
            
            # Update continuous level
            if new_boolean_states[gene]:
                target_level = min(2.0, regulatory_input)
                rate = 1.0  # per hour
            else:
                target_level = 0.0
                rate = 0.5  # per hour
            
            level_change = rate * (target_level - current_continuous) * dt
            new_continuous_levels[gene] = max(0.0, current_continuous + level_change)
        
        # Update states
        self.boolean_states = new_boolean_states
        self.continuous_levels = new_continuous_levels
        
        return self.continuous_levels.copy()
    
    def _calculate_regulatory_input(self, target_gene: str, 
                                  external_signals: Dict[str, float]) -> float:
        """Calculate total regulatory input for a gene."""
        if target_gene not in self.topology:
            return external_signals.get(target_gene, 0.0)
        
        total_input = 0.0
        regulators = self.topology[target_gene]
        
        for regulator, regulation_type in regulators.items():
            # Get regulator level (from network or external)
            if regulator in self.continuous_levels:
                regulator_level = self.continuous_levels[regulator]
            else:
                regulator_level = external_signals.get(regulator, 0.0)
            
            # Apply regulation
            if regulation_type == 'activation':
                total_input += regulator_level
            elif regulation_type == 'inhibition':
                total_input -= regulator_level
        
        # Add external signal
        total_input += external_signals.get(target_gene, 0.0)
        
        return max(0.0, total_input)


class SynergyEffectCalculator:
    """
    Calculator for synergistic effects in multi-drug treatments.
    """
    
    def __init__(self, method: str = "bliss_independence"):
        """
        Initialize synergy calculator.
        
        Args:
            method: "bliss_independence", "loewe_additivity", or "highest_single_agent"
        """
        self.method = method
    
    def calculate_synergy(self, drug_A_effect: float, drug_B_effect: float,
                         combination_effect: float) -> Dict[str, float]:
        """
        Calculate synergy metrics.
        
        Args:
            drug_A_effect: Effect of drug A alone (0-1)
            drug_B_effect: Effect of drug B alone (0-1)
            combination_effect: Effect of combination (0-1)
            
        Returns:
            Dictionary with synergy metrics
        """
        if self.method == "bliss_independence":
            expected_effect = drug_A_effect + drug_B_effect - drug_A_effect * drug_B_effect
        elif self.method == "loewe_additivity":
            # Simplified Loewe model assuming equal potencies
            expected_effect = max(drug_A_effect, drug_B_effect)
        elif self.method == "highest_single_agent":
            expected_effect = max(drug_A_effect, drug_B_effect)
        else:
            raise ValueError(f"Unknown synergy method: {self.method}")
        
        # Calculate synergy metrics
        synergy_score = combination_effect - expected_effect
        synergy_ratio = combination_effect / expected_effect if expected_effect > 0 else 1.0
        
        # Classification
        if synergy_ratio > 1.1:
            synergy_class = "synergistic"
        elif synergy_ratio < 0.9:
            synergy_class = "antagonistic"
        else:
            synergy_class = "additive"
        
        return {
            'expected_effect': expected_effect,
            'observed_effect': combination_effect,
            'synergy_score': synergy_score,
            'synergy_ratio': synergy_ratio,
            'synergy_class': synergy_class
        }
    
    def dose_response_synergy(self, dose_A_range: np.ndarray, dose_B_range: np.ndarray,
                            effect_function_A: Callable, effect_function_B: Callable,
                            combination_function: Callable) -> np.ndarray:
        """
        Calculate synergy across dose ranges.
        
        Args:
            dose_A_range: Range of doses for drug A
            dose_B_range: Range of doses for drug B
            effect_function_A: Function mapping dose A to effect
            effect_function_B: Function mapping dose B to effect
            combination_function: Function mapping (dose_A, dose_B) to effect
            
        Returns:
            Synergy matrix
        """
        synergy_matrix = np.zeros((len(dose_A_range), len(dose_B_range)))
        
        for i, dose_A in enumerate(dose_A_range):
            for j, dose_B in enumerate(dose_B_range):
                effect_A = effect_function_A(dose_A)
                effect_B = effect_function_B(dose_B)
                combination_effect = combination_function(dose_A, dose_B)
                
                synergy_metrics = self.calculate_synergy(effect_A, effect_B, combination_effect)
                synergy_matrix[i, j] = synergy_metrics['synergy_score']
        
        return synergy_matrix


class MasterEquationSolver:
    """
    Solver for master equations in stochastic gene expression.
    """
    
    def __init__(self, max_molecules: int = 100):
        """
        Initialize master equation solver.
        
        Args:
            max_molecules: Maximum number of molecules to consider
        """
        self.max_molecules = max_molecules
        self.state_space = np.arange(max_molecules + 1)
        
    def solve_steady_state(self, birth_rate: float, death_rate: float) -> np.ndarray:
        """
        Solve for steady-state distribution of birth-death process.
        
        Args:
            birth_rate: Molecule production rate
            death_rate: Molecule degradation rate
            
        Returns:
            Steady-state probability distribution
        """
        # For birth-death process: P(n) = (birth_rate/death_rate)^n * P(0) / n!
        # This is a Poisson distribution with lambda = birth_rate/death_rate
        
        if death_rate == 0:
            raise ValueError("Death rate cannot be zero")
        
        lambda_param = birth_rate / death_rate
        
        # Calculate Poisson probabilities
        probabilities = np.zeros(self.max_molecules + 1)
        for n in self.state_space:
            probabilities[n] = (lambda_param**n * np.exp(-lambda_param) / 
                              np.math.factorial(n))
        
        # Normalize (in case of truncation)
        probabilities /= np.sum(probabilities)
        
        return probabilities
    
    def solve_time_dependent(self, initial_distribution: np.ndarray,
                           birth_rate: float, death_rate: float,
                           time_points: np.ndarray) -> np.ndarray:
        """
        Solve time-dependent master equation.
        
        Args:
            initial_distribution: Initial probability distribution
            birth_rate: Molecule production rate
            death_rate: Molecule degradation rate  
            time_points: Time points for solution
            
        Returns:
            Time-dependent probability distributions
        """
        # Build transition rate matrix
        Q = self._build_rate_matrix(birth_rate, death_rate)
        
        # Solve matrix exponential
        from scipy.linalg import expm
        
        solutions = np.zeros((len(time_points), len(initial_distribution)))
        
        for i, t in enumerate(time_points):
            # P(t) = exp(Q*t) * P(0)
            transition_matrix = expm(Q * t)
            solutions[i] = transition_matrix @ initial_distribution
        
        return solutions
    
    def _build_rate_matrix(self, birth_rate: float, death_rate: float) -> np.ndarray:
        """Build transition rate matrix for birth-death process."""
        size = self.max_molecules + 1
        Q = np.zeros((size, size))
        
        for n in range(size):
            # Birth transitions (n -> n+1)
            if n < self.max_molecules:
                Q[n, n+1] = birth_rate
                Q[n, n] -= birth_rate
            
            # Death transitions (n -> n-1)
            if n > 0:
                Q[n, n-1] = n * death_rate
                Q[n, n] -= n * death_rate
        
        return Q


class SignalTransductionPathway:
    """
    Generic signal transduction pathway model with cascades and feedback.
    """
    
    def __init__(self, pathway_structure: List[Dict]):
        """
        Initialize signal transduction pathway.
        
        Args:
            pathway_structure: List of pathway components with parameters
        """
        self.components = pathway_structure
        self.component_states = {comp['name']: 0.0 for comp in pathway_structure}
        
    def update_pathway(self, dt: float, input_signal: float) -> Dict[str, float]:
        """
        Update entire pathway state.
        
        Args:
            dt: Time step
            input_signal: Input signal strength
            
        Returns:
            Updated component states
        """
        new_states = self.component_states.copy()
        
        # Process each component in order
        current_signal = input_signal
        
        for component in self.components:
            name = component['name']
            component_type = component['type']
            params = component['parameters']
            
            if component_type == 'receptor':
                new_states[name] = self._update_receptor(
                    self.component_states[name], current_signal, params, dt)
                
            elif component_type == 'kinase':
                upstream_signal = current_signal
                new_states[name] = self._update_kinase(
                    self.component_states[name], upstream_signal, params, dt)
                
            elif component_type == 'transcription_factor':
                upstream_signal = current_signal
                new_states[name] = self._update_transcription_factor(
                    self.component_states[name], upstream_signal, params, dt)
            
            # Update signal for next component
            current_signal = new_states[name]
        
        self.component_states = new_states
        return self.component_states.copy()
    
    def _update_receptor(self, current_state: float, ligand: float,
                        params: Dict, dt: float) -> float:
        """Update receptor activation state."""
        binding_rate = params.get('binding_rate', 1.0)
        dissociation_rate = params.get('dissociation_rate', 0.1)
        
        # Simple binding kinetics
        binding = binding_rate * ligand * (1.0 - current_state)
        dissociation = dissociation_rate * current_state
        
        change = (binding - dissociation) * dt
        return max(0.0, min(1.0, current_state + change))
    
    def _update_kinase(self, current_state: float, upstream_signal: float,
                      params: Dict, dt: float) -> float:
        """Update kinase activation state."""
        activation_rate = params.get('activation_rate', 2.0)
        deactivation_rate = params.get('deactivation_rate', 0.5)
        
        # Michaelis-Menten-like kinetics
        Km = params.get('Km', 0.5)
        activation = activation_rate * upstream_signal / (Km + upstream_signal)
        deactivation = deactivation_rate * current_state
        
        change = (activation - deactivation) * dt
        return max(0.0, min(2.0, current_state + change))
    
    def _update_transcription_factor(self, current_state: float, upstream_signal: float,
                                   params: Dict, dt: float) -> float:
        """Update transcription factor activation state."""
        nuclear_import_rate = params.get('nuclear_import_rate', 1.0)
        nuclear_export_rate = params.get('nuclear_export_rate', 0.2)
        
        # Nuclear translocation
        import_flux = nuclear_import_rate * upstream_signal
        export_flux = nuclear_export_rate * current_state
        
        change = (import_flux - export_flux) * dt
        return max(0.0, min(1.0, current_state + change))
