"""
Membrane Transport Module

Implements comprehensive membrane transport systems including endocytosis,
exocytosis, intracellular trafficking, and autophagy mechanisms for
cellular uptake and processing of therapeutic agents.

Key Features:
- Receptor-mediated endocytosis modeling
- Constitutive and regulated exocytosis
- Intracellular vesicle trafficking networks
- Membrane receptor dynamics
- Vesicle transport and fusion
- Autophagy pathway modeling
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable
from dataclasses import dataclass
from enum import Enum
import warnings


class VesicleType(Enum):
    """Types of intracellular vesicles."""
    EARLY_ENDOSOME = "early_endosome"
    LATE_ENDOSOME = "late_endosome"
    LYSOSOME = "lysosome"
    SECRETORY_VESICLE = "secretory_vesicle"
    AUTOPHAGOSOME = "autophagosome"
    AUTOLYSOSOME = "autolysosome"


class TransportDirection(Enum):
    """Direction of membrane transport."""
    ENDOCYTIC = "endocytic"
    EXOCYTIC = "exocytic"
    RECYCLING = "recycling"


@dataclass
class VesicleCompartment:
    """Container for vesicle compartment properties."""
    vesicle_type: VesicleType
    count: float = 0.0
    volume: float = 1.0  # femtoliters
    pH: float = 7.0
    cargo_concentration: Dict[str, float] = None
    membrane_proteins: Dict[str, float] = None
    
    def __post_init__(self):
        if self.cargo_concentration is None:
            self.cargo_concentration = {}
        if self.membrane_proteins is None:
            self.membrane_proteins = {}


@dataclass
class MembraneReceptor:
    """Container for membrane receptor properties."""
    name: str
    surface_count: float = 100.0
    internalized_count: float = 0.0
    synthesis_rate: float = 10.0  # per hour
    degradation_rate: float = 0.1  # per hour
    endocytosis_rate: float = 0.2  # per hour
    recycling_rate: float = 0.15  # per hour
    ligand_affinity: float = 1.0  # nM^-1
    clustering_factor: float = 1.0


class EndocytosisModel:
    """
    Comprehensive endocytosis model including clathrin-mediated,
    caveolin-mediated, and macropinocytosis pathways.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize endocytosis model."""
        default_params = {
            # Clathrin-mediated endocytosis
            'clathrin_pit_formation_rate': 0.5,    # per hour
            'clathrin_pit_maturation_time': 0.02,  # hours (1.2 min)
            'clathrin_cargo_capacity': 100.0,      # molecules per vesicle
            
            # Caveolin-mediated endocytosis
            'caveolin_pit_formation_rate': 0.1,    # per hour
            'caveolin_cargo_capacity': 50.0,       # molecules per vesicle
            
            # Macropinocytosis
            'macropinocytosis_rate': 0.05,         # per hour
            'macropinosome_volume': 10.0,          # femtoliters
            
            # Receptor clustering
            'clustering_threshold': 10.0,          # receptors per cluster
            'clustering_enhancement': 3.0,         # fold increase in endocytosis
            
            # Energy dependence
            'ATP_requirement': 0.1,                # ATP per vesicle
            'temperature_factor': 1.0              # Q10 = 2
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_endocytosis(self, receptors: Dict[str, MembraneReceptor],
                          ligand_concentrations: Dict[str, float],
                          ATP_level: float, temperature: float,
                          dt: float) -> Tuple[Dict[str, MembraneReceptor], Dict[str, VesicleCompartment]]:
        """
        Update endocytosis processes.
        
        Args:
            receptors: Dictionary of membrane receptors
            ligand_concentrations: Extracellular ligand concentrations
            ATP_level: Cellular ATP level (normalized)
            temperature: Temperature (Celsius)
            dt: Time step (hours)
            
        Returns:
            Updated receptors and newly formed vesicles
        """
        updated_receptors = {}
        new_vesicles = {}
        
        # Temperature correction (Q10 = 2)
        temp_factor = 2.0**((temperature - 37.0) / 10.0)
        
        # ATP availability factor
        atp_factor = ATP_level / (0.1 + ATP_level)
        
        for receptor_name, receptor in receptors.items():
            updated_receptor = self._update_single_receptor_endocytosis(
                receptor, ligand_concentrations.get(receptor_name, 0.0),
                atp_factor, temp_factor, dt)
            
            updated_receptors[receptor_name] = updated_receptor
            
            # Calculate vesicle formation
            vesicles = self._calculate_vesicle_formation(
                receptor, ligand_concentrations.get(receptor_name, 0.0),
                atp_factor, temp_factor, dt)
            
            if vesicles:
                new_vesicles[receptor_name] = vesicles
        
        return updated_receptors, new_vesicles
    
    def _update_single_receptor_endocytosis(self, receptor: MembraneReceptor,
                                          ligand_concentration: float,
                                          atp_factor: float, temp_factor: float,
                                          dt: float) -> MembraneReceptor:
        """Update endocytosis for a single receptor type."""
        new_receptor = MembraneReceptor(
            name=receptor.name,
            surface_count=receptor.surface_count,
            internalized_count=receptor.internalized_count,
            synthesis_rate=receptor.synthesis_rate,
            degradation_rate=receptor.degradation_rate,
            endocytosis_rate=receptor.endocytosis_rate,
            recycling_rate=receptor.recycling_rate,
            ligand_affinity=receptor.ligand_affinity,
            clustering_factor=receptor.clustering_factor
        )
        
        # Calculate ligand binding
        bound_fraction = self._calculate_ligand_binding(receptor, ligand_concentration)
        
        # Calculate clustering enhancement
        clustering_enhancement = self._calculate_clustering_enhancement(
            receptor, ligand_concentration)
        
        # Endocytosis rate (ligand-dependent)
        base_endocytosis_rate = receptor.endocytosis_rate * temp_factor * atp_factor
        ligand_enhanced_rate = base_endocytosis_rate * (1.0 + 5.0 * bound_fraction)
        cluster_enhanced_rate = ligand_enhanced_rate * clustering_enhancement
        
        # Calculate endocytosis
        endocytosed_receptors = min(receptor.surface_count,
                                  cluster_enhanced_rate * receptor.surface_count * dt)
        
        # Receptor recycling
        recycled_receptors = receptor.recycling_rate * receptor.internalized_count * dt
        
        # Receptor synthesis and degradation
        synthesized_receptors = receptor.synthesis_rate * dt
        degraded_surface = receptor.degradation_rate * receptor.surface_count * dt
        degraded_internal = 0.05 * receptor.internalized_count * dt  # Lysosomal degradation
        
        # Update receptor counts
        new_receptor.surface_count = max(0.0, 
            receptor.surface_count - endocytosed_receptors + recycled_receptors + 
            synthesized_receptors - degraded_surface)
        
        new_receptor.internalized_count = max(0.0,
            receptor.internalized_count + endocytosed_receptors - recycled_receptors - 
            degraded_internal)
        
        return new_receptor
    
    def _calculate_ligand_binding(self, receptor: MembraneReceptor,
                                ligand_concentration: float) -> float:
        """Calculate fraction of receptors bound to ligand."""
        # Simple Langmuir binding
        return (receptor.ligand_affinity * ligand_concentration / 
                (1.0 + receptor.ligand_affinity * ligand_concentration))
    
    def _calculate_clustering_enhancement(self, receptor: MembraneReceptor,
                                        ligand_concentration: float) -> float:
        """Calculate clustering enhancement factor."""
        bound_receptors = (receptor.surface_count * 
                         self._calculate_ligand_binding(receptor, ligand_concentration))
        
        if bound_receptors > self.params['clustering_threshold']:
            cluster_factor = min(self.params['clustering_enhancement'],
                               1.0 + (bound_receptors / self.params['clustering_threshold'] - 1.0))
            return cluster_factor * receptor.clustering_factor
        else:
            return 1.0
    
    def _calculate_vesicle_formation(self, receptor: MembraneReceptor,
                                   ligand_concentration: float,
                                   atp_factor: float, temp_factor: float,
                                   dt: float) -> Optional[VesicleCompartment]:
        """Calculate formation of endocytic vesicles."""
        bound_fraction = self._calculate_ligand_binding(receptor, ligand_concentration)
        
        # Estimate number of vesicles formed
        vesicle_formation_rate = (self.params['clathrin_pit_formation_rate'] * 
                                temp_factor * atp_factor * (1.0 + bound_fraction))
        
        vesicles_formed = vesicle_formation_rate * dt
        
        if vesicles_formed > 0.01:  # Threshold for vesicle formation
            # Create early endosome
            vesicle = VesicleCompartment(
                vesicle_type=VesicleType.EARLY_ENDOSOME,
                count=vesicles_formed,
                volume=0.1,  # femtoliters
                pH=6.5,
                cargo_concentration={receptor.name: bound_fraction * 10.0},
                membrane_proteins={receptor.name: vesicles_formed * 5.0}
            )
            return vesicle
        
        return None


class ExocytosisModel:
    """
    Exocytosis model including constitutive and regulated secretion pathways.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize exocytosis model."""
        default_params = {
            # Constitutive exocytosis
            'constitutive_rate': 0.1,              # per hour
            'constitutive_cargo_release': 0.8,     # fraction released
            
            # Regulated exocytosis
            'calcium_threshold': 0.5,              # μM
            'calcium_sensitivity': 2.0,            # Hill coefficient
            'regulated_rate_max': 2.0,             # per hour
            'regulated_cargo_release': 0.95,       # fraction released
            
            # Vesicle docking and priming
            'docking_rate': 1.0,                   # per hour
            'priming_rate': 0.5,                   # per hour
            'fusion_probability': 0.8,             # per calcium spike
            
            # SNARE complex formation
            'SNARE_assembly_rate': 3.0,            # per hour
            'SNARE_disassembly_rate': 0.2,         # per hour
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_exocytosis(self, vesicles: Dict[str, VesicleCompartment],
                         calcium_concentration: float,
                         stress_signals: float,
                         dt: float) -> Tuple[Dict[str, VesicleCompartment], Dict[str, float]]:
        """
        Update exocytosis processes.
        
        Args:
            vesicles: Dictionary of vesicle compartments
            calcium_concentration: Intracellular calcium level (μM)
            stress_signals: Cellular stress level (0-1)
            dt: Time step (hours)
            
        Returns:
            Updated vesicles and secreted cargo amounts
        """
        updated_vesicles = {}
        secreted_cargo = {}
        
        for vesicle_name, vesicle in vesicles.items():
            if vesicle.vesicle_type == VesicleType.SECRETORY_VESICLE:
                updated_vesicle, cargo = self._process_secretory_vesicle(
                    vesicle, calcium_concentration, stress_signals, dt)
                
                updated_vesicles[vesicle_name] = updated_vesicle
                if cargo:
                    secreted_cargo.update(cargo)
        
        return updated_vesicles, secreted_cargo
    
    def _process_secretory_vesicle(self, vesicle: VesicleCompartment,
                                 calcium: float, stress: float,
                                 dt: float) -> Tuple[VesicleCompartment, Dict[str, float]]:
        """Process exocytosis for secretory vesicles."""
        # Constitutive exocytosis
        constitutive_fusion_rate = self.params['constitutive_rate'] * dt
        
        # Regulated exocytosis (calcium-dependent)
        calcium_factor = (calcium**self.params['calcium_sensitivity'] / 
                         (self.params['calcium_threshold']**self.params['calcium_sensitivity'] + 
                          calcium**self.params['calcium_sensitivity']))
        
        regulated_fusion_rate = (self.params['regulated_rate_max'] * 
                               calcium_factor * (1.0 + stress) * dt)
        
        total_fusion_rate = constitutive_fusion_rate + regulated_fusion_rate
        
        # Calculate vesicles that undergo fusion
        fused_vesicles = min(vesicle.count, total_fusion_rate * vesicle.count)
        
        # Calculate cargo release
        if fused_vesicles > 0:
            if regulated_fusion_rate > constitutive_fusion_rate:
                # Predominantly regulated
                cargo_release_efficiency = self.params['regulated_cargo_release']
            else:
                # Predominantly constitutive
                cargo_release_efficiency = self.params['constitutive_cargo_release']
            
            secreted_cargo = {}
            for cargo_name, concentration in vesicle.cargo_concentration.items():
                amount_secreted = (fused_vesicles * vesicle.volume * 
                                 concentration * cargo_release_efficiency)
                secreted_cargo[cargo_name] = amount_secreted
        else:
            secreted_cargo = {}
        
        # Update vesicle count
        updated_vesicle = VesicleCompartment(
            vesicle_type=vesicle.vesicle_type,
            count=max(0.0, vesicle.count - fused_vesicles),
            volume=vesicle.volume,
            pH=vesicle.pH,
            cargo_concentration=vesicle.cargo_concentration.copy(),
            membrane_proteins=vesicle.membrane_proteins.copy()
        )
        
        return updated_vesicle, secreted_cargo


class IntracellularTraffickingModel:
    """
    Model for intracellular vesicle trafficking including endosome maturation,
    lysosomal fusion, and recycling pathways.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize intracellular trafficking model."""
        default_params = {
            # Endosome maturation
            'early_to_late_endosome_rate': 0.3,    # per hour
            'late_endosome_maturation_time': 0.5,  # hours
            
            # Lysosomal fusion
            'endosome_lysosome_fusion_rate': 0.4,  # per hour
            'lysosomal_degradation_rate': 0.8,     # per hour
            
            # Recycling pathways
            'fast_recycling_rate': 0.6,            # per hour (direct to surface)
            'slow_recycling_rate': 0.2,            # per hour (via recycling endosomes)
            'recycling_endosome_residence_time': 1.0, # hours
            
            # pH changes
            'pH_acidification_rate': 2.0,          # pH units per hour
            'target_early_endosome_pH': 6.5,
            'target_late_endosome_pH': 5.5,
            'target_lysosome_pH': 4.5,
            
            # Cargo sorting
            'cargo_degradation_efficiency': 0.9,
            'receptor_recycling_efficiency': 0.8
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_trafficking(self, vesicles: Dict[str, VesicleCompartment],
                          dt: float) -> Dict[str, VesicleCompartment]:
        """
        Update intracellular trafficking processes.
        
        Args:
            vesicles: Dictionary of vesicle compartments
            dt: Time step (hours)
            
        Returns:
            Updated vesicle compartments
        """
        updated_vesicles = {}
        
        # Process each vesicle type
        for vesicle_name, vesicle in vesicles.items():
            if vesicle.vesicle_type == VesicleType.EARLY_ENDOSOME:
                updated_vesicles[vesicle_name] = self._process_early_endosome(vesicle, dt)
                
            elif vesicle.vesicle_type == VesicleType.LATE_ENDOSOME:
                updated_vesicles[vesicle_name] = self._process_late_endosome(vesicle, dt)
                
            elif vesicle.vesicle_type == VesicleType.LYSOSOME:
                updated_vesicles[vesicle_name] = self._process_lysosome(vesicle, dt)
                
            else:
                updated_vesicles[vesicle_name] = vesicle  # No change for other types
        
        return updated_vesicles
    
    def _process_early_endosome(self, vesicle: VesicleCompartment,
                              dt: float) -> VesicleCompartment:
        """Process early endosome maturation and recycling."""
        # Maturation to late endosome
        maturation_rate = self.params['early_to_late_endosome_rate'] * dt
        maturing_vesicles = min(vesicle.count, maturation_rate * vesicle.count)
        
        # Fast recycling (direct to surface)
        fast_recycling_rate = self.params['fast_recycling_rate'] * dt
        recycling_vesicles = min(vesicle.count - maturing_vesicles,
                               fast_recycling_rate * vesicle.count)
        
        # pH acidification
        current_pH = vesicle.pH
        target_pH = self.params['target_early_endosome_pH']
        pH_change = self.params['pH_acidification_rate'] * (target_pH - current_pH) * dt
        new_pH = max(target_pH, current_pH + pH_change)
        
        # Update vesicle
        updated_vesicle = VesicleCompartment(
            vesicle_type=vesicle.vesicle_type,
            count=max(0.0, vesicle.count - maturing_vesicles - recycling_vesicles),
            volume=vesicle.volume,
            pH=new_pH,
            cargo_concentration=vesicle.cargo_concentration.copy(),
            membrane_proteins=vesicle.membrane_proteins.copy()
        )
        
        return updated_vesicle
    
    def _process_late_endosome(self, vesicle: VesicleCompartment,
                             dt: float) -> VesicleCompartment:
        """Process late endosome lysosomal fusion."""
        # Fusion with lysosomes
        fusion_rate = self.params['endosome_lysosome_fusion_rate'] * dt
        fusing_vesicles = min(vesicle.count, fusion_rate * vesicle.count)
        
        # pH acidification
        current_pH = vesicle.pH
        target_pH = self.params['target_late_endosome_pH']
        pH_change = self.params['pH_acidification_rate'] * (target_pH - current_pH) * dt
        new_pH = max(target_pH, current_pH + pH_change)
        
        # Cargo degradation begins
        degraded_cargo = {}
        for cargo_name, concentration in vesicle.cargo_concentration.items():
            degradation_rate = 0.1 * dt  # Slow degradation in late endosomes
            remaining_concentration = concentration * (1.0 - degradation_rate)
            degraded_cargo[cargo_name] = remaining_concentration
        
        # Update vesicle
        updated_vesicle = VesicleCompartment(
            vesicle_type=vesicle.vesicle_type,
            count=max(0.0, vesicle.count - fusing_vesicles),
            volume=vesicle.volume,
            pH=new_pH,
            cargo_concentration=degraded_cargo,
            membrane_proteins=vesicle.membrane_proteins.copy()
        )
        
        return updated_vesicle
    
    def _process_lysosome(self, vesicle: VesicleCompartment,
                        dt: float) -> VesicleCompartment:
        """Process lysosomal degradation."""
        # Cargo degradation
        degraded_cargo = {}
        degradation_rate = self.params['lysosomal_degradation_rate'] * dt
        
        for cargo_name, concentration in vesicle.cargo_concentration.items():
            remaining_concentration = concentration * (1.0 - degradation_rate)
            degraded_cargo[cargo_name] = max(0.0, remaining_concentration)
        
        # pH maintenance
        target_pH = self.params['target_lysosome_pH']
        
        # Update vesicle
        updated_vesicle = VesicleCompartment(
            vesicle_type=vesicle.vesicle_type,
            count=vesicle.count,
            volume=vesicle.volume,
            pH=target_pH,
            cargo_concentration=degraded_cargo,
            membrane_proteins=vesicle.membrane_proteins.copy()
        )
        
        return updated_vesicle


class MembraneReceptorModel:
    """
    Comprehensive membrane receptor model including synthesis, trafficking,
    and degradation pathways.
    """
    
    def __init__(self, receptor_types: List[str], parameters: Optional[Dict] = None):
        """
        Initialize membrane receptor model.
        
        Args:
            receptor_types: List of receptor type names
            parameters: Model parameters
        """
        self.receptor_types = receptor_types
        
        default_params = {
            'synthesis_rate': 10.0,        # receptors per hour
            'degradation_rate': 0.1,       # per hour
            'surface_delivery_rate': 2.0,  # per hour
            'internalization_rate': 0.2,   # per hour
            'recycling_rate': 0.15,        # per hour
            'lysosomal_targeting_rate': 0.05, # per hour
            'activation_upregulation': 2.0, # fold increase
            'stress_downregulation': 0.5   # fold decrease
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
        
        # Initialize receptors
        self.receptors = {}
        for receptor_type in receptor_types:
            self.receptors[receptor_type] = MembraneReceptor(
                name=receptor_type,
                synthesis_rate=self.params['synthesis_rate'],
                degradation_rate=self.params['degradation_rate']
            )
    
    def update_receptor_dynamics(self, activation_level: float, stress_level: float,
                               ligand_concentrations: Dict[str, float],
                               dt: float) -> Dict[str, MembraneReceptor]:
        """
        Update receptor dynamics for all receptor types.
        
        Args:
            activation_level: Cell activation level (0-1)
            stress_level: Cellular stress level (0-1)
            ligand_concentrations: Ligand concentrations
            dt: Time step (hours)
            
        Returns:
            Updated receptor states
        """
        updated_receptors = {}
        
        for receptor_name, receptor in self.receptors.items():
            updated_receptor = self._update_single_receptor(
                receptor, activation_level, stress_level,
                ligand_concentrations.get(receptor_name, 0.0), dt)
            
            updated_receptors[receptor_name] = updated_receptor
        
        self.receptors = updated_receptors
        return self.receptors.copy()
    
    def _update_single_receptor(self, receptor: MembraneReceptor,
                              activation_level: float, stress_level: float,
                              ligand_concentration: float,
                              dt: float) -> MembraneReceptor:
        """Update dynamics for a single receptor type."""
        # Synthesis rate modulation
        synthesis_factor = (1.0 + self.params['activation_upregulation'] * activation_level) * \
                          (1.0 - self.params['stress_downregulation'] * stress_level)
        
        effective_synthesis_rate = receptor.synthesis_rate * synthesis_factor
        
        # Ligand-induced internalization
        ligand_factor = 1.0 + 5.0 * ligand_concentration / (0.1 + ligand_concentration)
        effective_internalization_rate = receptor.endocytosis_rate * ligand_factor
        
        # Calculate fluxes
        synthesis_flux = effective_synthesis_rate * dt
        surface_delivery_flux = self.params['surface_delivery_rate'] * receptor.internalized_count * dt
        internalization_flux = effective_internalization_rate * receptor.surface_count * dt
        recycling_flux = receptor.recycling_rate * receptor.internalized_count * dt
        degradation_flux = (self.params['lysosomal_targeting_rate'] * 
                          receptor.internalized_count * dt)
        surface_degradation_flux = receptor.degradation_rate * receptor.surface_count * dt
        
        # Update receptor counts
        new_surface_count = max(0.0, 
            receptor.surface_count + synthesis_flux + recycling_flux - 
            internalization_flux - surface_degradation_flux)
        
        new_internalized_count = max(0.0,
            receptor.internalized_count + internalization_flux - 
            recycling_flux - degradation_flux - surface_delivery_flux)
        
        # Create updated receptor
        updated_receptor = MembraneReceptor(
            name=receptor.name,
            surface_count=new_surface_count,
            internalized_count=new_internalized_count,
            synthesis_rate=receptor.synthesis_rate,
            degradation_rate=receptor.degradation_rate,
            endocytosis_rate=receptor.endocytosis_rate,
            recycling_rate=receptor.recycling_rate,
            ligand_affinity=receptor.ligand_affinity,
            clustering_factor=receptor.clustering_factor
        )
        
        return updated_receptor


class VesicleTransportModel:
    """
    Model for vesicle transport along cytoskeletal networks including
    motor protein-driven transport and diffusion.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize vesicle transport model."""
        default_params = {
            # Motor protein transport
            'kinesin_velocity': 0.8,        # μm/s
            'dynein_velocity': 1.2,         # μm/s
            'motor_processivity': 1.0,      # μm average run length
            'motor_attachment_rate': 5.0,   # per second
            'motor_detachment_rate': 1.0,   # per second
            
            # Cytoskeletal network
            'microtubule_density': 0.1,     # μm/μm³
            'actin_density': 0.5,           # μm/μm³
            'network_organization': 0.8,    # 0-1, higher = more organized
            
            # Vesicle properties
            'vesicle_diffusion_coeff': 0.01, # μm²/s
            'vesicle_radius': 0.05,         # μm
            'cargo_load_factor': 1.0,       # effect on transport
            
            # Cellular geometry
            'cell_radius': 10.0,            # μm
            'nuclear_radius': 3.0,          # μm
            'perinuclear_region': 5.0       # μm from nucleus
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def calculate_transport_velocity(self, vesicle_position: np.ndarray,
                                   target_position: np.ndarray,
                                   vesicle_type: VesicleType,
                                   cargo_load: float = 1.0) -> np.ndarray:
        """
        Calculate vesicle transport velocity.
        
        Args:
            vesicle_position: Current vesicle position (μm)
            target_position: Target position (μm)
            vesicle_type: Type of vesicle
            cargo_load: Cargo load factor
            
        Returns:
            Transport velocity vector (μm/s)
        """
        # Direction vector
        direction = target_position - vesicle_position
        distance = np.linalg.norm(direction)
        
        if distance < 0.1:  # Very close to target
            return np.zeros_like(direction)
        
        unit_direction = direction / distance
        
        # Determine transport mode based on vesicle type and position
        distance_from_nucleus = np.linalg.norm(vesicle_position)
        
        if vesicle_type in [VesicleType.EARLY_ENDOSOME, VesicleType.LATE_ENDOSOME]:
            # Endosomes: generally move toward nucleus (retrograde)
            if distance_from_nucleus > self.params['perinuclear_region']:
                motor_velocity = self.params['dynein_velocity']  # Retrograde
            else:
                motor_velocity = 0.2 * self.params['kinesin_velocity']  # Slow anterograde
        
        elif vesicle_type == VesicleType.SECRETORY_VESICLE:
            # Secretory vesicles: move toward cell periphery (anterograde)
            motor_velocity = self.params['kinesin_velocity']
        
        else:
            # Default: mixed transport
            motor_velocity = 0.5 * (self.params['kinesin_velocity'] + 
                                  self.params['dynein_velocity'])
        
        # Apply cargo load factor
        effective_velocity = motor_velocity / (1.0 + 0.5 * (cargo_load - 1.0))
        
        # Add diffusion component
        diffusion_velocity = np.random.normal(0, 
            np.sqrt(2 * self.params['vesicle_diffusion_coeff']), 2)
        
        # Combine directed and random components
        total_velocity = effective_velocity * unit_direction + diffusion_velocity
        
        return total_velocity
    
    def update_vesicle_positions(self, vesicles: Dict[str, VesicleCompartment],
                               vesicle_positions: Dict[str, np.ndarray],
                               target_positions: Dict[str, np.ndarray],
                               dt: float) -> Dict[str, np.ndarray]:
        """
        Update positions of all vesicles.
        
        Args:
            vesicles: Vesicle compartments
            vesicle_positions: Current vesicle positions
            target_positions: Target positions for each vesicle type
            dt: Time step (seconds)
            
        Returns:
            Updated vesicle positions
        """
        updated_positions = {}
        
        for vesicle_name, vesicle in vesicles.items():
            if vesicle_name in vesicle_positions:
                current_pos = vesicle_positions[vesicle_name]
                target_pos = target_positions.get(vesicle_name, current_pos)
                
                # Calculate cargo load
                total_cargo = sum(vesicle.cargo_concentration.values())
                cargo_load = 1.0 + 0.1 * total_cargo
                
                # Calculate velocity
                velocity = self.calculate_transport_velocity(
                    current_pos, target_pos, vesicle.vesicle_type, cargo_load)
                
                # Update position
                new_position = current_pos + velocity * dt
                
                # Boundary conditions (keep within cell)
                distance_from_center = np.linalg.norm(new_position)
                if distance_from_center > self.params['cell_radius']:
                    # Reflect off cell boundary
                    new_position = (new_position / distance_from_center * 
                                  self.params['cell_radius'] * 0.95)
                
                updated_positions[vesicle_name] = new_position
            else:
                # Initialize position if not present
                updated_positions[vesicle_name] = np.random.uniform(
                    -self.params['cell_radius']/2, self.params['cell_radius']/2, 2)
        
        return updated_positions


class AutophagyModel:
    """
    Autophagy pathway model including autophagosome formation,
    maturation, and lysosomal fusion.
    """
    
    def __init__(self, parameters: Optional[Dict] = None):
        """Initialize autophagy model."""
        default_params = {
            # Autophagy induction
            'basal_autophagy_rate': 0.02,      # per hour
            'stress_induction_factor': 5.0,    # fold increase under stress
            'nutrient_inhibition_factor': 0.3, # inhibition under nutrient abundance
            'mTOR_inhibition_threshold': 0.5,  # mTOR activity threshold
            
            # Autophagosome formation
            'nucleation_rate': 0.5,            # per hour
            'elongation_rate': 2.0,            # per hour
            'closure_rate': 1.0,               # per hour
            'autophagosome_capacity': 5.0,     # femtoliters
            
            # Cargo selection
            'bulk_autophagy_fraction': 0.7,    # fraction of bulk cytoplasm
            'selective_autophagy_fraction': 0.3, # fraction of specific targets
            'damaged_organelle_targeting': 3.0, # preference for damaged organelles
            
            # Autolysosome formation
            'autophagosome_lysosome_fusion_rate': 0.8, # per hour
            'autolysosome_degradation_rate': 1.5,      # per hour
        }
        
        self.params = default_params
        if parameters:
            self.params.update(parameters)
    
    def update_autophagy(self, stress_level: float, nutrient_level: float,
                        mTOR_activity: float, damaged_organelles: float,
                        dt: float) -> Dict[str, VesicleCompartment]:
        """
        Update autophagy pathway.
        
        Args:
            stress_level: Cellular stress level (0-1)
            nutrient_level: Nutrient availability (0-1)
            mTOR_activity: mTOR pathway activity (0-1)
            damaged_organelles: Level of damaged organelles (0-1)
            dt: Time step (hours)
            
        Returns:
            Autophagy-related vesicle compartments
        """
        # Calculate autophagy induction
        autophagy_rate = self._calculate_autophagy_induction_rate(
            stress_level, nutrient_level, mTOR_activity)
        
        # Autophagosome formation
        autophagosome_formation = autophagy_rate * self.params['nucleation_rate'] * dt
        
        # Cargo selection and packaging
        bulk_cargo = self.params['bulk_autophagy_fraction'] * autophagosome_formation
        selective_cargo = (self.params['selective_autophagy_fraction'] * 
                         autophagosome_formation * 
                         (1.0 + self.params['damaged_organelle_targeting'] * damaged_organelles))
        
        # Create autophagosome
        autophagosome = VesicleCompartment(
            vesicle_type=VesicleType.AUTOPHAGOSOME,
            count=autophagosome_formation,
            volume=self.params['autophagosome_capacity'],
            pH=7.0,  # Initially neutral
            cargo_concentration={
                'bulk_cytoplasm': bulk_cargo,
                'damaged_organelles': selective_cargo,
                'protein_aggregates': selective_cargo * 0.5
            }
        )
        
        # Autolysosome formation (fusion with lysosomes)
        fusion_rate = self.params['autophagosome_lysosome_fusion_rate'] * dt
        autolysosome_formation = autophagosome.count * fusion_rate
        
        autolysosome = VesicleCompartment(
            vesicle_type=VesicleType.AUTOLYSOSOME,
            count=autolysosome_formation,
            volume=autophagosome.volume,
            pH=4.5,  # Acidic after lysosomal fusion
            cargo_concentration=autophagosome.cargo_concentration.copy()
        )
        
        # Cargo degradation in autolysosomes
        degradation_rate = self.params['autolysosome_degradation_rate'] * dt
        degraded_cargo = {}
        for cargo_type, amount in autolysosome.cargo_concentration.items():
            remaining_amount = amount * (1.0 - degradation_rate)
            degraded_cargo[cargo_type] = max(0.0, remaining_amount)
        
        autolysosome.cargo_concentration = degraded_cargo
        
        return {
            'autophagosome': autophagosome,
            'autolysosome': autolysosome
        }
    
    def _calculate_autophagy_induction_rate(self, stress_level: float,
                                          nutrient_level: float,
                                          mTOR_activity: float) -> float:
        """Calculate autophagy induction rate based on cellular conditions."""
        # Basal autophagy
        basal_rate = self.params['basal_autophagy_rate']
        
        # Stress induction
        stress_factor = 1.0 + self.params['stress_induction_factor'] * stress_level
        
        # Nutrient inhibition
        nutrient_factor = 1.0 / (1.0 + self.params['nutrient_inhibition_factor'] * nutrient_level)
        
        # mTOR inhibition
        if mTOR_activity < self.params['mTOR_inhibition_threshold']:
            mTOR_factor = 2.0  # Enhanced autophagy when mTOR is low
        else:
            mTOR_factor = 1.0 / (1.0 + mTOR_activity)
        
        total_rate = basal_rate * stress_factor * nutrient_factor * mTOR_factor
        
        return total_rate
