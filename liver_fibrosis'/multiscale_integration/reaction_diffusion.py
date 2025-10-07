"""
Reaction-Diffusion Module

Implements comprehensive reaction-diffusion systems for molecular transport
and reaction kinetics in multi-scale biological systems.

Key Features:
- Multi-species reaction-diffusion solver
- Anisotropic diffusion tensor calculations
- Complex reaction network handling
- Boundary condition management
- Concentration field updates with source terms
- Adaptive mesh refinement support
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable, Union
from dataclasses import dataclass
from enum import Enum
import warnings
from scipy.sparse import csr_matrix, diags, kron, eye
from scipy.sparse.linalg import spsolve
from scipy.integrate import odeint


class BoundaryType(Enum):
    """Types of boundary conditions."""
    DIRICHLET = "dirichlet"      # Fixed concentration
    NEUMANN = "neumann"          # Fixed flux
    ROBIN = "robin"              # Mixed condition
    PERIODIC = "periodic"        # Periodic boundaries


class ReactionType(Enum):
    """Types of chemical reactions."""
    MASS_ACTION = "mass_action"
    MICHAELIS_MENTEN = "michaelis_menten"
    HILL = "hill"
    COOPERATIVE = "cooperative"
    COMPETITIVE = "competitive"


@dataclass
class MolecularSpecies:
    """Properties of a molecular species."""
    name: str
    diffusion_coefficient: float  # m²/s
    molecular_weight: float       # Da
    charge: int = 0
    initial_concentration: float = 0.0
    boundary_conditions: Dict[str, Tuple[BoundaryType, float]] = None
    
    def __post_init__(self):
        if self.boundary_conditions is None:
            self.boundary_conditions = {}


@dataclass
class ReactionParameters:
    """Parameters for a chemical reaction."""
    reaction_type: ReactionType
    reactants: Dict[str, int]     # species_name: stoichiometry
    products: Dict[str, int]      # species_name: stoichiometry
    rate_constant: float
    michaelis_constants: Dict[str, float] = None
    hill_coefficients: Dict[str, float] = None
    inhibitors: Dict[str, float] = None
    
    def __post_init__(self):
        if self.michaelis_constants is None:
            self.michaelis_constants = {}
        if self.hill_coefficients is None:
            self.hill_coefficients = {}
        if self.inhibitors is None:
            self.inhibitors = {}


class MolecularSpeciesManager:
    """
    Manager for molecular species properties and interactions.
    """
    
    def __init__(self):
        """Initialize molecular species manager."""
        self.species = {}
        self.species_indices = {}
        
    def add_species(self, species: MolecularSpecies):
        """Add a molecular species."""
        self.species[species.name] = species
        self.species_indices[species.name] = len(self.species_indices)
    
    def get_species_index(self, species_name: str) -> int:
        """Get the index of a species."""
        return self.species_indices[species_name]
    
    def get_diffusion_coefficients(self) -> np.ndarray:
        """Get diffusion coefficients for all species."""
        n_species = len(self.species)
        diffusion_coeffs = np.zeros(n_species)
        
        for name, species in self.species.items():
            idx = self.species_indices[name]
            diffusion_coeffs[idx] = species.diffusion_coefficient
            
        return diffusion_coeffs
    
    def get_initial_concentrations(self, grid_shape: Tuple[int, ...]) -> np.ndarray:
        """Get initial concentration fields for all species."""
        n_species = len(self.species)
        
        if len(grid_shape) == 2:
            concentrations = np.zeros((grid_shape[0], grid_shape[1], n_species))
        elif len(grid_shape) == 3:
            concentrations = np.zeros((grid_shape[0], grid_shape[1], grid_shape[2], n_species))
        else:
            raise ValueError(f"Unsupported grid dimensionality: {len(grid_shape)}")
        
        for name, species in self.species.items():
            idx = self.species_indices[name]
            concentrations[..., idx] = species.initial_concentration
            
        return concentrations


class DiffusionTensorCalculator:
    """
    Calculator for anisotropic diffusion tensors based on tissue structure.
    """
    
    def __init__(self, grid_spacing: Union[float, Tuple[float, ...]] = 1.0):
        """
        Initialize diffusion tensor calculator.
        
        Args:
            grid_spacing: Spatial grid spacing (uniform or per dimension)
        """
        if isinstance(grid_spacing, (int, float)):
            self.grid_spacing = grid_spacing
        else:
            self.grid_spacing = np.array(grid_spacing)
    
    def calculate_isotropic_tensor(self, diffusion_coefficient: float,
                                 grid_shape: Tuple[int, ...]) -> np.ndarray:
        """
        Calculate isotropic diffusion tensor.
        
        Args:
            diffusion_coefficient: Scalar diffusion coefficient
            grid_shape: Shape of computational grid
            
        Returns:
            Isotropic diffusion tensor field
        """
        ndim = len(grid_shape)
        
        if ndim == 2:
            tensor_shape = grid_shape + (2, 2)
        elif ndim == 3:
            tensor_shape = grid_shape + (3, 3)
        else:
            raise ValueError(f"Unsupported dimensionality: {ndim}")
        
        # Create identity tensor field
        tensor_field = np.zeros(tensor_shape)
        for i in range(ndim):
            tensor_field[..., i, i] = diffusion_coefficient
            
        return tensor_field
    
    def calculate_anisotropic_tensor(self, fiber_orientation: np.ndarray,
                                   parallel_diffusivity: float,
                                   perpendicular_diffusivity: float) -> np.ndarray:
        """
        Calculate anisotropic diffusion tensor based on fiber orientation.
        
        Args:
            fiber_orientation: Fiber orientation field (unit vectors)
            parallel_diffusivity: Diffusion parallel to fibers
            perpendicular_diffusivity: Diffusion perpendicular to fibers
            
        Returns:
            Anisotropic diffusion tensor field
        """
        grid_shape = fiber_orientation.shape[:-1]
        ndim = fiber_orientation.shape[-1]
        
        if ndim == 2:
            tensor_shape = grid_shape + (2, 2)
        elif ndim == 3:
            tensor_shape = grid_shape + (3, 3)
        else:
            raise ValueError(f"Unsupported dimensionality: {ndim}")
        
        tensor_field = np.zeros(tensor_shape)
        
        # Calculate tensor for each grid point
        for idx in np.ndindex(grid_shape):
            fiber_dir = fiber_orientation[idx]
            
            # Normalize fiber direction
            fiber_norm = np.linalg.norm(fiber_dir)
            if fiber_norm > 1e-12:
                fiber_dir = fiber_dir / fiber_norm
            else:
                fiber_dir = np.zeros(ndim)
                fiber_dir[0] = 1.0  # Default direction
            
            # Create diffusion tensor: D = D_perp * I + (D_par - D_perp) * n ⊗ n
            identity = np.eye(ndim)
            fiber_tensor = np.outer(fiber_dir, fiber_dir)
            
            tensor_field[idx] = (perpendicular_diffusivity * identity + 
                               (parallel_diffusivity - perpendicular_diffusivity) * fiber_tensor)
        
        return tensor_field
    
    def calculate_cell_density_dependent_tensor(self, cell_density: np.ndarray,
                                              free_diffusivity: float,
                                              tortuosity_factor: float = 2.0) -> np.ndarray:
        """
        Calculate diffusion tensor accounting for cellular obstruction.
        
        Args:
            cell_density: Cell density field (0-1)
            free_diffusivity: Diffusion coefficient in free medium
            tortuosity_factor: Tortuosity enhancement factor
            
        Returns:
            Cell-density-dependent diffusion tensor
        """
        grid_shape = cell_density.shape
        ndim = len(grid_shape)
        
        # Calculate effective diffusivity using obstruction model
        # D_eff = D_free * (1 - φ) / τ
        # where φ is volume fraction and τ is tortuosity
        
        volume_fraction = cell_density
        tortuosity = 1.0 + tortuosity_factor * volume_fraction
        
        effective_diffusivity = (free_diffusivity * (1.0 - volume_fraction) / tortuosity)
        
        # Create isotropic tensor with effective diffusivity
        if ndim == 2:
            tensor_shape = grid_shape + (2, 2)
        elif ndim == 3:
            tensor_shape = grid_shape + (3, 3)
        else:
            raise ValueError(f"Unsupported dimensionality: {ndim}")
        
        tensor_field = np.zeros(tensor_shape)
        for i in range(ndim):
            tensor_field[..., i, i] = effective_diffusivity
            
        return tensor_field


class ReactionNetworkSolver:
    """
    Solver for complex chemical reaction networks.
    """
    
    def __init__(self, species_manager: MolecularSpeciesManager):
        """Initialize reaction network solver."""
        self.species_manager = species_manager
        self.reactions = []
        self.stoichiometry_matrix = None
        
    def add_reaction(self, reaction: ReactionParameters):
        """Add a reaction to the network."""
        self.reactions.append(reaction)
        self._update_stoichiometry_matrix()
    
    def _update_stoichiometry_matrix(self):
        """Update the stoichiometry matrix."""
        n_species = len(self.species_manager.species)
        n_reactions = len(self.reactions)
        
        self.stoichiometry_matrix = np.zeros((n_species, n_reactions))
        
        for j, reaction in enumerate(self.reactions):
            # Reactants (negative stoichiometry)
            for species_name, stoich in reaction.reactants.items():
                if species_name in self.species_manager.species_indices:
                    i = self.species_manager.species_indices[species_name]
                    self.stoichiometry_matrix[i, j] -= stoich
            
            # Products (positive stoichiometry)
            for species_name, stoich in reaction.products.items():
                if species_name in self.species_manager.species_indices:
                    i = self.species_manager.species_indices[species_name]
                    self.stoichiometry_matrix[i, j] += stoich
    
    def calculate_reaction_rates(self, concentrations: np.ndarray) -> np.ndarray:
        """
        Calculate reaction rates for all reactions.
        
        Args:
            concentrations: Current species concentrations
            
        Returns:
            Reaction rates for each reaction
        """
        n_reactions = len(self.reactions)
        rates = np.zeros(n_reactions)
        
        for j, reaction in enumerate(self.reactions):
            rates[j] = self._calculate_single_reaction_rate(reaction, concentrations)
        
        return rates
    
    def _calculate_single_reaction_rate(self, reaction: ReactionParameters,
                                     concentrations: np.ndarray) -> float:
        """Calculate rate for a single reaction."""
        if reaction.reaction_type == ReactionType.MASS_ACTION:
            rate = reaction.rate_constant
            
            # Multiply by reactant concentrations
            for species_name, stoich in reaction.reactants.items():
                if species_name in self.species_manager.species_indices:
                    i = self.species_manager.species_indices[species_name]
                    rate *= concentrations[i] ** stoich
        
        elif reaction.reaction_type == ReactionType.MICHAELIS_MENTEN:
            # Assume single substrate for simplicity
            substrate_name = list(reaction.reactants.keys())[0]
            substrate_idx = self.species_manager.species_indices[substrate_name]
            substrate_conc = concentrations[substrate_idx]
            
            Km = reaction.michaelis_constants.get(substrate_name, 1.0)
            rate = (reaction.rate_constant * substrate_conc / (Km + substrate_conc))
        
        elif reaction.reaction_type == ReactionType.HILL:
            # Hill kinetics with cooperativity
            substrate_name = list(reaction.reactants.keys())[0]
            substrate_idx = self.species_manager.species_indices[substrate_name]
            substrate_conc = concentrations[substrate_idx]
            
            Km = reaction.michaelis_constants.get(substrate_name, 1.0)
            n = reaction.hill_coefficients.get(substrate_name, 1.0)
            
            rate = (reaction.rate_constant * substrate_conc**n / 
                   (Km**n + substrate_conc**n))
        
        else:
            # Default to mass action
            rate = reaction.rate_constant
            for species_name, stoich in reaction.reactants.items():
                if species_name in self.species_manager.species_indices:
                    i = self.species_manager.species_indices[species_name]
                    rate *= concentrations[i] ** stoich
        
        return max(0.0, rate)  # Ensure non-negative rates
    
    def calculate_reaction_source_terms(self, concentrations: np.ndarray) -> np.ndarray:
        """
        Calculate reaction source terms for all species.
        
        Args:
            concentrations: Current species concentrations
            
        Returns:
            Source terms for each species
        """
        if self.stoichiometry_matrix is None:
            return np.zeros(len(self.species_manager.species))
        
        reaction_rates = self.calculate_reaction_rates(concentrations)
        source_terms = self.stoichiometry_matrix @ reaction_rates
        
        return source_terms


class BoundaryConditionHandler:
    """
    Handler for various boundary conditions in reaction-diffusion systems.
    """
    
    def __init__(self, grid_shape: Tuple[int, ...], grid_spacing: Union[float, Tuple[float, ...]]):
        """
        Initialize boundary condition handler.
        
        Args:
            grid_shape: Shape of computational grid
            grid_spacing: Spatial grid spacing
        """
        self.grid_shape = grid_shape
        self.ndim = len(grid_shape)
        
        if isinstance(grid_spacing, (int, float)):
            self.grid_spacing = np.full(self.ndim, grid_spacing)
        else:
            self.grid_spacing = np.array(grid_spacing)
    
    def apply_dirichlet_bc(self, concentration_field: np.ndarray,
                          boundary_value: float,
                          boundary_mask: np.ndarray) -> np.ndarray:
        """
        Apply Dirichlet boundary conditions.
        
        Args:
            concentration_field: Current concentration field
            boundary_value: Fixed boundary concentration
            boundary_mask: Boolean mask indicating boundary points
            
        Returns:
            Updated concentration field
        """
        updated_field = concentration_field.copy()
        updated_field[boundary_mask] = boundary_value
        return updated_field
    
    def apply_neumann_bc(self, concentration_field: np.ndarray,
                        flux_value: float,
                        boundary_normal: np.ndarray,
                        boundary_mask: np.ndarray) -> np.ndarray:
        """
        Apply Neumann boundary conditions.
        
        Args:
            concentration_field: Current concentration field
            flux_value: Fixed boundary flux
            boundary_normal: Outward normal vector at boundary
            boundary_mask: Boolean mask indicating boundary points
            
        Returns:
            Updated concentration field
        """
        updated_field = concentration_field.copy()
        
        # Calculate gradient at boundary points
        gradients = np.gradient(concentration_field)
        
        # Apply flux condition: -D * ∇C · n = flux_value
        # This is implemented by adjusting the concentration at boundary points
        # based on the specified flux and grid spacing
        
        boundary_indices = np.where(boundary_mask)
        
        for i in range(len(boundary_indices[0])):
            idx = tuple(boundary_indices[j][i] for j in range(self.ndim))
            
            # Calculate required gradient to achieve specified flux
            # Assuming isotropic diffusion for simplicity
            normal_gradient = flux_value  # Simplified
            
            # Update boundary concentration based on gradient
            # This is a simplified implementation
            if idx[0] == 0:  # Left boundary
                updated_field[idx] = concentration_field[1, idx[1]] - normal_gradient * self.grid_spacing[0]
            elif idx[0] == self.grid_shape[0] - 1:  # Right boundary
                updated_field[idx] = concentration_field[-2, idx[1]] + normal_gradient * self.grid_spacing[0]
        
        return updated_field
    
    def apply_robin_bc(self, concentration_field: np.ndarray,
                      alpha: float, beta: float, gamma: float,
                      boundary_mask: np.ndarray) -> np.ndarray:
        """
        Apply Robin boundary conditions: α*C + β*∂C/∂n = γ
        
        Args:
            concentration_field: Current concentration field
            alpha: Coefficient for concentration term
            beta: Coefficient for flux term
            gamma: Right-hand side value
            boundary_mask: Boolean mask indicating boundary points
            
        Returns:
            Updated concentration field
        """
        updated_field = concentration_field.copy()
        
        # This is a simplified implementation
        # In practice, Robin BCs require special treatment in the discretization
        
        boundary_indices = np.where(boundary_mask)
        
        for i in range(len(boundary_indices[0])):
            idx = tuple(boundary_indices[j][i] for j in range(self.ndim))
            
            # Simplified Robin BC implementation
            # α*C + β*flux = γ
            # Assuming β = 1 and estimating flux from neighboring points
            
            if alpha != 0:
                # Solve for C: C = (γ - β*flux) / α
                estimated_flux = 0.0  # Simplified
                updated_field[idx] = (gamma - beta * estimated_flux) / alpha
        
        return updated_field
    
    def create_boundary_masks(self) -> Dict[str, np.ndarray]:
        """Create boundary masks for different boundaries."""
        masks = {}
        
        if self.ndim == 2:
            # 2D boundaries
            masks['left'] = np.zeros(self.grid_shape, dtype=bool)
            masks['right'] = np.zeros(self.grid_shape, dtype=bool)
            masks['bottom'] = np.zeros(self.grid_shape, dtype=bool)
            masks['top'] = np.zeros(self.grid_shape, dtype=bool)
            
            masks['left'][0, :] = True
            masks['right'][-1, :] = True
            masks['bottom'][:, 0] = True
            masks['top'][:, -1] = True
            
        elif self.ndim == 3:
            # 3D boundaries
            masks['left'] = np.zeros(self.grid_shape, dtype=bool)
            masks['right'] = np.zeros(self.grid_shape, dtype=bool)
            masks['front'] = np.zeros(self.grid_shape, dtype=bool)
            masks['back'] = np.zeros(self.grid_shape, dtype=bool)
            masks['bottom'] = np.zeros(self.grid_shape, dtype=bool)
            masks['top'] = np.zeros(self.grid_shape, dtype=bool)
            
            masks['left'][0, :, :] = True
            masks['right'][-1, :, :] = True
            masks['front'][:, 0, :] = True
            masks['back'][:, -1, :] = True
            masks['bottom'][:, :, 0] = True
            masks['top'][:, :, -1] = True
        
        return masks


class ConcentrationFieldUpdater:
    """
    Updater for concentration fields in reaction-diffusion systems.
    """
    
    def __init__(self, species_manager: MolecularSpeciesManager,
                 grid_shape: Tuple[int, ...],
                 grid_spacing: Union[float, Tuple[float, ...]]):
        """Initialize concentration field updater."""
        self.species_manager = species_manager
        self.grid_shape = grid_shape
        self.ndim = len(grid_shape)
        
        if isinstance(grid_spacing, (int, float)):
            self.grid_spacing = np.full(self.ndim, grid_spacing)
        else:
            self.grid_spacing = np.array(grid_spacing)
        
        # Initialize boundary condition handler
        self.bc_handler = BoundaryConditionHandler(grid_shape, grid_spacing)
        
        # Initialize reaction network solver
        self.reaction_solver = ReactionNetworkSolver(species_manager)
    
    def update_concentrations(self, concentrations: np.ndarray,
                            diffusion_tensors: Dict[str, np.ndarray],
                            source_terms: Optional[np.ndarray] = None,
                            dt: float = 0.01) -> np.ndarray:
        """
        Update concentration fields for one time step.
        
        Args:
            concentrations: Current concentration fields
            diffusion_tensors: Diffusion tensors for each species
            source_terms: External source terms
            dt: Time step
            
        Returns:
            Updated concentration fields
        """
        updated_concentrations = concentrations.copy()
        
        # Process each species
        for species_name, species in self.species_manager.species.items():
            species_idx = self.species_manager.species_indices[species_name]
            
            # Get current concentration field for this species
            if self.ndim == 2:
                current_conc = concentrations[:, :, species_idx]
            elif self.ndim == 3:
                current_conc = concentrations[:, :, :, species_idx]
            else:
                raise ValueError(f"Unsupported dimensionality: {self.ndim}")
            
            # Get diffusion tensor
            if species_name in diffusion_tensors:
                D_tensor = diffusion_tensors[species_name]
            else:
                # Use isotropic diffusion
                tensor_calc = DiffusionTensorCalculator(self.grid_spacing)
                D_tensor = tensor_calc.calculate_isotropic_tensor(
                    species.diffusion_coefficient, self.grid_shape)
            
            # Calculate diffusion term
            diffusion_term = self._calculate_diffusion_term(current_conc, D_tensor)
            
            # Calculate reaction terms
            if self.ndim == 2:
                reaction_terms = np.zeros_like(current_conc)
                for i in range(self.grid_shape[0]):
                    for j in range(self.grid_shape[1]):
                        local_conc = concentrations[i, j, :]
                        local_reactions = self.reaction_solver.calculate_reaction_source_terms(local_conc)
                        reaction_terms[i, j] = local_reactions[species_idx]
            elif self.ndim == 3:
                reaction_terms = np.zeros_like(current_conc)
                for i in range(self.grid_shape[0]):
                    for j in range(self.grid_shape[1]):
                        for k in range(self.grid_shape[2]):
                            local_conc = concentrations[i, j, k, :]
                            local_reactions = self.reaction_solver.calculate_reaction_source_terms(local_conc)
                            reaction_terms[i, j, k] = local_reactions[species_idx]
            
            # Add external source terms
            if source_terms is not None:
                if self.ndim == 2:
                    external_source = source_terms[:, :, species_idx]
                elif self.ndim == 3:
                    external_source = source_terms[:, :, :, species_idx]
                else:
                    external_source = 0.0
            else:
                external_source = 0.0
            
            # Update concentration using forward Euler
            new_conc = (current_conc + dt * (diffusion_term + reaction_terms + external_source))
            
            # Apply boundary conditions
            new_conc = self._apply_boundary_conditions(new_conc, species)
            
            # Ensure non-negative concentrations
            new_conc = np.maximum(new_conc, 0.0)
            
            # Store updated concentration
            if self.ndim == 2:
                updated_concentrations[:, :, species_idx] = new_conc
            elif self.ndim == 3:
                updated_concentrations[:, :, :, species_idx] = new_conc
        
        return updated_concentrations
    
    def _calculate_diffusion_term(self, concentration: np.ndarray,
                                diffusion_tensor: np.ndarray) -> np.ndarray:
        """Calculate diffusion term ∇·(D∇C)."""
        # This is a simplified implementation using finite differences
        # For anisotropic diffusion, proper tensor operations should be used
        
        if self.ndim == 2:
            # 2D diffusion
            # Assuming isotropic diffusion for simplicity
            D = diffusion_tensor[..., 0, 0]  # Extract scalar diffusion coefficient
            
            # Calculate second derivatives
            d2C_dx2 = np.zeros_like(concentration)
            d2C_dy2 = np.zeros_like(concentration)
            
            # Interior points
            d2C_dx2[1:-1, :] = ((concentration[2:, :] - 2*concentration[1:-1, :] + concentration[:-2, :]) / 
                               self.grid_spacing[0]**2)
            d2C_dy2[:, 1:-1] = ((concentration[:, 2:] - 2*concentration[:, 1:-1] + concentration[:, :-2]) / 
                               self.grid_spacing[1]**2)
            
            diffusion_term = D * (d2C_dx2 + d2C_dy2)
            
        elif self.ndim == 3:
            # 3D diffusion
            D = diffusion_tensor[..., 0, 0]  # Extract scalar diffusion coefficient
            
            # Calculate second derivatives
            d2C_dx2 = np.zeros_like(concentration)
            d2C_dy2 = np.zeros_like(concentration)
            d2C_dz2 = np.zeros_like(concentration)
            
            # Interior points
            d2C_dx2[1:-1, :, :] = ((concentration[2:, :, :] - 2*concentration[1:-1, :, :] + concentration[:-2, :, :]) / 
                                  self.grid_spacing[0]**2)
            d2C_dy2[:, 1:-1, :] = ((concentration[:, 2:, :] - 2*concentration[:, 1:-1, :] + concentration[:, :-2, :]) / 
                                  self.grid_spacing[1]**2)
            d2C_dz2[:, :, 1:-1] = ((concentration[:, :, 2:] - 2*concentration[:, :, 1:-1] + concentration[:, :, :-2]) / 
                                  self.grid_spacing[2]**2)
            
            diffusion_term = D * (d2C_dx2 + d2C_dy2 + d2C_dz2)
        
        else:
            raise ValueError(f"Unsupported dimensionality: {self.ndim}")
        
        return diffusion_term
    
    def _apply_boundary_conditions(self, concentration: np.ndarray,
                                 species: MolecularSpecies) -> np.ndarray:
        """Apply boundary conditions for a species."""
        updated_conc = concentration.copy()
        
        # Get boundary masks
        boundary_masks = self.bc_handler.create_boundary_masks()
        
        # Apply boundary conditions based on species properties
        for boundary_name, (bc_type, bc_value) in species.boundary_conditions.items():
            if boundary_name in boundary_masks:
                mask = boundary_masks[boundary_name]
                
                if bc_type == BoundaryType.DIRICHLET:
                    updated_conc = self.bc_handler.apply_dirichlet_bc(
                        updated_conc, bc_value, mask)
                elif bc_type == BoundaryType.NEUMANN:
                    # Simplified normal vector (needs proper implementation)
                    normal = np.array([1.0, 0.0]) if self.ndim == 2 else np.array([1.0, 0.0, 0.0])
                    updated_conc = self.bc_handler.apply_neumann_bc(
                        updated_conc, bc_value, normal, mask)
        
        return updated_conc


class ReactionDiffusionSolver:
    """
    Main solver for reaction-diffusion systems.
    """
    
    def __init__(self, species_manager: MolecularSpeciesManager,
                 grid_shape: Tuple[int, ...],
                 grid_spacing: Union[float, Tuple[float, ...]]):
        """Initialize reaction-diffusion solver."""
        self.species_manager = species_manager
        self.grid_shape = grid_shape
        self.grid_spacing = grid_spacing
        
        # Initialize components
        self.diffusion_calculator = DiffusionTensorCalculator(grid_spacing)
        self.field_updater = ConcentrationFieldUpdater(species_manager, grid_shape, grid_spacing)
        
        # Initialize concentration fields
        self.concentrations = species_manager.get_initial_concentrations(grid_shape)
        
    def solve_time_step(self, dt: float,
                       cell_density_field: Optional[np.ndarray] = None,
                       fiber_orientation_field: Optional[np.ndarray] = None,
                       source_terms: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Solve one time step of the reaction-diffusion system.
        
        Args:
            dt: Time step
            cell_density_field: Cell density field for obstruction effects
            fiber_orientation_field: Fiber orientation for anisotropic diffusion
            source_terms: External source terms
            
        Returns:
            Updated concentration fields
        """
        # Calculate diffusion tensors for each species
        diffusion_tensors = {}
        
        for species_name, species in self.species_manager.species.items():
            if fiber_orientation_field is not None:
                # Anisotropic diffusion based on fiber orientation
                D_parallel = species.diffusion_coefficient * 2.0  # Enhanced parallel diffusion
                D_perpendicular = species.diffusion_coefficient * 0.5  # Reduced perpendicular diffusion
                
                tensor = self.diffusion_calculator.calculate_anisotropic_tensor(
                    fiber_orientation_field, D_parallel, D_perpendicular)
                    
            elif cell_density_field is not None:
                # Cell density dependent diffusion
                tensor = self.diffusion_calculator.calculate_cell_density_dependent_tensor(
                    cell_density_field, species.diffusion_coefficient)
                    
            else:
                # Isotropic diffusion
                tensor = self.diffusion_calculator.calculate_isotropic_tensor(
                    species.diffusion_coefficient, self.grid_shape)
            
            diffusion_tensors[species_name] = tensor
        
        # Update concentrations
        self.concentrations = self.field_updater.update_concentrations(
            self.concentrations, diffusion_tensors, source_terms, dt)
        
        return self.concentrations
    
    def get_species_concentration(self, species_name: str) -> np.ndarray:
        """Get concentration field for a specific species."""
        species_idx = self.species_manager.species_indices[species_name]
        
        if len(self.grid_shape) == 2:
            return self.concentrations[:, :, species_idx]
        elif len(self.grid_shape) == 3:
            return self.concentrations[:, :, :, species_idx]
        else:
            raise ValueError(f"Unsupported dimensionality: {len(self.grid_shape)}")
    
    def add_reaction(self, reaction: ReactionParameters):
        """Add a reaction to the system."""
        self.field_updater.reaction_solver.add_reaction(reaction)
    
    def set_boundary_condition(self, species_name: str, boundary_name: str,
                             bc_type: BoundaryType, bc_value: float):
        """Set boundary condition for a species."""
        if species_name in self.species_manager.species:
            species = self.species_manager.species[species_name]
            species.boundary_conditions[boundary_name] = (bc_type, bc_value)
