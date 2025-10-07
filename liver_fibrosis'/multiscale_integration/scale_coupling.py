"""
Scale Coupling Module

Implements multi-scale coupling frameworks for integrating molecular, cellular,
and tissue-level processes in biological systems.

Key Features:
- Multi-scale coupling framework with scale separation
- Cross-scale information transfer mechanisms
- Homogenization techniques for upscaling
- Downscaling operators for local refinement
- Scale-dependent time stepping
- Coupling stability analysis
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable, Union, Any
from dataclasses import dataclass
from enum import Enum
import warnings
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import spsolve


class ScaleLevel(Enum):
    """Enumeration of scale levels."""
    MOLECULAR = "molecular"
    CELLULAR = "cellular"
    TISSUE = "tissue"
    ORGAN = "organ"


class CouplingType(Enum):
    """Types of scale coupling."""
    ONE_WAY = "one_way"
    TWO_WAY = "two_way"
    HIERARCHICAL = "hierarchical"
    CONCURRENT = "concurrent"


@dataclass
class ScaleProperties:
    """Properties defining a scale level."""
    level: ScaleLevel
    characteristic_length: float  # meters
    characteristic_time: float    # seconds
    characteristic_mass: float    # kg
    spatial_resolution: float     # meters
    temporal_resolution: float    # seconds
    domain_size: Tuple[float, float, float]  # meters


@dataclass
class CouplingInterface:
    """Interface between two scale levels."""
    source_scale: ScaleLevel
    target_scale: ScaleLevel
    coupling_type: CouplingType
    transfer_variables: List[str]
    coupling_strength: float = 1.0
    coupling_frequency: float = 1.0  # Hz
    interpolation_method: str = "linear"


class MultiScaleCouplingFramework:
    """
    Framework for managing multi-scale coupling between different
    scale levels in biological systems.
    """
    
    def __init__(self, scale_properties: Dict[ScaleLevel, ScaleProperties]):
        """
        Initialize multi-scale coupling framework.
        
        Args:
            scale_properties: Dictionary mapping scale levels to their properties
        """
        self.scale_properties = scale_properties
        self.coupling_interfaces = {}
        self.scale_models = {}
        self.coupling_operators = {}
        
        # Initialize scale separation handler
        self.scale_separator = ScaleSeparationHandler(scale_properties)
        
        # Initialize information transfer
        self.info_transfer = CrossScaleInformationTransfer()
        
        # Initialize homogenization techniques
        self.homogenizer = HomogenizationTechniques()
        
    def add_coupling_interface(self, interface: CouplingInterface):
        """Add a coupling interface between scale levels."""
        interface_key = (interface.source_scale, interface.target_scale)
        self.coupling_interfaces[interface_key] = interface
        
        # Create appropriate coupling operators
        if interface.source_scale.value < interface.target_scale.value:
            # Upscaling (fine to coarse)
            operator = UpscalingOperator(interface)
        else:
            # Downscaling (coarse to fine)
            operator = DownscalingOperator(interface)
            
        self.coupling_operators[interface_key] = operator
    
    def register_scale_model(self, scale: ScaleLevel, model: Any):
        """Register a model for a specific scale level."""
        self.scale_models[scale] = model
    
    def update_coupled_system(self, dt: float, 
                            scale_states: Dict[ScaleLevel, Dict[str, np.ndarray]]) -> Dict[ScaleLevel, Dict[str, np.ndarray]]:
        """
        Update the entire coupled multi-scale system.
        
        Args:
            dt: Global time step
            scale_states: Current states at each scale level
            
        Returns:
            Updated states at each scale level
        """
        updated_states = {}
        
        # Determine coupling order based on scale hierarchy
        coupling_order = self._determine_coupling_order()
        
        # Process coupling in hierarchical order
        for source_scale, target_scale in coupling_order:
            interface_key = (source_scale, target_scale)
            
            if interface_key in self.coupling_interfaces:
                interface = self.coupling_interfaces[interface_key]
                operator = self.coupling_operators[interface_key]
                
                # Check if coupling should occur at this time step
                if self._should_couple(interface, dt):
                    # Transfer information between scales
                    transferred_data = operator.transfer(
                        scale_states[source_scale],
                        scale_states.get(target_scale, {}),
                        dt
                    )
                    
                    # Update target scale state
                    if target_scale not in updated_states:
                        updated_states[target_scale] = scale_states[target_scale].copy()
                    
                    for var_name, data in transferred_data.items():
                        updated_states[target_scale][var_name] = data
        
        # Update states that weren't modified by coupling
        for scale in scale_states:
            if scale not in updated_states:
                updated_states[scale] = scale_states[scale].copy()
        
        return updated_states
    
    def _determine_coupling_order(self) -> List[Tuple[ScaleLevel, ScaleLevel]]:
        """Determine the order of scale coupling operations."""
        # Sort by scale hierarchy (molecular -> cellular -> tissue -> organ)
        scale_order = [ScaleLevel.MOLECULAR, ScaleLevel.CELLULAR, 
                      ScaleLevel.TISSUE, ScaleLevel.ORGAN]
        
        coupling_order = []
        
        # First pass: upscaling (fine to coarse)
        for i, source_scale in enumerate(scale_order[:-1]):
            for target_scale in scale_order[i+1:]:
                if (source_scale, target_scale) in self.coupling_interfaces:
                    coupling_order.append((source_scale, target_scale))
        
        # Second pass: downscaling (coarse to fine)
        for i, source_scale in enumerate(scale_order[1:], 1):
            for target_scale in scale_order[:i]:
                if (source_scale, target_scale) in self.coupling_interfaces:
                    coupling_order.append((source_scale, target_scale))
        
        return coupling_order
    
    def _should_couple(self, interface: CouplingInterface, dt: float) -> bool:
        """Determine if coupling should occur at current time step."""
        coupling_period = 1.0 / interface.coupling_frequency
        return dt >= coupling_period or np.random.rand() < dt * interface.coupling_frequency
    
    def analyze_coupling_stability(self, scale_states: Dict[ScaleLevel, Dict[str, np.ndarray]]) -> Dict[str, float]:
        """
        Analyze stability of the coupled multi-scale system.
        
        Args:
            scale_states: Current system states
            
        Returns:
            Stability metrics
        """
        stability_metrics = {}
        
        # Calculate scale separation ratios
        for (source_scale, target_scale), interface in self.coupling_interfaces.items():
            source_props = self.scale_properties[source_scale]
            target_props = self.scale_properties[target_scale]
            
            # Length scale separation
            length_ratio = target_props.characteristic_length / source_props.characteristic_length
            
            # Time scale separation
            time_ratio = target_props.characteristic_time / source_props.characteristic_time
            
            # Coupling strength analysis
            coupling_strength = interface.coupling_strength
            
            stability_metrics[f"{source_scale.value}_to_{target_scale.value}_length_ratio"] = length_ratio
            stability_metrics[f"{source_scale.value}_to_{target_scale.value}_time_ratio"] = time_ratio
            stability_metrics[f"{source_scale.value}_to_{target_scale.value}_coupling_strength"] = coupling_strength
            
            # Stability criterion (empirical)
            stability_criterion = (length_ratio > 10.0 and time_ratio > 10.0 and 
                                 coupling_strength < 1.0)
            stability_metrics[f"{source_scale.value}_to_{target_scale.value}_stable"] = float(stability_criterion)
        
        return stability_metrics


class ScaleSeparationHandler:
    """
    Handler for managing scale separation in multi-scale systems.
    """
    
    def __init__(self, scale_properties: Dict[ScaleLevel, ScaleProperties]):
        """Initialize scale separation handler."""
        self.scale_properties = scale_properties
        
    def calculate_scale_separation_parameter(self, scale1: ScaleLevel, 
                                           scale2: ScaleLevel) -> float:
        """
        Calculate scale separation parameter between two scales.
        
        Args:
            scale1: First scale level
            scale2: Second scale level
            
        Returns:
            Scale separation parameter (epsilon)
        """
        props1 = self.scale_properties[scale1]
        props2 = self.scale_properties[scale2]
        
        # Use characteristic length ratio as separation parameter
        epsilon = min(props1.characteristic_length, props2.characteristic_length) / \
                 max(props1.characteristic_length, props2.characteristic_length)
        
        return epsilon
    
    def is_scale_separated(self, scale1: ScaleLevel, scale2: ScaleLevel,
                          threshold: float = 0.1) -> bool:
        """
        Check if two scales are well-separated.
        
        Args:
            scale1: First scale level
            scale2: Second scale level
            threshold: Separation threshold
            
        Returns:
            True if scales are well-separated
        """
        epsilon = self.calculate_scale_separation_parameter(scale1, scale2)
        return epsilon < threshold
    
    def calculate_homogenization_validity(self, scale: ScaleLevel,
                                        field_data: np.ndarray,
                                        grid_spacing: float) -> float:
        """
        Calculate validity of homogenization approximation.
        
        Args:
            scale: Scale level
            field_data: Field data to analyze
            grid_spacing: Spatial grid spacing
            
        Returns:
            Homogenization validity measure (0-1)
        """
        props = self.scale_properties[scale]
        
        # Calculate field gradients
        gradients = np.gradient(field_data)
        gradient_magnitude = np.sqrt(sum(g**2 for g in gradients))
        
        # Calculate characteristic gradient scale
        characteristic_gradient = np.mean(gradient_magnitude)
        field_scale = np.mean(np.abs(field_data)) / (characteristic_gradient + 1e-12)
        
        # Homogenization validity: field varies slowly compared to characteristic scale
        validity = min(1.0, field_scale / props.characteristic_length)
        
        return validity


class CrossScaleInformationTransfer:
    """
    Manager for cross-scale information transfer mechanisms.
    """
    
    def __init__(self):
        """Initialize cross-scale information transfer."""
        self.transfer_operators = {}
        
    def register_transfer_operator(self, source_scale: ScaleLevel,
                                 target_scale: ScaleLevel,
                                 operator: Callable):
        """Register a custom transfer operator."""
        key = (source_scale, target_scale)
        self.transfer_operators[key] = operator
    
    def transfer_molecular_to_cellular(self, molecular_data: Dict[str, np.ndarray],
                                     cellular_grid: np.ndarray,
                                     cell_positions: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Transfer molecular-scale information to cellular scale.
        
        Args:
            molecular_data: Molecular concentration fields
            cellular_grid: Cellular grid coordinates
            cell_positions: Individual cell positions
            
        Returns:
            Cellular-scale data
        """
        cellular_data = {}
        
        for species_name, concentration_field in molecular_data.items():
            # Interpolate molecular concentrations to cell positions
            if len(concentration_field.shape) == 2:
                # 2D field
                x_coords = np.linspace(0, concentration_field.shape[1]-1, concentration_field.shape[1])
                y_coords = np.linspace(0, concentration_field.shape[0]-1, concentration_field.shape[0])
                
                interpolator = RegularGridInterpolator(
                    (y_coords, x_coords), concentration_field,
                    method='linear', bounds_error=False, fill_value=0.0
                )
                
                cellular_concentrations = interpolator(cell_positions)
                
            elif len(concentration_field.shape) == 3:
                # 3D field
                x_coords = np.linspace(0, concentration_field.shape[2]-1, concentration_field.shape[2])
                y_coords = np.linspace(0, concentration_field.shape[1]-1, concentration_field.shape[1])
                z_coords = np.linspace(0, concentration_field.shape[0]-1, concentration_field.shape[0])
                
                interpolator = RegularGridInterpolator(
                    (z_coords, y_coords, x_coords), concentration_field,
                    method='linear', bounds_error=False, fill_value=0.0
                )
                
                cellular_concentrations = interpolator(cell_positions)
            
            else:
                raise ValueError(f"Unsupported field dimensionality: {concentration_field.shape}")
            
            cellular_data[species_name] = cellular_concentrations
        
        return cellular_data
    
    def transfer_cellular_to_tissue(self, cellular_data: Dict[str, np.ndarray],
                                  cell_positions: np.ndarray,
                                  tissue_grid: Tuple[np.ndarray, ...]) -> Dict[str, np.ndarray]:
        """
        Transfer cellular-scale information to tissue scale.
        
        Args:
            cellular_data: Cellular-scale data
            cell_positions: Cell positions
            tissue_grid: Tissue-scale grid
            
        Returns:
            Tissue-scale fields
        """
        tissue_data = {}
        
        # Create tissue grid points
        if len(tissue_grid) == 2:
            # 2D tissue
            grid_points = np.stack(np.meshgrid(tissue_grid[0], tissue_grid[1], indexing='ij'), axis=-1)
            grid_shape = grid_points.shape[:-1]
            grid_points_flat = grid_points.reshape(-1, 2)
        elif len(tissue_grid) == 3:
            # 3D tissue
            grid_points = np.stack(np.meshgrid(tissue_grid[0], tissue_grid[1], tissue_grid[2], indexing='ij'), axis=-1)
            grid_shape = grid_points.shape[:-1]
            grid_points_flat = grid_points.reshape(-1, 3)
        else:
            raise ValueError(f"Unsupported grid dimensionality: {len(tissue_grid)}")
        
        for field_name, cell_values in cellular_data.items():
            # Interpolate cellular data to tissue grid using inverse distance weighting
            tissue_field_flat = griddata(
                cell_positions, cell_values, grid_points_flat,
                method='linear', fill_value=0.0
            )
            
            # Handle NaN values with nearest neighbor interpolation
            nan_mask = np.isnan(tissue_field_flat)
            if np.any(nan_mask):
                tissue_field_flat[nan_mask] = griddata(
                    cell_positions, cell_values, grid_points_flat[nan_mask],
                    method='nearest'
                )
            
            tissue_field = tissue_field_flat.reshape(grid_shape)
            tissue_data[field_name] = tissue_field
        
        return tissue_data
    
    def transfer_tissue_to_cellular(self, tissue_data: Dict[str, np.ndarray],
                                  tissue_grid: Tuple[np.ndarray, ...],
                                  cell_positions: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Transfer tissue-scale information to cellular scale.
        
        Args:
            tissue_data: Tissue-scale fields
            tissue_grid: Tissue grid coordinates
            cell_positions: Cell positions
            
        Returns:
            Cellular-scale data
        """
        cellular_data = {}
        
        for field_name, tissue_field in tissue_data.items():
            if len(tissue_grid) == 2:
                # 2D interpolation
                interpolator = RegularGridInterpolator(
                    tissue_grid, tissue_field,
                    method='linear', bounds_error=False, fill_value=0.0
                )
            elif len(tissue_grid) == 3:
                # 3D interpolation
                interpolator = RegularGridInterpolator(
                    tissue_grid, tissue_field,
                    method='linear', bounds_error=False, fill_value=0.0
                )
            else:
                raise ValueError(f"Unsupported grid dimensionality: {len(tissue_grid)}")
            
            cellular_values = interpolator(cell_positions)
            cellular_data[field_name] = cellular_values
        
        return cellular_data


class HomogenizationTechniques:
    """
    Implementation of homogenization techniques for upscaling.
    """
    
    def __init__(self):
        """Initialize homogenization techniques."""
        pass
    
    def volume_averaging(self, fine_scale_field: np.ndarray,
                        coarse_grid_shape: Tuple[int, ...],
                        averaging_kernel: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Perform volume averaging homogenization.
        
        Args:
            fine_scale_field: Fine-scale field data
            coarse_grid_shape: Target coarse grid shape
            averaging_kernel: Optional averaging kernel
            
        Returns:
            Coarse-scale homogenized field
        """
        if averaging_kernel is None:
            # Use uniform averaging kernel
            kernel_size = tuple(fine_scale_field.shape[i] // coarse_grid_shape[i] 
                              for i in range(len(coarse_grid_shape)))
            averaging_kernel = np.ones(kernel_size) / np.prod(kernel_size)
        
        # Perform convolution-based averaging
        from scipy.ndimage import convolve
        
        # Calculate stride for downsampling
        stride = tuple(fine_scale_field.shape[i] // coarse_grid_shape[i] 
                      for i in range(len(coarse_grid_shape)))
        
        # Apply averaging kernel
        averaged_field = convolve(fine_scale_field, averaging_kernel, mode='constant')
        
        # Downsample to coarse grid
        if len(coarse_grid_shape) == 2:
            coarse_field = averaged_field[::stride[0], ::stride[1]]
        elif len(coarse_grid_shape) == 3:
            coarse_field = averaged_field[::stride[0], ::stride[1], ::stride[2]]
        else:
            raise ValueError(f"Unsupported dimensionality: {len(coarse_grid_shape)}")
        
        return coarse_field
    
    def effective_property_calculation(self, fine_scale_properties: np.ndarray,
                                     volume_fractions: np.ndarray,
                                     homogenization_method: str = "arithmetic") -> float:
        """
        Calculate effective properties using homogenization.
        
        Args:
            fine_scale_properties: Fine-scale property values
            volume_fractions: Volume fractions of each phase
            homogenization_method: "arithmetic", "harmonic", or "geometric"
            
        Returns:
            Effective property value
        """
        if homogenization_method == "arithmetic":
            # Arithmetic mean (Voigt bound)
            effective_property = np.sum(fine_scale_properties * volume_fractions)
            
        elif homogenization_method == "harmonic":
            # Harmonic mean (Reuss bound)
            effective_property = 1.0 / np.sum(volume_fractions / fine_scale_properties)
            
        elif homogenization_method == "geometric":
            # Geometric mean
            effective_property = np.prod(fine_scale_properties ** volume_fractions)
            
        else:
            raise ValueError(f"Unknown homogenization method: {homogenization_method}")
        
        return effective_property
    
    def multiscale_asymptotic_expansion(self, fine_scale_solution: np.ndarray,
                                      epsilon: float,
                                      expansion_order: int = 2) -> Dict[str, np.ndarray]:
        """
        Perform multiscale asymptotic expansion.
        
        Args:
            fine_scale_solution: Fine-scale solution
            epsilon: Scale separation parameter
            expansion_order: Order of asymptotic expansion
            
        Returns:
            Expansion coefficients
        """
        # This is a simplified implementation
        # In practice, this would involve solving cell problems
        
        expansion_coefficients = {}
        
        # Zeroth order (homogenized solution)
        expansion_coefficients['u0'] = self.volume_averaging(
            fine_scale_solution, 
            tuple(s//4 for s in fine_scale_solution.shape)
        )
        
        if expansion_order >= 1:
            # First order correction
            gradient = np.gradient(expansion_coefficients['u0'])
            expansion_coefficients['u1'] = epsilon * gradient[0]  # Simplified
        
        if expansion_order >= 2:
            # Second order correction
            laplacian = np.gradient(np.gradient(expansion_coefficients['u0'])[0])[0]
            expansion_coefficients['u2'] = epsilon**2 * laplacian  # Simplified
        
        return expansion_coefficients


class UpscalingOperator:
    """
    Operator for transferring information from fine to coarse scales.
    """
    
    def __init__(self, coupling_interface: CouplingInterface):
        """Initialize upscaling operator."""
        self.interface = coupling_interface
        self.homogenizer = HomogenizationTechniques()
    
    def transfer(self, source_data: Dict[str, np.ndarray],
                target_data: Dict[str, np.ndarray],
                dt: float) -> Dict[str, np.ndarray]:
        """
        Transfer data from fine to coarse scale.
        
        Args:
            source_data: Fine-scale data
            target_data: Coarse-scale data
            dt: Time step
            
        Returns:
            Updated coarse-scale data
        """
        transferred_data = {}
        
        for var_name in self.interface.transfer_variables:
            if var_name in source_data:
                fine_field = source_data[var_name]
                
                # Determine target grid shape
                if var_name in target_data:
                    target_shape = target_data[var_name].shape
                else:
                    # Default coarsening factor of 4
                    target_shape = tuple(s//4 for s in fine_field.shape)
                
                # Apply homogenization
                coarse_field = self.homogenizer.volume_averaging(
                    fine_field, target_shape
                )
                
                # Apply coupling strength
                if var_name in target_data:
                    # Weighted combination with existing data
                    alpha = self.interface.coupling_strength
                    coarse_field = alpha * coarse_field + (1.0 - alpha) * target_data[var_name]
                
                transferred_data[var_name] = coarse_field
        
        return transferred_data


class DownscalingOperator:
    """
    Operator for transferring information from coarse to fine scales.
    """
    
    def __init__(self, coupling_interface: CouplingInterface):
        """Initialize downscaling operator."""
        self.interface = coupling_interface
    
    def transfer(self, source_data: Dict[str, np.ndarray],
                target_data: Dict[str, np.ndarray],
                dt: float) -> Dict[str, np.ndarray]:
        """
        Transfer data from coarse to fine scale.
        
        Args:
            source_data: Coarse-scale data
            target_data: Fine-scale data
            dt: Time step
            
        Returns:
            Updated fine-scale data
        """
        transferred_data = {}
        
        for var_name in self.interface.transfer_variables:
            if var_name in source_data:
                coarse_field = source_data[var_name]
                
                # Determine target grid shape
                if var_name in target_data:
                    target_shape = target_data[var_name].shape
                else:
                    # Default refinement factor of 4
                    target_shape = tuple(s*4 for s in coarse_field.shape)
                
                # Apply interpolation for refinement
                fine_field = self._interpolate_to_fine_grid(
                    coarse_field, target_shape
                )
                
                # Apply coupling strength
                if var_name in target_data:
                    # Weighted combination with existing data
                    alpha = self.interface.coupling_strength
                    fine_field = alpha * fine_field + (1.0 - alpha) * target_data[var_name]
                
                transferred_data[var_name] = fine_field
        
        return transferred_data
    
    def _interpolate_to_fine_grid(self, coarse_field: np.ndarray,
                                target_shape: Tuple[int, ...]) -> np.ndarray:
        """Interpolate coarse field to fine grid."""
        from scipy.ndimage import zoom
        
        # Calculate zoom factors
        zoom_factors = tuple(target_shape[i] / coarse_field.shape[i] 
                           for i in range(len(target_shape)))
        
        # Apply interpolation
        fine_field = zoom(coarse_field, zoom_factors, order=1)
        
        return fine_field
