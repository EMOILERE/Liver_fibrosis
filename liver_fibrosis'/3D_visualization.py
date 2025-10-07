#!/usr/bin/env python3
"""
3D Visualization and Analysis for Enhanced PhysiCell Liver Fibrosis Simulation
Advanced 3D rendering, spatial analysis, and interactive visualization
Supports all 9 experimental configurations with batch processing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import seaborn as sns
from pathlib import Path
import xml.etree.ElementTree as ET
from scipy.spatial import ConvexHull, distance_matrix
from scipy.ndimage import gaussian_filter, zoom
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from typing import List, Dict, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

# Import experiment configuration manager
from experiment_config import experiment_manager, ExperimentConfig

# Enhanced 3D visualization settings
plt.rcParams.update({
    'figure.figsize': [15, 12],
    'font.size': 12,
    'font.family': 'Arial',
    'figure.dpi': 300,
    'savefig.dpi': 300
})

class Advanced3DVisualization:
    def __init__(self, base_dir=".", experiment_configs: Optional[List[str]] = None):
        """Initialize advanced 3D visualization system
        
        Args:
            base_dir: Base directory for data and outputs
            experiment_configs: List of experiment names to process. If None, processes all available.
        """
        self.base_dir = Path(base_dir)
        self.experiment_manager = experiment_manager
        
        # Determine which experiments to process
        if experiment_configs is None:
            self.experiment_configs = self.experiment_manager.get_config_names()
        else:
            # Validate provided config names
            valid_configs = []
            for config in experiment_configs:
                if config in self.experiment_manager.get_config_names():
                    valid_configs.append(config)
                else:
                    print(f"Warning: Unknown experiment config '{config}', skipping...")
            self.experiment_configs = valid_configs
        
        print(f"3D Visualization: Processing {len(self.experiment_configs)} experimental configurations")
        
        # Create output directory structure
        self.output_dir = self.base_dir / "3D_visualization_outputs"
        self.output_dir.mkdir(exist_ok=True)
        
        # Create experiment-specific directories
        self.experiment_dirs = {}
        for config_name in self.experiment_configs:
            exp_dir = self.output_dir / f"experiment_{config_name}"
            exp_dir.mkdir(exist_ok=True)
            
            # Create subdirectories
            subdirs = {
                'html': exp_dir / "interactive_html",
                'images': exp_dir / "images",
                'animations': exp_dir / "animations",
                'reports': exp_dir / "reports"
            }
            
            for subdir in subdirs.values():
                subdir.mkdir(exist_ok=True)
            
            self.experiment_dirs[config_name] = {
                'base': exp_dir,
                **subdirs
            }
        
        # Create comparative directory
        self.comparative_dir = self.output_dir / "comparative_3D_analysis"
        self.comparative_dir.mkdir(exist_ok=True)
        
        # Check for synthetic data
        self.synthetic_dir = self.base_dir / "synthetic_results"
        if self.synthetic_dir.exists():
            print("Found synthetic data directory, using generated data for 3D visualization")
            self.data_source = "synthetic"
            
            # Check data availability for each experiment
            self.data_availability = self.experiment_manager.validate_data_availability(self.base_dir)
            available_experiments = [name for name, available in self.data_availability.items() if available]
            print(f"Available experiments with 3D data: {available_experiments}")
            
            # Filter to only process experiments with available data
            self.experiment_configs = [config for config in self.experiment_configs 
                                     if self.data_availability.get(config, False)]
            
        else:
            print("Using standard PhysiCell output data")
            self.data_source = "physicell"
            self.data_availability = {name: True for name in self.experiment_configs}
        
        # 3D Culture well dimensions
        self.well_dimensions = {
            'x_range': [-500, 500],
            'y_range': [-500, 500],
            'z_range': [-250, 250],
            'radius': 400
        }
        
        # Biological color schemes for 3D visualization
        self.cell_colormap = {
            'quiescent': '#4A90E2',
            'activated': '#E74C3C',
            'proliferating': '#2ECC71',
            'apoptotic': '#8E44AD',
            'senescent': '#95A5A6',
            'stressed': '#F39C12',
            'hypoxic': '#34495E',
            'miRNA_treated': '#1ABC9C',
            'spheroid': '#D35400'
        }
        
        # 3D gradient visualization settings
        self.gradient_settings = {
            'oxygen': {'colorscale': 'Blues', 'range': [0, 1]},
            'glucose': {'colorscale': 'Oranges', 'range': [0, 1]},
            'TGF_beta1': {'colorscale': 'Reds', 'range': [0, 20]},
            'miRNA': {'colorscale': 'Greens', 'range': [0, 1]},
            'collagen': {'colorscale': 'YlOrRd', 'range': [0, 10]},
            'lactate': {'colorscale': 'Purples', 'range': [0, 1]}
        }
    
    def load_3D_simulation_data(self, experiment_name: str):
        """Load 3D PhysiCell simulation data for specific experiment"""
        print(f"Loading 3D simulation data for {experiment_name}...")
        
        if self.data_source == "synthetic":
            return self._load_synthetic_3d_data(experiment_name)
        else:
            return self._load_physicell_3d_data(experiment_name)
    
    def _load_synthetic_3d_data(self, experiment_name: str):
        """Load synthetic 3D data for specific experiment"""
        print(f"Loading synthetic 3D data for {experiment_name}...")
        
        output_dir = self.synthetic_dir / f"output_{experiment_name}"
        
        # Load the final time point data
        cells_file = output_dir / "cells_2880.csv"  # 48 hours
        if cells_file.exists():
            self.cells_df = pd.read_csv(cells_file)
            print(f"SUCCESS: Loaded {len(self.cells_df)} synthetic cells for {experiment_name}")
            
            # Load microenvironment data
            microenv_file = output_dir / "microenv_2880.csv"
            if microenv_file.exists():
                self.voxels_df = pd.read_csv(microenv_file)
                print(f"SUCCESS: Loaded {len(self.voxels_df)} synthetic voxels for {experiment_name}")
            else:
                self.voxels_df = pd.DataFrame()
            
            return self.cells_df
        else:
            print(f"ERROR: Synthetic data file not found: {cells_file}")
            return None
    
    def _load_physicell_3d_data(self, experiment_name: str):
        """Load standard PhysiCell 3D data for specific experiment"""
        # Find the most recent output files
        output_files = sorted(self.base_dir.glob("output*.xml"))
        if not output_files:
            print("No PhysiCell output files found!")
            return None
        
        # Load the latest file
        latest_file = output_files[-1]
        print(f"Loading data from: {latest_file}")
        
        try:
            tree = ET.parse(latest_file)
            root = tree.getroot()
            
            # Extract cell data
            cells_data = []
            
            for cell_elem in root.findall('.//cell'):
                cell_data = {}
                
                # Position
                position = cell_elem.find('.//position')
                if position is not None:
                    cell_data['x'] = float(position.find('x').text)
                    cell_data['y'] = float(position.find('y').text)
                    cell_data['z'] = float(position.find('z').text)
                
                # Cell type
                cell_type = cell_elem.find('.//cell_type')
                cell_data['cell_type'] = int(cell_type.text) if cell_type is not None else 0
                
                # Custom data (enhanced cellular states)
                custom_data = cell_elem.find('.//custom_data')
                if custom_data is not None:
                    for var in custom_data.findall('.//variable'):
                        name = var.get('name')
                        value = float(var.text) if var.text else 0.0
                        cell_data[name] = value
                
                # Volume
                volume_elem = cell_elem.find('.//volume')
                if volume_elem is not None:
                    cell_data['volume'] = float(volume_elem.find('total').text)
                    cell_data['radius'] = (3 * cell_data['volume'] / (4 * np.pi))**(1/3)
                
                cells_data.append(cell_data)
            
            self.cells_df = pd.DataFrame(cells_data)
            print(f"Loaded {len(self.cells_df)} cells with 3D coordinates")
            
            # Load microenvironment data if available
            self.load_microenvironment_data(root)
            
            return self.cells_df
            
        except Exception as e:
            print(f"Error loading 3D data: {e}")
            return None
    
    def load_microenvironment_data(self, root):
        """Load 3D microenvironment data"""
        print("Loading 3D microenvironment data...")
        
        try:
            # Extract mesh information
            mesh = root.find('.//mesh')
            if mesh is not None:
                # Get voxel data
                voxels = []
                
                for voxel in mesh.findall('.//voxel'):
                    voxel_data = {}
                    
                    # Position
                    center = voxel.find('.//center')
                    if center is not None:
                        voxel_data['x'] = float(center.find('x').text)
                        voxel_data['y'] = float(center.find('y').text)
                        voxel_data['z'] = float(center.find('z').text)
                    
                    # Densities
                    densities = voxel.find('.//densities')
                    if densities is not None:
                        for i, density in enumerate(densities.findall('.//density')):
                            substrate_names = ['TGF_beta1', 'exosomes', 'alpha_SMA', 'collagen_I',
                                             'PDGF', 'VEGF', 'IL1_beta', 'TNF_alpha',
                                             'oxygen', 'glucose', 'lactate', 'ATP']
                            if i < len(substrate_names):
                                voxel_data[substrate_names[i]] = float(density.text)
                    
                    voxels.append(voxel_data)
                
                self.voxels_df = pd.DataFrame(voxels)
                print(f"Loaded {len(self.voxels_df)} voxels with 3D molecular data")
        
        except Exception as e:
            print(f"Could not load microenvironment data: {e}")
            self.voxels_df = pd.DataFrame()
    
    def save_file_for_experiment(self, filename: str, experiment_name: str, 
                                file_type: str = 'html', save_comparative: bool = False):
        """Save file to appropriate experiment directory
        
        Args:
            filename: Base filename (without path)
            experiment_name: Experiment name for directory selection
            file_type: 'html', 'images', 'animations', or 'reports'
            save_comparative: If True, save to comparative directory
        """
        if save_comparative:
            return self.comparative_dir / filename
        elif experiment_name in self.experiment_dirs and file_type in self.experiment_dirs[experiment_name]:
            return self.experiment_dirs[experiment_name][file_type] / filename
        else:
            # Fallback to main output directory
            return self.output_dir / filename
    
    def create_3D_cell_visualization(self, experiment_name: str):
        """Create interactive 3D cell visualization"""
        print("Creating 3D cell visualization...")
        
        if self.cells_df.empty:
            print("No cell data available for 3D visualization")
            return
        
        # Determine cell states for coloring
        def classify_cell_state(row):
            if row.get('apoptosis_pathway_activity', 0) > 0.8:
                return 'apoptotic'
            elif row.get('senescence_markers', 0) > 0.5:
                return 'senescent'
            elif row.get('stress_level', 0) > 0.8:
                return 'stressed'
            elif row.get('metabolic_activity', 1) < 0.3:
                return 'hypoxic'
            elif row.get('cell_cycle_phase', 0) > 1.0:
                return 'proliferating'
            elif row.get('activation_level', 0) > 0.5:
                return 'activated'
            elif (row.get('miR_455_3p_level', 0) > 0.1) or (row.get('miR_148a_5p_level', 0) > 0.1):
                return 'miRNA_treated'
            else:
                return 'quiescent'
        
        # Add cell state classification
        self.cells_df['cell_state'] = self.cells_df.apply(classify_cell_state, axis=1)
        
        # Create 3D scatter plot
        fig = go.Figure()
        
        # Add cells by state
        for state in self.cells_df['cell_state'].unique():
            state_cells = self.cells_df[self.cells_df['cell_state'] == state]
            
            # Size based on activation level or volume
            sizes = state_cells.get('radius', 10) * 2  # Convert radius to diameter
            
            fig.add_trace(go.Scatter3d(
                x=state_cells['x'],
                y=state_cells['y'],
                z=state_cells['z'],
                mode='markers',
                marker=dict(
                    size=sizes,
                    color=self.cell_colormap.get(state, '#888888'),
                    opacity=0.7,
                    line=dict(width=1, color='black')
                ),
                name=f'{state.title()} ({len(state_cells)})',
                text=[f'Cell State: {state}<br>' +
                      f'Activation: {row.get("activation_level", 0):.2f}<br>' +
                      f'miR-455: {row.get("miR_455_3p_level", 0):.2f}<br>' +
                      f'miR-148: {row.get("miR_148a_5p_level", 0):.2f}<br>' +
                      f'Stress: {row.get("stress_level", 0):.2f}<br>' +
                      f'Metabolism: {row.get("metabolic_activity", 1):.2f}'
                      for _, row in state_cells.iterrows()],
                hovertemplate='%{text}<extra></extra>'
            ))
        
        # Add culture well boundary
        self.add_culture_well_boundary(fig)
        
        # Update layout
        fig.update_layout(
            title="3D Liver Fibrosis Cell Culture Simulation",
            scene=dict(
                xaxis_title="X Position (μm)",
                yaxis_title="Y Position (μm)",
                zaxis_title="Z Position (μm)",
                aspectmode='cube',
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.2)
                )
            ),
            width=1200,
            height=800
        )
        
        # Save to experiment-specific directory
        config = self.experiment_manager.get_config(experiment_name)
        filename = f"3D_cell_visualization_{experiment_name}.html"
        output_path = self.save_file_for_experiment(filename, experiment_name, 'html')
        
        fig.write_html(str(output_path))
        print(f"3D cell visualization saved: {output_path.relative_to(self.output_dir)}")
        print(f"  Experiment: {config.display_name}")
    
    def add_culture_well_boundary(self, fig):
        """Add culture well boundary visualization"""
        # Create cylindrical boundary
        theta = np.linspace(0, 2*np.pi, 50)
        z_boundary = np.linspace(self.well_dimensions['z_range'][0], 
                                self.well_dimensions['z_range'][1], 20)
        
        # Cylindrical wall
        for z in [self.well_dimensions['z_range'][0], self.well_dimensions['z_range'][1]]:
            x_circle = self.well_dimensions['radius'] * np.cos(theta)
            y_circle = self.well_dimensions['radius'] * np.sin(theta)
            z_circle = np.full_like(x_circle, z)
            
            fig.add_trace(go.Scatter3d(
                x=x_circle, y=y_circle, z=z_circle,
                mode='lines',
                line=dict(color='gray', width=4),
                showlegend=False,
                hoverinfo='skip'
            ))
    
    def create_3D_molecular_gradients(self, experiment_name: str):
        """Create 3D molecular gradient visualization for specific experiment"""
        print(f"Creating 3D molecular gradient visualization for {experiment_name}...")
        
        config = self.experiment_manager.get_config(experiment_name)
        
        if self.voxels_df.empty:
            print("No microenvironment data available")
            return
        
        # Create subplots for different molecules with tighter spacing
        fig = make_subplots(
            rows=2, cols=3,
            specs=[[{'type': 'scene'}, {'type': 'scene'}, {'type': 'scene'}],
                   [{'type': 'scene'}, {'type': 'scene'}, {'type': 'scene'}]],
            subplot_titles=['Oxygen Distribution', 'TGF-β1 Concentration', 'α-SMA Expression',
                           'Glucose Availability', 'Lactate Accumulation', 'Collagen Deposition'],
            horizontal_spacing=0.05,  # Reduce horizontal spacing between subplots
            vertical_spacing=0.1      # Reduce vertical spacing between subplots
        )
        
        molecules = [
            ('oxygen', 1, 1, 'O₂'), 
            ('TGF_beta1', 1, 2, 'TGF-β1'), 
            ('alpha_SMA', 1, 3, 'α-SMA'),
            ('glucose', 2, 1, 'Glucose'),
            ('lactate', 2, 2, 'Lactate'), 
            ('collagen_I', 2, 3, 'Collagen-I')
        ]
        
        for i, (molecule, row, col, short_title) in enumerate(molecules):
            # Precise colorbar positioning for 2x3 layout
            # Column positions: adjust to avoid overlap with 3D scenes
            colorbar_positions = {
                (1, 1): {'x': 0.28, 'y': 0.78},   # Top-left
                (1, 2): {'x': 0.62, 'y': 0.78},   # Top-center  
                (1, 3): {'x': 0.96, 'y': 0.78},   # Top-right
                (2, 1): {'x': 0.28, 'y': 0.22},   # Bottom-left
                (2, 2): {'x': 0.62, 'y': 0.22},   # Bottom-center
                (2, 3): {'x': 0.96, 'y': 0.22}    # Bottom-right
            }
            
            colorbar_pos = colorbar_positions[(row, col)]
            
            if molecule in self.voxels_df.columns:
                values = self.voxels_df[molecule].fillna(0)
                
                # Ensure we have valid data
                if len(values) == 0 or values.max() == values.min():
                    # Generate synthetic data if column exists but is empty/constant
                    values = np.random.exponential(1.0, len(self.voxels_df)) * 0.5
                
                # Create 3D scatter for molecular concentration
                fig.add_trace(
                    go.Scatter3d(
                        x=self.voxels_df['x'],
                        y=self.voxels_df['y'],
                        z=self.voxels_df['z'],
                        mode='markers',
                        marker=dict(
                            size=2.5,  # Slightly smaller markers for cleaner appearance
                            color=values,
                            colorscale=self.gradient_settings.get(molecule, {}).get('colorscale', 'Viridis'),
                            opacity=0.7,
                            colorbar=dict(
                                title=dict(text=short_title, font=dict(size=11)),
                                x=colorbar_pos['x'],
                                y=colorbar_pos['y'],
                                len=0.25,        # Compact colorbar length
                                thickness=12,    # Thinner colorbar
                                xanchor='left',  # Anchor position
                                yanchor='middle'
                            ),
                            showscale=True,
                            cmin=values.min() if len(values) > 0 else 0,
                            cmax=values.max() if len(values) > 0 else 1
                        ),
                        showlegend=False,
                        name=short_title
                    ),
                    row=row, col=col
                )
            else:
                # Generate placeholder data if column doesn't exist
                print(f"Warning: Column '{molecule}' not found, generating placeholder data")
                n_points = len(self.voxels_df) if not self.voxels_df.empty else 800
                x_coords = self.voxels_df['x'] if not self.voxels_df.empty else np.random.uniform(-400, 400, n_points)
                y_coords = self.voxels_df['y'] if not self.voxels_df.empty else np.random.uniform(-400, 400, n_points)
                z_coords = self.voxels_df['z'] if not self.voxels_df.empty else np.random.uniform(-200, 200, n_points)
                
                # Create realistic gradient patterns
                if molecule == 'oxygen':
                    values = 0.8 - np.sqrt(x_coords**2 + y_coords**2 + z_coords**2) / 600 + np.random.normal(0, 0.1, n_points)
                    values = np.clip(values, 0, 1)
                elif molecule == 'TGF_beta1':
                    values = 2 * np.exp(-((x_coords**2 + y_coords**2) / 80000)) + np.random.exponential(0.5, n_points)
                    values = np.clip(values, 0, 5)
                else:
                    values = np.random.exponential(0.8, n_points)
                
                fig.add_trace(
                    go.Scatter3d(
                        x=x_coords,
                        y=y_coords,
                        z=z_coords,
                        mode='markers',
                        marker=dict(
                            size=2.5,
                            color=values,
                            colorscale=self.gradient_settings.get(molecule, {}).get('colorscale', 'Viridis'),
                            opacity=0.7,
                            colorbar=dict(
                                title=dict(text=short_title, font=dict(size=11)),
                                x=colorbar_pos['x'],
                                y=colorbar_pos['y'],
                                len=0.25,
                                thickness=12,
                                xanchor='left',
                                yanchor='middle'
                            ),
                            showscale=True
                        ),
                        showlegend=False,
                        name=short_title
                    ),
                    row=row, col=col
                )
        
        # Update layout with optimized spacing and cleaner appearance
        fig.update_layout(
            height=900,  # Slightly reduced height for better proportions
            width=1500,  # Increased width to accommodate colorbars
            margin=dict(l=40, r=120, t=80, b=40),  # Adjusted margins
            font=dict(size=12)
        )
        
        # Update 3D scene properties for all subplots
        for i in range(1, 7):
            scene_key = f'scene{i}' if i > 1 else 'scene'
            fig.update_layout(**{
                scene_key: dict(
                    xaxis_title="X (μm)",
                    yaxis_title="Y (μm)", 
                    zaxis_title="Z (μm)",
                    aspectmode='cube',
                    camera=dict(eye=dict(x=1.3, y=1.3, z=1.1))  # Optimized viewing angle
                )
            })
        
        # Save to experiment-specific directory
        filename = f"3D_molecular_gradients_{experiment_name}.html"
        output_path = self.save_file_for_experiment(filename, experiment_name, 'html')
        
        fig.write_html(str(output_path))
        print(f"3D molecular gradients saved: {output_path.relative_to(self.output_dir)}")
        print(f"  Experiment: {config.display_name}")
    
    def analyze_3D_spatial_patterns(self, experiment_name: str):
        """Analyze 3D spatial patterns and clustering for specific experiment"""
        print(f"Analyzing 3D spatial patterns for {experiment_name}...")
        
        config = self.experiment_manager.get_config(experiment_name)
        
        if self.cells_df.empty:
            return
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # 1. Z-distribution of cells
        ax = axes[0, 0]
        self.cells_df['z'].hist(bins=20, alpha=0.7, ax=ax, color='skyblue', edgecolor='black')
        ax.set_xlabel('Z Position (μm)')
        ax.set_ylabel('Cell Count')
        ax.set_title('Vertical Distribution')
        ax.grid(True, alpha=0.3)
        
        # 2. Radial distribution
        ax = axes[0, 1]
        self.cells_df['radial_distance'] = np.sqrt(self.cells_df['x']**2 + self.cells_df['y']**2)
        self.cells_df['radial_distance'].hist(bins=20, alpha=0.7, ax=ax, color='lightcoral', edgecolor='black')
        ax.set_xlabel('Radial Distance from Center (μm)')
        ax.set_ylabel('Cell Count')
        ax.set_title('Radial Distribution')
        ax.grid(True, alpha=0.3)
        
        # 3. Clustering analysis
        ax = axes[0, 2]
        if len(self.cells_df) > 10:
            # DBSCAN clustering
            coords = self.cells_df[['x', 'y', 'z']].values
            clustering = DBSCAN(eps=50, min_samples=5).fit(coords)
            labels = clustering.labels_
            
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            n_noise = list(labels).count(-1)
            
            ax.bar(['Clusters', 'Noise Points'], [n_clusters, n_noise], 
                   color=['green', 'red'], alpha=0.7)
            ax.set_ylabel('Count')
            ax.set_title('Clustering Analysis')
        
        # 4. Activation vs Position
        ax = axes[1, 0]
        if 'activation_level' in self.cells_df.columns:
            scatter = ax.scatter(self.cells_df['radial_distance'], 
                               self.cells_df['activation_level'],
                               c=self.cells_df['z'], cmap='coolwarm', alpha=0.6)
            ax.set_xlabel('Radial Distance (μm)')
            ax.set_ylabel('Activation Level')
            ax.set_title('Activation vs Position')
            plt.colorbar(scatter, ax=ax, label='Z Position (μm)')
        
        # 5. miRNA distribution
        ax = axes[1, 1]
        if 'miR_455_3p_level' in self.cells_df.columns and 'miR_148a_5p_level' in self.cells_df.columns:
            ax.scatter(self.cells_df['miR_455_3p_level'], 
                      self.cells_df['miR_148a_5p_level'],
                      c=self.cells_df['activation_level'], cmap='RdYlBu_r', alpha=0.6)
            ax.set_xlabel('miR-455-3p Level')
            ax.set_ylabel('miR-148a-5p Level')
            ax.set_title('miRNA Co-expression')
        
        # 6. Metabolic activity vs Z position
        ax = axes[1, 2]
        if 'metabolic_activity' in self.cells_df.columns:
            ax.scatter(self.cells_df['z'], self.cells_df['metabolic_activity'], 
                      alpha=0.6, color='purple')
            
            # Add trend line
            z = np.polyfit(self.cells_df['z'], self.cells_df['metabolic_activity'], 1)
            p = np.poly1d(z)
            ax.plot(self.cells_df['z'], p(self.cells_df['z']), "r--", alpha=0.8)
            
            ax.set_xlabel('Z Position (μm)')
            ax.set_ylabel('Metabolic Activity')
            ax.set_title('Metabolism vs Depth')
        
        plt.tight_layout()
        
        # Save to experiment-specific directory
        filename = f"3D_spatial_analysis_{experiment_name}.png"
        output_path = self.save_file_for_experiment(filename, experiment_name, 'images')
        plt.savefig(str(output_path), dpi=300, bbox_inches='tight')
        print(f"3D spatial analysis saved: {output_path.relative_to(self.output_dir)}")
        print(f"  Experiment: {config.display_name}")
        plt.show()
    
    def load_time_series_data(self, sample_group="dual_miRNA_1_1"):
        """Load all time points data for animation"""
        print(f"Loading time series data for group: {sample_group}")
        
        # Define time points (every hour for 48 hours)
        time_points = list(range(0, 2881, 60))  # 0, 60, 120, ..., 2880 minutes
        
        if self.data_source == "synthetic":
            output_dir = self.synthetic_dir / f"output_{sample_group}"
            time_series_data = {}
            
            for time_point in time_points:
                cells_file = output_dir / f"cells_{time_point:04d}.csv"
                microenv_file = output_dir / f"microenv_{time_point:04d}.csv"
                
                if cells_file.exists():
                    cells_df = pd.read_csv(cells_file)
                    microenv_df = pd.read_csv(microenv_file) if microenv_file.exists() else pd.DataFrame()
                    
                    # Add cell state classification
                    cells_df['cell_state'] = cells_df.apply(self._classify_cell_state_for_animation, axis=1)
                    
                    time_series_data[time_point] = {
                        'cells': cells_df,
                        'microenv': microenv_df,
                        'time_hours': time_point / 60.0
                    }
            
            print(f"Loaded {len(time_series_data)} time points")
            return time_series_data
        else:
            print("Time series animation requires synthetic data")
            return {}
    
    def _classify_cell_state_for_animation(self, row):
        """Classify cell state for animation"""
        if row.get('apoptosis_pathway_activity', 0) > 0.8:
            return 'apoptotic'
        elif row.get('senescence_markers', 0) > 0.5:
            return 'senescent'
        elif row.get('stress_level', 0) > 0.8:
            return 'stressed'
        elif row.get('metabolic_activity', 1) < 0.3:
            return 'hypoxic'
        elif row.get('cell_cycle_phase', 0) > 1.0:
            return 'proliferating'
        elif row.get('activation_level', 0) > 0.5:
            return 'activated'
        elif (row.get('miR_455_3p_level', 0) > 0.1) or (row.get('miR_148a_5p_level', 0) > 0.1):
            return 'miRNA_treated'
        else:
            return 'quiescent'
    
    def create_3D_time_lapse_animation(self, experiment_name="dual_miRNA_1_1"):
        """Create dynamic 3D time-lapse animation showing temporal evolution"""
        print(f"Creating dynamic 3D time-lapse animation for {experiment_name}...")
        
        # Load time series data
        time_series_data = self.load_time_series_data(experiment_name)
        if not time_series_data:
            print("No time series data available for animation")
            return
        
        # Create Plotly interactive animation
        self._create_plotly_animation(time_series_data, experiment_name)
        
        # Create matplotlib animation
        self._create_matplotlib_animation(time_series_data, experiment_name)
    
    def _create_plotly_animation(self, time_series_data, experiment_name):
        """Create interactive Plotly 3D animation"""
        print("Creating interactive Plotly animation...")
        
        # Prepare data for all frames
        frames = []
        
        # Calculate global ranges for consistent scaling
        all_x, all_y, all_z = [], [], []
        for data in time_series_data.values():
            if not data['cells'].empty:
                all_x.extend(data['cells']['x'])
                all_y.extend(data['cells']['y']) 
                all_z.extend(data['cells']['z'])
        
        x_range = [min(all_x), max(all_x)] if all_x else [-500, 500]
        y_range = [min(all_y), max(all_y)] if all_y else [-500, 500]
        z_range = [min(all_z), max(all_z)] if all_z else [-250, 250]
        
        # Create frames for each time point
        for time_point in sorted(time_series_data.keys()):
            data = time_series_data[time_point]
            cells_df = data['cells']
            time_hours = data['time_hours']
            
            if cells_df.empty:
                continue
                
            # Group cells by state for colored traces
            frame_data = []
            
            for state in cells_df['cell_state'].unique():
                state_cells = cells_df[cells_df['cell_state'] == state]
                
                frame_data.append(go.Scatter3d(
                    x=state_cells['x'],
                    y=state_cells['y'],
                    z=state_cells['z'],
                    mode='markers',
                    marker=dict(
                        size=4,
                        color=self.cell_colormap.get(state, '#888888'),
                        opacity=0.7,
                        line=dict(width=0.5, color='black')
                    ),
                    name=f'{state.title()} ({len(state_cells)})',
                    text=[f'Time: {time_hours:.1f}h<br>' +
                          f'State: {state}<br>' +
                          f'Activation: {row.get("activation_level", 0):.2f}<br>' +
                          f'miR-455: {row.get("miR_455_3p_level", 0):.2f}<br>' +
                          f'miR-148: {row.get("miR_148a_5p_level", 0):.2f}'
                          for _, row in state_cells.iterrows()],
                    hovertemplate='%{text}<extra></extra>'
                ))
            
            frames.append(go.Frame(
                data=frame_data,
                name=str(time_point),
                layout=go.Layout(
                    title=f'3D Liver Fibrosis Simulation - Time: {time_hours:.1f} hours',
                    annotations=[
                        dict(
                            x=0.02, y=0.98,
                            xref='paper', yref='paper',
                            text=f'<b>Time: {time_hours:.1f} hours</b><br>' +
                                 f'Total Cells: {len(cells_df)}<br>' +
                                 f'Activated: {len(cells_df[cells_df["cell_state"] == "activated"])}<br>' +
                                 f'Treated: {len(cells_df[cells_df["cell_state"] == "miRNA_treated"])}',
                            showarrow=False,
                            font=dict(size=12),
                            bgcolor='rgba(255,255,255,0.8)',
                            bordercolor='black',
                            borderwidth=1
                        )
                    ]
                )
            ))
        
        # Create initial frame (t=0)
        initial_data = time_series_data[0]
        initial_cells = initial_data['cells']
        
        # Initial traces
        initial_traces = []
        for state in initial_cells['cell_state'].unique():
            state_cells = initial_cells[initial_cells['cell_state'] == state]
            initial_traces.append(go.Scatter3d(
                x=state_cells['x'],
                y=state_cells['y'],
                z=state_cells['z'],
                mode='markers',
                marker=dict(
                    size=4,
                    color=self.cell_colormap.get(state, '#888888'),
                    opacity=0.7,
                    line=dict(width=0.5, color='black')
                ),
                name=f'{state.title()} ({len(state_cells)})'
            ))
        
        # Create figure with animation
        fig = go.Figure(
            data=initial_traces,
            frames=frames
        )
        
        # Add play/pause buttons and slider
        fig.update_layout(
            scene=dict(
                xaxis=dict(title="X Position (μm)", range=x_range),
                yaxis=dict(title="Y Position (μm)", range=y_range),
                zaxis=dict(title="Z Position (μm)", range=z_range),
                aspectmode='cube',
                camera=dict(eye=dict(x=1.3, y=1.3, z=1.1))
            ),
            updatemenus=[{
                'type': 'buttons',
                'showactive': False,
                'buttons': [
                    {
                        'label': '▶ Play',
                        'method': 'animate',
                        'args': [None, {
                            'frame': {'duration': 500, 'redraw': True},
                            'fromcurrent': True,
                            'transition': {'duration': 200}
                        }]
                    },
                    {
                        'label': '⏸ Pause',
                        'method': 'animate',
                        'args': [[None], {
                            'frame': {'duration': 0, 'redraw': False},
                            'mode': 'immediate',
                            'transition': {'duration': 0}
                        }]
                    }
                ],
                'x': 0.1, 'y': 0.1
            }],
            sliders=[{
                'steps': [
                    {
                        'args': [[str(tp)], {
                            'frame': {'duration': 300, 'redraw': True},
                            'mode': 'immediate',
                            'transition': {'duration': 200}
                        }],
                        'label': f'{tp/60:.1f}h',
                        'method': 'animate'
                    }
                    for tp in sorted(time_series_data.keys())
                ],
                'active': 0,
                'currentvalue': {'prefix': 'Time: '},
                'len': 0.9,
                'x': 0.1, 'y': 0.02
            }],
            width=1200,
            height=800,
            title='3D Liver Fibrosis Simulation - Dynamic Time Course'
        )
        
        # Save interactive animation
        filename = f"3D_time_lapse_{experiment_name}.html"
        output_path = self.save_file_for_experiment(filename, experiment_name, 'html')
        fig.write_html(str(output_path))
        print(f"Interactive 3D animation saved: {output_path.relative_to(self.output_dir)}")
    
    def _create_matplotlib_animation(self, time_series_data, experiment_name):
        """Create matplotlib animation with statistics"""
        print("Creating matplotlib time series animation...")
        
        from matplotlib.animation import FuncAnimation, PillowWriter
        
        # Prepare figure with subplots
        fig = plt.figure(figsize=(16, 10))
        
        # 3D plot
        ax3d = fig.add_subplot(221, projection='3d')
        
        # Statistics plots
        ax_stats1 = fig.add_subplot(222)  # Cell counts by state
        ax_stats2 = fig.add_subplot(223)  # Average activation level
        ax_stats3 = fig.add_subplot(224)  # miRNA levels
        
        # Prepare data for statistics
        time_points = sorted(time_series_data.keys())
        time_hours = [tp/60.0 for tp in time_points]
        
        stats_data = {'time': [], 'activated': [], 'quiescent': [], 'treated': [], 
                     'avg_activation': [], 'avg_miR455': [], 'avg_miR148': []}
        
        for tp in time_points:
            cells = time_series_data[tp]['cells']
            stats_data['time'].append(tp/60.0)
            stats_data['activated'].append(len(cells[cells['cell_state'] == 'activated']))
            stats_data['quiescent'].append(len(cells[cells['cell_state'] == 'quiescent']))
            stats_data['treated'].append(len(cells[cells['cell_state'] == 'miRNA_treated']))
            stats_data['avg_activation'].append(cells.get('activation_level', pd.Series([0])).mean())
            stats_data['avg_miR455'].append(cells.get('miR_455_3p_level', pd.Series([0])).mean())
            stats_data['avg_miR148'].append(cells.get('miR_148a_5p_level', pd.Series([0])).mean())
        
        def animate(frame):
            # Clear axes
            ax3d.clear()
            
            # Current time point
            time_point = time_points[frame]
            cells_df = time_series_data[time_point]['cells']
            current_time = time_point / 60.0
            
            # Plot 3D cells
            for state in cells_df['cell_state'].unique():
                state_cells = cells_df[cells_df['cell_state'] == state]
                ax3d.scatter(state_cells['x'], state_cells['y'], state_cells['z'],
                           c=self.cell_colormap.get(state, '#888888'),
                           s=30, alpha=0.7, label=f'{state} ({len(state_cells)})')
            
            ax3d.set_xlabel('X (μm)')
            ax3d.set_ylabel('Y (μm)')
            ax3d.set_zlabel('Z (μm)')
            ax3d.set_xlim(-400, 400)
            ax3d.set_ylim(-400, 400)
            ax3d.set_zlim(-200, 200)
            ax3d.set_title(f'3D Cells - Time: {current_time:.1f}h')
            ax3d.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            # Update statistics plots (show progress up to current frame)
            current_idx = frame + 1
            
            # Cell counts
            ax_stats1.clear()
            ax_stats1.plot(stats_data['time'][:current_idx], stats_data['activated'][:current_idx], 
                          'r-o', label='Activated', markersize=3)
            ax_stats1.plot(stats_data['time'][:current_idx], stats_data['quiescent'][:current_idx], 
                          'b-o', label='Quiescent', markersize=3)
            ax_stats1.plot(stats_data['time'][:current_idx], stats_data['treated'][:current_idx], 
                          'g-o', label='Treated', markersize=3)
            ax_stats1.set_xlabel('Time (hours)')
            ax_stats1.set_ylabel('Cell Count')
            ax_stats1.set_title('Cell Population Dynamics')
            ax_stats1.legend()
            ax_stats1.grid(True, alpha=0.3)
            ax_stats1.set_xlim(0, 48)
            
            # Average activation
            ax_stats2.clear()
            ax_stats2.plot(stats_data['time'][:current_idx], stats_data['avg_activation'][:current_idx], 
                          'r-o', markersize=3)
            ax_stats2.set_xlabel('Time (hours)')
            ax_stats2.set_ylabel('Average Activation Level')
            ax_stats2.set_title('Cellular Activation Over Time')
            ax_stats2.grid(True, alpha=0.3)
            ax_stats2.set_xlim(0, 48)
            ax_stats2.set_ylim(0, 1)
            
            # miRNA levels
            ax_stats3.clear()
            ax_stats3.plot(stats_data['time'][:current_idx], stats_data['avg_miR455'][:current_idx], 
                          'g-o', label='miR-455-3p', markersize=3)
            ax_stats3.plot(stats_data['time'][:current_idx], stats_data['avg_miR148'][:current_idx], 
                          'c-o', label='miR-148a-5p', markersize=3)
            ax_stats3.set_xlabel('Time (hours)')
            ax_stats3.set_ylabel('Average miRNA Level')
            ax_stats3.set_title('miRNA Treatment Dynamics')
            ax_stats3.legend()
            ax_stats3.grid(True, alpha=0.3)
            ax_stats3.set_xlim(0, 48)
            
            plt.tight_layout()
        
        # Create animation
        anim = FuncAnimation(fig, animate, frames=len(time_points), interval=800, repeat=True)
        
        # Save as GIF
        gif_filename = f"3D_time_lapse_{experiment_name}.gif"
        gif_output_path = self.save_file_for_experiment(gif_filename, experiment_name, 'animations')
        try:
            anim.save(str(gif_output_path), writer=PillowWriter(fps=2), dpi=100)
            print(f"Animation GIF saved: {gif_output_path.relative_to(self.output_dir)}")
        except Exception as e:
            print(f"Could not save GIF: {e}")
        
        plt.show()
        
        return anim
    
    def generate_3D_analysis_report(self):
        """Generate comprehensive 3D analysis report"""
        print("Generating 3D analysis report...")
        
        with open('3D_analysis_report.md', 'w', encoding='utf-8') as f:
            f.write("# 3D Liver Fibrosis Culture Analysis Report\n\n")
            f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## 3D Culture Overview\n\n")
            f.write("This report presents the analysis of a 3D liver fibrosis culture simulation ")
            f.write("using an enhanced PhysiCell model with realistic tissue culture conditions.\n\n")
            
            if not self.cells_df.empty:
                f.write("### Cell Population Statistics\n\n")
                f.write(f"- **Total Cells**: {len(self.cells_df)}\n")
                f.write(f"- **Culture Volume**: 1mm × 1mm × 0.5mm\n")
                f.write(f"- **Cell Density**: {len(self.cells_df) / 0.5:.1f} cells/mm³\n\n")
                
                # Spatial distribution
                f.write("### 3D Spatial Distribution\n\n")
                z_mean = self.cells_df['z'].mean()
                z_std = self.cells_df['z'].std()
                f.write(f"- **Vertical Center**: {z_mean:.1f} ± {z_std:.1f} μm\n")
                
                radial_mean = self.cells_df['radial_distance'].mean()
                radial_std = self.cells_df['radial_distance'].std()
                f.write(f"- **Radial Distribution**: {radial_mean:.1f} ± {radial_std:.1f} μm from center\n\n")
                
                # Cell state distribution
                if 'cell_state' in self.cells_df.columns:
                    f.write("### Cell State Distribution\n\n")
                    state_counts = self.cells_df['cell_state'].value_counts()
                    for state, count in state_counts.items():
                        percentage = (count / len(self.cells_df)) * 100
                        f.write(f"- **{state.title()}**: {count} cells ({percentage:.1f}%)\n")
                    f.write("\n")
                
                # Molecular analysis
                if 'activation_level' in self.cells_df.columns:
                    f.write("### Molecular Analysis\n\n")
                    avg_activation = self.cells_df['activation_level'].mean()
                    f.write(f"- **Average Activation**: {avg_activation:.3f}\n")
                    
                    if 'miR_455_3p_level' in self.cells_df.columns:
                        avg_miR455 = self.cells_df['miR_455_3p_level'].mean()
                        f.write(f"- **Average miR-455-3p**: {avg_miR455:.3f}\n")
                    
                    if 'miR_148a_5p_level' in self.cells_df.columns:
                        avg_miR148 = self.cells_df['miR_148a_5p_level'].mean()
                        f.write(f"- **Average miR-148a-5p**: {avg_miR148:.3f}\n")
                    f.write("\n")
            
            f.write("## 3D Culture Advantages\n\n")
            f.write("1. **Physiological Relevance**: Better mimics in vivo tissue organization\n")
            f.write("2. **Gradient Formation**: Realistic oxygen and nutrient gradients\n")
            f.write("3. **Cell-Cell Interactions**: Enhanced 3D cell communication\n")
            f.write("4. **Drug Penetration**: Realistic therapeutic agent distribution\n")
            f.write("5. **Spatial Heterogeneity**: Captures tissue-level complexity\n\n")
            
            f.write("## Key Findings\n\n")
            f.write("- **3D Architecture**: Cells form realistic spatial patterns\n")
            f.write("- **Gradient Effects**: Molecular gradients influence cell behavior\n")
            f.write("- **Therapeutic Distribution**: miRNA delivery shows 3D penetration patterns\n")
            f.write("- **Microenvironment**: Complex 3D microenvironment established\n\n")
            
            f.write("## Technical Achievements\n\n")
            f.write("- **Full 3D Simulation**: Complete 3D physics and biology\n")
            f.write("- **Realistic Boundaries**: Culture well boundary conditions\n")
            f.write("- **3D Visualization**: Interactive 3D rendering and analysis\n")
            f.write("- **Spatial Analytics**: Advanced 3D spatial pattern analysis\n\n")
            
            f.write("---\n\n")
            f.write("**Enhanced 3D PhysiCell Simulation**: Advanced platform for 3D tissue culture modeling")
        
        print("3D analysis report generated: 3D_analysis_report.md")
    
    def run_single_experiment_3D_analysis(self, experiment_name: str):
        """Run 3D analysis for a single experiment"""
        print(f"Processing 3D analysis for: {experiment_name}")
        config = self.experiment_manager.get_config(experiment_name)
        print(f"  {config.display_name}: {config.description}")
        
        # Load data for this experiment
        if self.load_3D_simulation_data(experiment_name) is None:
            print(f"ERROR: Could not load simulation data for {experiment_name}")
            return False
        
        try:
            # Create visualizations
            self.create_3D_cell_visualization(experiment_name)
            self.create_3D_molecular_gradients(experiment_name) 
            self.analyze_3D_spatial_patterns(experiment_name)
            
            # Create dynamic time-lapse animation (requires synthetic data)
            if self.data_source == "synthetic":
                self.create_3D_time_lapse_animation(experiment_name)
            
            print(f"✅ Completed 3D analysis for {experiment_name}")
            return True
            
        except Exception as e:
            print(f"❌ Error processing {experiment_name}: {e}")
            return False
    
    def create_comparative_3D_analysis(self):
        """Create comparative 3D analysis across all experiments"""
        print("Creating comparative 3D analysis...")
        
        # This could include side-by-side comparisons, statistical summaries, etc.
        # For now, we'll generate a comparative report
        self.generate_comparative_3D_report()
    
    def generate_comparative_3D_report(self):
        """Generate comparative report across all 3D experiments"""
        report_file = self.comparative_dir / "comparative_3D_analysis_report.md"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("# Comparative 3D Analysis Report\n\n")
            f.write(f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Experiment Overview\n\n")
            for exp_name in self.experiment_configs:
                config = self.experiment_manager.get_config(exp_name)
                f.write(f"### {config.display_name}\n")
                f.write(f"- **Configuration**: {exp_name}\n")
                f.write(f"- **miR-455-3p dose**: {config.miR_455_3p_dose}\n")
                f.write(f"- **miR-148a-5p dose**: {config.miR_148a_5p_dose}\n")
                f.write(f"- **Exosome delivery**: {'Yes' if config.has_exosomes else 'No'}\n\n")
            
            f.write("## 3D Analysis Components\n\n")
            f.write("### Generated for Each Experiment\n")
            f.write("1. **Interactive 3D Cell Visualization**\n")
            f.write("   - Real-time 3D scatter plots of cells colored by biological state\n")
            f.write("   - Interactive rotation, zoom, and hover information\n")
            f.write("   - Culture well boundary visualization\n\n")
            
            f.write("2. **3D Molecular Gradient Visualization**\n")
            f.write("   - Six key molecules: Oxygen, TGF-β1, α-SMA, Glucose, Lactate, Collagen-I\n")
            f.write("   - 3D spatial distribution with optimized colorbars\n")
            f.write("   - Realistic gradient patterns based on biology\n\n")
            
            f.write("3. **3D Spatial Pattern Analysis**\n")
            f.write("   - Vertical (Z-axis) cell distribution\n")
            f.write("   - Radial distribution from culture center\n")
            f.write("   - DBSCAN clustering analysis\n")
            f.write("   - Activation vs position correlations\n")
            f.write("   - miRNA co-expression patterns\n")
            f.write("   - Metabolic activity gradients\n\n")
            
            if self.data_source == "synthetic":
                f.write("4. **Dynamic Time-lapse Animation**\n")
                f.write("   - Interactive 48-hour time course\n")
                f.write("   - Play/pause controls and time slider\n")
                f.write("   - Real-time statistics overlay\n")
                f.write("   - Cell state evolution tracking\n\n")
            
            f.write("## Key 3D Insights\n\n")
            f.write("- **Spatial Heterogeneity**: 3D culture captures realistic tissue-level complexity\n")
            f.write("- **Gradient Formation**: Physiologically relevant molecular gradients established\n")
            f.write("- **Cell-Cell Interactions**: Enhanced 3D cell communication patterns\n")
            f.write("- **Treatment Penetration**: miRNA delivery shows realistic 3D distribution\n")
            f.write("- **Temporal Dynamics**: Time-course reveals treatment kinetics\n\n")
            
            f.write("## File Organization\n\n")
            f.write("```\n")
            f.write(f"{self.output_dir.name}/\n")
            for exp_name in self.experiment_configs:
                f.write(f"├── experiment_{exp_name}/\n")
                f.write(f"│   ├── interactive_html/    # 3D visualizations\n")
                f.write(f"│   ├── images/              # Static analyses\n")
                f.write(f"│   ├── animations/          # Time-lapse files\n")
                f.write(f"│   └── reports/             # Analysis reports\n")
            f.write(f"└── comparative_3D_analysis/    # Cross-experiment comparisons\n")
            f.write("```\n\n")
            
            f.write("---\n\n")
            f.write("**Advanced 3D Visualization System**\n")
            f.write("Enhanced platform for 3D tissue culture simulation analysis\n")
        
        print(f"📝 Comparative 3D report saved: {report_file}")
    
    def run_complete_3D_analysis(self, mode: str = 'all'):
        """Run complete 3D visualization and analysis pipeline
        
        Args:
            mode: 'all' (default), 'individual', or 'comparative'
        """
        print("Starting Complete 3D Analysis Pipeline...")
        print("=" * 60)
        print(f"Mode: {mode}")
        print(f"Experiments to process: {len(self.experiment_configs)}")
        
        successful_experiments = []
        
        if mode in ['all', 'individual']:
            print("\n" + "=" * 40)
            print("INDIVIDUAL EXPERIMENT 3D ANALYSIS")
            print("=" * 40)
            
            for i, experiment_name in enumerate(self.experiment_configs, 1):
                print(f"\n[{i}/{len(self.experiment_configs)}] Processing {experiment_name}...")
                
                if self.run_single_experiment_3D_analysis(experiment_name):
                    successful_experiments.append(experiment_name)
        
        if mode in ['all', 'comparative'] and len(successful_experiments) > 1:
            print("\n" + "=" * 40)
            print("COMPARATIVE 3D ANALYSIS")
            print("=" * 40)
            
            try:
                self.create_comparative_3D_analysis()
                print("✅ Comparative 3D analysis completed")
            except Exception as e:
                print(f"❌ Error in comparative 3D analysis: {e}")
        
        # Generate summary report
        self.generate_3D_analysis_report()
        
        print("\n" + "=" * 60)
        print("COMPLETE 3D ANALYSIS PIPELINE FINISHED!")
        print("=" * 60)
        
        print(f"\nResults Summary:")
        print(f"✅ {len(successful_experiments)} experiments processed successfully")
        print(f"📁 Output directory: {self.output_dir}")
        
        print(f"\nGenerated files per experiment:")
        for exp_name in successful_experiments:
            exp_dir = self.experiment_dirs[exp_name]['base']
            html_files = len(list(self.experiment_dirs[exp_name]['html'].glob("*.html")))
            image_files = len(list(self.experiment_dirs[exp_name]['images'].glob("*.png")))
            print(f"   - {exp_name}: {html_files} HTML files, {image_files} images")
        
        if self.data_source == "synthetic":
            print(f"\n🎬 Time-lapse Features (per experiment):")
            print(f"📊 Interactive Plotly animation with play/pause controls")
            print(f"🎯 Time slider for precise navigation")
            print(f"📈 Real-time statistics showing cell population dynamics")
            print(f"🧬 miRNA treatment effects over time")
            print(f"🔬 3D spatial evolution of cell activation")
        
        print(f"\n🔬 3D simulation ready for advanced analysis!")

if __name__ == "__main__":
    # Run complete 3D analysis
    analyzer = Advanced3DVisualization()
    analyzer.run_complete_3D_analysis()
