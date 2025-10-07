#!/usr/bin/env python3
"""
Advanced Visualization Suite for Liver Fibrosis Simulation
Professional publication-quality visualizations including 3D heatmaps and animated GIFs
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

# Set high-quality defaults
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10

class AdvancedVisualization:
    def __init__(self, base_dir="."):
        """Initialize advanced visualization system"""
        self.base_dir = Path(base_dir)
        self.output_dir = self.base_dir / "visualization_outputs"
        self.output_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        self.figures_dir = self.output_dir / "figures"
        self.pdfs_dir = self.output_dir / "pdfs"
        self.animations_dir = self.output_dir / "animations"
        self.interactive_dir = self.output_dir / "interactive"
        
        for dir_path in [self.figures_dir, self.pdfs_dir, self.animations_dir, self.interactive_dir]:
            dir_path.mkdir(exist_ok=True)
        
        # Check for data sources
        self.synthetic_dir = self.base_dir / "synthetic_results"
        if self.synthetic_dir.exists():
            print("Using synthetic data for advanced visualization")
            self.data_source = "synthetic"
        else:
            print("Using PhysiCell output data for advanced visualization")
            self.data_source = "physicell"
        
        # Professional color schemes
        self.cell_colors = {
            'quiescent': '#2E86AB',      # Professional blue
            'activated': '#A23B72',      # Professional magenta
            'proliferating': '#F18F01',  # Professional orange
            'senescent': '#C73E1D',      # Professional red
            'apoptotic': '#592720'       # Professional dark red
        }
        
        # Load data
        self.load_simulation_data()
    
    def load_simulation_data(self):
        """Load simulation data from available sources"""
        self.cells_data = {}
        self.microenv_data = {}
        self.time_series_data = {}
        
        if self.data_source == "synthetic":
            self._load_synthetic_data()
        else:
            self._load_physicell_data()
    
    def _load_synthetic_data(self):
        """Load synthetic data"""
        print("Loading synthetic data for advanced visualization...")
        
        # Load time series data for all groups
        for group_file in self.synthetic_dir.glob("time_series_*.csv"):
            group_name = group_file.stem.replace("time_series_", "")
            self.time_series_data[group_name] = pd.read_csv(group_file)
        
        # Load detailed cell and microenvironment data
        sample_group = "dual_miRNA_1_1"  # Use representative group
        sample_dir = self.synthetic_dir / f"output_{sample_group}"
        
        if sample_dir.exists():
            # Load multiple time points
            for time_point in [0, 720, 1440, 2160, 2880]:  # 0, 12h, 24h, 36h, 48h
                cells_file = sample_dir / f"cells_{time_point:04d}.csv"
                microenv_file = sample_dir / f"microenv_{time_point:04d}.csv"
                
                if cells_file.exists():
                    self.cells_data[time_point] = pd.read_csv(cells_file)
                if microenv_file.exists():
                    self.microenv_data[time_point] = pd.read_csv(microenv_file)
        
        print(f"Loaded data for {len(self.time_series_data)} experimental groups")
        print(f"Loaded detailed data for {len(self.cells_data)} time points")
    
    def _load_physicell_data(self):
        """Load standard PhysiCell data"""
        # Implementation for PhysiCell XML data loading
        print("PhysiCell data loading not yet implemented")
        
        # Generate demo data for now
        self._generate_demo_data()
    
    def _generate_demo_data(self):
        """Generate demonstration data"""
        print("Generating demonstration data...")
        
        # Create time series data
        time_points = np.linspace(0, 2880, 49)
        groups = ['control', 'positive_model', 'miR_455_3p', 'miR_148a_5p', 'dual_miRNA_1_1']
        
        for group in groups:
            data = []
            for t in time_points:
                # Generate realistic time-dependent data
                base_activation = 0.1 if group == 'control' else 0.7
                if 'miR' in group:
                    treatment_effect = 1 - 0.6 * (1 - np.exp(-t/1440))
                    if 'dual' in group:
                        treatment_effect *= 1.3  # Synergy effect
                    activation = base_activation * treatment_effect
                else:
                    activation = base_activation
                
                data.append({
                    'time': t,
                    'avg_activation': activation + np.random.normal(0, 0.02),
                    'total_cells': 120 + int(np.random.normal(0, 5)),
                    'avg_miR455': 0.4 if '455' in group or 'dual' in group else 0.0,
                    'avg_miR148': 0.4 if '148' in group or 'dual' in group else 0.0,
                    'treatment_efficacy': (1 - activation/base_activation) * 100 if base_activation > 0 else 0
                })
            
            self.time_series_data[group] = pd.DataFrame(data)
        
        # Generate spatial cell data
        for time_point in [0, 720, 1440, 2160, 2880]:
            n_cells = 120
            x = np.random.uniform(-400, 400, n_cells)
            y = np.random.uniform(-400, 400, n_cells)
            z = np.random.uniform(-200, 200, n_cells)
            
            # Generate cell states based on time
            activation_levels = np.random.beta(2, 5, n_cells) * (1 + time_point/2880)
            activation_levels = np.clip(activation_levels, 0, 1)
            
            self.cells_data[time_point] = pd.DataFrame({
                'x': x, 'y': y, 'z': z,
                'activation_level': activation_levels,
                'miR_455_3p_level': np.random.exponential(0.2, n_cells),
                'miR_148a_5p_level': np.random.exponential(0.2, n_cells),
                'stress_level': np.random.beta(2, 3, n_cells),
                'cell_cycle_phase': np.random.randint(0, 4, n_cells)
            })
            
            # Generate microenvironment data
            n_voxels = 1000
            x_vox = np.random.uniform(-400, 400, n_voxels)
            y_vox = np.random.uniform(-400, 400, n_voxels)
            z_vox = np.random.uniform(-200, 200, n_voxels)
            
            self.microenv_data[time_point] = pd.DataFrame({
                'x': x_vox, 'y': y_vox, 'z': z_vox,
                'TGF_beta1': np.random.exponential(5, n_voxels),
                'oxygen': np.random.normal(0.8, 0.1, n_voxels),
                'glucose': np.random.normal(0.7, 0.1, n_voxels),
                'alpha_SMA': np.random.exponential(2, n_voxels),
                'collagen_I': np.random.exponential(1.5, n_voxels)
            })
    
    def create_3d_heatmaps(self):
        """Create professional 3D heatmaps of molecular concentrations"""
        print("Creating 3D molecular distribution heatmaps...")
        
        # Select final time point
        final_time = max(self.microenv_data.keys())
        microenv_data = self.microenv_data[final_time]
        
        molecules = ['TGF_beta1', 'oxygen', 'glucose', 'alpha_SMA']
        molecule_titles = {
            'TGF_beta1': 'TGF-β1 Concentration',
            'oxygen': 'Oxygen Distribution', 
            'glucose': 'Glucose Distribution',
            'alpha_SMA': 'α-SMA Expression'
        }
        
        # Create 2x2 subplot for multiple molecules
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[molecule_titles[mol] for mol in molecules],
            specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}],
                   [{'type': 'scatter3d'}, {'type': 'scatter3d'}]]
        )
        
        positions = [(1,1), (1,2), (2,1), (2,2)]
        
        for i, molecule in enumerate(molecules):
            row, col = positions[i]
            
            # Create 3D scatter plot with molecular concentrations
            values = microenv_data[molecule].values
            
            scatter = go.Scatter3d(
                x=microenv_data['x'],
                y=microenv_data['y'], 
                z=microenv_data['z'],
                mode='markers',
                marker=dict(
                    size=3,
                    color=values,
                    colorscale='Viridis',
                    showscale=True if i == 0 else False,
                    colorbar=dict(title=f"{molecule_titles[molecule]} (au)", x=0.02 + i*0.25, len=0.4),
                    opacity=0.8
                ),
                name=molecule_titles[molecule],
                showlegend=False
            )
            
            fig.add_trace(scatter, row=row, col=col)
        
        # Update layout
        fig.update_layout(
            title=dict(
                text="Three-Dimensional Molecular Distribution Analysis",
                x=0.5,
                font=dict(size=16)
            ),
            scene=dict(
                xaxis_title="X Position (μm)",
                yaxis_title="Y Position (μm)",
                zaxis_title="Z Position (μm)",
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            scene2=dict(
                xaxis_title="X Position (μm)",
                yaxis_title="Y Position (μm)", 
                zaxis_title="Z Position (μm)",
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            scene3=dict(
                xaxis_title="X Position (μm)",
                yaxis_title="Y Position (μm)",
                zaxis_title="Z Position (μm)",
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            scene4=dict(
                xaxis_title="X Position (μm)",
                yaxis_title="Y Position (μm)",
                zaxis_title="Z Position (μm)",
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            width=1200,
            height=900,
            margin=dict(l=0, r=0, t=60, b=0)
        )
        
        # Save interactive version
        interactive_file = self.interactive_dir / "3D_molecular_heatmaps.html"
        fig.write_html(str(interactive_file))
        
        # Save static high-resolution image
        img_file = self.figures_dir / "3D_molecular_heatmaps.png"
        fig.write_image(str(img_file), width=1200, height=900, scale=2)
        
        # Save PDF
        pdf_file = self.pdfs_dir / "3D_molecular_heatmaps.pdf"
        fig.write_image(str(pdf_file), width=1200, height=900, format='pdf')
        
        print(f"3D heatmaps saved to: {interactive_file}")
        print(f"High-resolution PNG: {img_file}")
        print(f"PDF version: {pdf_file}")
    
    def create_2d_cell_dynamics_gif(self):
        """Create 2D animated GIF of cell dynamics over time"""
        print("Creating 2D cell dynamics animation...")
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Prepare data for all time points
        time_points = sorted(self.cells_data.keys())
        
        def animate(frame_idx):
            ax.clear()
            
            time_point = time_points[frame_idx]
            cells_df = self.cells_data[time_point]
            
            # Color cells by activation level
            scatter = ax.scatter(
                cells_df['x'], cells_df['y'],
                c=cells_df['activation_level'],
                cmap='RdYlBu_r',
                s=60,
                alpha=0.8,
                edgecolors='black',
                linewidth=0.5
            )
            
            # Styling
            ax.set_xlim(-450, 450)
            ax.set_ylim(-450, 450)
            ax.set_xlabel('X Position (μm)', fontsize=12)
            ax.set_ylabel('Y Position (μm)', fontsize=12)
            ax.set_title(f'Hepatic Stellate Cell Activation Dynamics\nTime: {time_point/60:.1f} hours', 
                        fontsize=14, fontweight='bold')
            
            # Add colorbar
            if frame_idx == 0:
                cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label('Activation Level', fontsize=12)
            
            # Add time annotation
            ax.text(0.02, 0.98, f't = {time_point/60:.1f} h', 
                   transform=ax.transAxes, fontsize=12, 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                   verticalalignment='top')
            
            # Add statistics
            mean_activation = cells_df['activation_level'].mean()
            std_activation = cells_df['activation_level'].std()
            ax.text(0.02, 0.85, f'Mean Activation: {mean_activation:.3f}±{std_activation:.3f}', 
                   transform=ax.transAxes, fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
        
        # Create animation
        anim = animation.FuncAnimation(
            fig, animate, frames=len(time_points), 
            interval=800, blit=False, repeat=True
        )
        
        # Save as GIF
        gif_file = self.animations_dir / "cell_activation_dynamics_2D.gif"
        anim.save(str(gif_file), writer='pillow', fps=1.25, dpi=150)
        
        # Save individual frames as high-res images
        for i, time_point in enumerate(time_points):
            animate(i)
            frame_file = self.figures_dir / f"cell_dynamics_2D_t{time_point:04d}.png"
            plt.savefig(str(frame_file), dpi=300, bbox_inches='tight')
            
            # Save PDF version
            pdf_file = self.pdfs_dir / f"cell_dynamics_2D_t{time_point:04d}.pdf" 
            plt.savefig(str(pdf_file), format='pdf', bbox_inches='tight')
        
        plt.close()
        print(f"2D animation saved: {gif_file}")
        print(f"Individual frames saved to: {self.figures_dir}")
    
    def create_3d_cell_dynamics_gif(self):
        """Create 3D animated GIF of cell dynamics"""
        print("Creating 3D cell dynamics animation...")
        
        # Create 3D visualization using plotly for better quality
        time_points = sorted(self.cells_data.keys())
        
        frames = []
        for time_point in time_points:
            cells_df = self.cells_data[time_point]
            
            frame = go.Frame(
                data=[
                    go.Scatter3d(
                        x=cells_df['x'],
                        y=cells_df['y'],
                        z=cells_df['z'],
                        mode='markers',
                        marker=dict(
                            size=6,
                            color=cells_df['activation_level'],
                            colorscale='RdYlBu_r',
                            showscale=True,
                            colorbar=dict(title="Activation Level", x=0.02),
                            opacity=0.8,
                            line=dict(width=0.5, color='black')
                        ),
                        name=f't = {time_point/60:.1f}h'
                    )
                ],
                name=str(time_point)
            )
            frames.append(frame)
        
        # Initial frame
        initial_cells = self.cells_data[time_points[0]]
        
        fig = go.Figure(
            data=[
                go.Scatter3d(
                    x=initial_cells['x'],
                    y=initial_cells['y'],
                    z=initial_cells['z'],
                    mode='markers',
                    marker=dict(
                        size=6,
                        color=initial_cells['activation_level'],
                        colorscale='RdYlBu_r',
                        showscale=True,
                        colorbar=dict(title="Activation Level"),
                        opacity=0.8,
                        line=dict(width=0.5, color='black')
                    )
                )
            ],
            frames=frames
        )
        
        # Add animation controls
        fig.update_layout(
            title=dict(
                text="Three-Dimensional Cell Activation Dynamics",
                x=0.5,
                font=dict(size=16)
            ),
            scene=dict(
                xaxis_title="X Position (μm)",
                yaxis_title="Y Position (μm)",
                zaxis_title="Z Position (μm)",
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)),
                xaxis=dict(range=[-450, 450]),
                yaxis=dict(range=[-450, 450]),
                zaxis=dict(range=[-250, 250])
            ),
            width=900,
            height=700,
            updatemenus=[{
                'type': 'buttons',
                'showactive': False,
                'buttons': [
                    {
                        'label': 'Play',
                        'method': 'animate',
                        'args': [None, {'frame': {'duration': 800, 'redraw': True},
                                       'fromcurrent': True, 'transition': {'duration': 300}}]
                    },
                    {
                        'label': 'Pause',
                        'method': 'animate',
                        'args': [[None], {'frame': {'duration': 0, 'redraw': False},
                                         'mode': 'immediate',
                                         'transition': {'duration': 0}}]
                    }
                ]
            }],
            sliders=[{
                'steps': [
                    {
                        'args': [[str(tp)], {'frame': {'duration': 300, 'redraw': True},
                                            'mode': 'immediate',
                                            'transition': {'duration': 300}}],
                        'label': f'{tp/60:.1f}h',
                        'method': 'animate'
                    }
                    for tp in time_points
                ],
                'active': 0,
                'currentvalue': {'prefix': 'Time: '},
                'len': 0.9,
                'x': 0.05,
                'xanchor': 'left',
                'y': 0,
                'yanchor': 'top'
            }]
        )
        
        # Save interactive version
        interactive_file = self.interactive_dir / "3D_cell_dynamics_interactive.html"
        fig.write_html(str(interactive_file))
        
        # Save static frames
        for i, time_point in enumerate(time_points):
            # Update to specific frame
            fig.update_traces(
                x=self.cells_data[time_point]['x'],
                y=self.cells_data[time_point]['y'],
                z=self.cells_data[time_point]['z'],
                marker_color=self.cells_data[time_point]['activation_level']
            )
            
            frame_file = self.figures_dir / f"cell_dynamics_3D_t{time_point:04d}.png"
            fig.write_image(str(frame_file), width=900, height=700, scale=2)
            
            pdf_file = self.pdfs_dir / f"cell_dynamics_3D_t{time_point:04d}.pdf"
            fig.write_image(str(pdf_file), width=900, height=700, format='pdf')
        
        print(f"3D interactive animation saved: {interactive_file}")
        print(f"3D static frames saved to: {self.figures_dir}")
    
    def create_treatment_efficacy_analysis(self):
        """Create comprehensive treatment efficacy analysis"""
        print("Creating treatment efficacy analysis...")
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Therapeutic Efficacy Analysis in Hepatic Stellate Cell Model', 
                    fontsize=16, fontweight='bold')
        
        # Prepare data
        groups = list(self.time_series_data.keys())
        colors = plt.cm.Set2(np.linspace(0, 1, len(groups)))
        
        # Plot 1: Time course of activation
        ax = axes[0, 0]
        for i, group in enumerate(groups):
            data = self.time_series_data[group]
            time_hours = data['time'] / 60
            ax.plot(time_hours, data['avg_activation'], 
                   label=group.replace('_', ' ').title(), 
                   color=colors[i], linewidth=2.5, alpha=0.8)
        
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Mean Activation Level')
        ax.set_title('Temporal Activation Dynamics')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Final efficacy comparison
        ax = axes[0, 1]
        final_efficacies = []
        group_names = []
        for group in groups:
            data = self.time_series_data[group]
            if 'treatment_efficacy' in data.columns:
                final_efficacy = data['treatment_efficacy'].iloc[-1]
            else:
                # Calculate based on activation reduction
                final_activation = data['avg_activation'].iloc[-1]
                control_activation = self.time_series_data.get('control', data)['avg_activation'].iloc[-1]
                final_efficacy = max(0, (control_activation - final_activation) / control_activation * 100)
            
            final_efficacies.append(final_efficacy)
            group_names.append(group.replace('_', ' ').title())
        
        bars = ax.bar(range(len(group_names)), final_efficacies, 
                     color=colors[:len(groups)], alpha=0.8, edgecolor='black')
        ax.set_xlabel('Treatment Group')
        ax.set_ylabel('Therapeutic Efficacy (%)')
        ax.set_title('Comparative Treatment Efficacy')
        ax.set_xticks(range(len(group_names)))
        ax.set_xticklabels(group_names, rotation=45, ha='right')
        
        # Add value labels on bars
        for bar, value in zip(bars, final_efficacies):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{value:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Plot 3: miRNA levels over time
        ax = axes[1, 0]
        mirna_groups = [g for g in groups if 'miR' in g]
        for group in mirna_groups:
            data = self.time_series_data[group] 
            time_hours = data['time'] / 60
            if 'avg_miR455' in data.columns and data['avg_miR455'].max() > 0.01:
                ax.plot(time_hours, data['avg_miR455'], 
                       label=f"{group} (miR-455-3p)", linestyle='-', linewidth=2)
            if 'avg_miR148' in data.columns and data['avg_miR148'].max() > 0.01:
                ax.plot(time_hours, data['avg_miR148'], 
                       label=f"{group} (miR-148a-5p)", linestyle='--', linewidth=2)
        
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('miRNA Expression Level')
        ax.set_title('miRNA Expression Dynamics')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 4: Dose-response relationship
        ax = axes[1, 1]
        if len(mirna_groups) >= 3:  # Need multiple doses for dose-response
            doses = []
            responses = []
            for group in mirna_groups:
                data = self.time_series_data[group]
                # Calculate "dose" as sum of miRNA levels
                dose = data['avg_miR455'].max() + data['avg_miR148'].max()
                response = final_efficacies[groups.index(group)]
                doses.append(dose)
                responses.append(response)
            
            # Sort by dose
            sorted_data = sorted(zip(doses, responses))
            doses, responses = zip(*sorted_data)
            
            ax.plot(doses, responses, 'o-', linewidth=2, markersize=8, color='red')
            ax.set_xlabel('Combined miRNA Level')
            ax.set_ylabel('Therapeutic Efficacy (%)')
            ax.set_title('Dose-Response Relationship')
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'Insufficient data\nfor dose-response\nanalysis', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=12)
        
        plt.tight_layout()
        
        # Save high-resolution versions
        png_file = self.figures_dir / "treatment_efficacy_analysis.png"
        pdf_file = self.pdfs_dir / "treatment_efficacy_analysis.pdf"
        
        plt.savefig(str(png_file), dpi=300, bbox_inches='tight')
        plt.savefig(str(pdf_file), format='pdf', bbox_inches='tight')
        plt.close()
        
        print(f"Treatment efficacy analysis saved: {png_file}")
        print(f"PDF version: {pdf_file}")
    
    def create_morphological_analysis(self):
        """Create advanced morphological analysis"""
        print("Creating morphological analysis...")
        
        # Select final time point
        final_time = max(self.cells_data.keys())
        cells_df = self.cells_data[final_time]
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Cellular Morphological Analysis', fontsize=16, fontweight='bold')
        
        # Cell distribution by activation
        ax = axes[0, 0]
        ax.hist(cells_df['activation_level'], bins=20, alpha=0.7, 
               color='steelblue', edgecolor='black')
        ax.set_xlabel('Activation Level')
        ax.set_ylabel('Cell Count')
        ax.set_title('Activation Level Distribution')
        ax.grid(True, alpha=0.3)
        
        # Spatial distribution colored by activation
        ax = axes[0, 1]
        scatter = ax.scatter(cells_df['x'], cells_df['y'], 
                           c=cells_df['activation_level'], 
                           cmap='RdYlBu_r', s=40, alpha=0.7, edgecolors='black', linewidth=0.5)
        ax.set_xlabel('X Position (μm)')
        ax.set_ylabel('Y Position (μm)')
        ax.set_title('Spatial Activation Pattern')
        plt.colorbar(scatter, ax=ax, label='Activation Level')
        
        # Stress vs Activation correlation
        ax = axes[0, 2]
        if 'stress_level' in cells_df.columns:
            ax.scatter(cells_df['activation_level'], cells_df['stress_level'], 
                      alpha=0.6, s=30, color='red')
            ax.set_xlabel('Activation Level')
            ax.set_ylabel('Stress Level')
            ax.set_title('Activation-Stress Correlation')
            
            # Add correlation coefficient
            corr = cells_df['activation_level'].corr(cells_df['stress_level'])
            ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes, 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            ax.text(0.5, 0.5, 'Stress data\nnot available', 
                   ha='center', va='center', transform=ax.transAxes)
        ax.grid(True, alpha=0.3)
        
        # Cell cycle distribution
        ax = axes[1, 0]
        if 'cell_cycle_phase' in cells_df.columns:
            cycle_counts = cells_df['cell_cycle_phase'].value_counts().sort_index()
            cycle_labels = ['G0/G1', 'S', 'G2', 'M'][:len(cycle_counts)]
            ax.pie(cycle_counts.values, labels=cycle_labels, autopct='%1.1f%%', 
                  colors=['lightblue', 'lightgreen', 'orange', 'pink'])
            ax.set_title('Cell Cycle Distribution')
        else:
            ax.text(0.5, 0.5, 'Cell cycle data\nnot available', 
                   ha='center', va='center', transform=ax.transAxes)
        
        # miRNA levels comparison
        ax = axes[1, 1]
        if 'miR_455_3p_level' in cells_df.columns and 'miR_148a_5p_level' in cells_df.columns:
            ax.scatter(cells_df['miR_455_3p_level'], cells_df['miR_148a_5p_level'], 
                      alpha=0.6, s=30, color='purple')
            ax.set_xlabel('miR-455-3p Level')
            ax.set_ylabel('miR-148a-5p Level')
            ax.set_title('miRNA Expression Correlation')
            
            # Add trend line
            z = np.polyfit(cells_df['miR_455_3p_level'], cells_df['miR_148a_5p_level'], 1)
            p = np.poly1d(z)
            ax.plot(cells_df['miR_455_3p_level'].sort_values(), 
                   p(cells_df['miR_455_3p_level'].sort_values()), "r--", alpha=0.8)
        else:
            ax.text(0.5, 0.5, 'miRNA data\nnot available', 
                   ha='center', va='center', transform=ax.transAxes)
        ax.grid(True, alpha=0.3)
        
        # Summary statistics
        ax = axes[1, 2]
        ax.axis('off')
        
        stats_text = f"""
Morphological Summary
{'='*20}
Total Cells: {len(cells_df)}
Mean Activation: {cells_df['activation_level'].mean():.3f}
Std Activation: {cells_df['activation_level'].std():.3f}
Max Activation: {cells_df['activation_level'].max():.3f}
Min Activation: {cells_df['activation_level'].min():.3f}

Spatial Statistics
{'='*15}
X Range: {cells_df['x'].max()-cells_df['x'].min():.1f} μm
Y Range: {cells_df['y'].max()-cells_df['y'].min():.1f} μm
Z Range: {cells_df['z'].max()-cells_df['z'].min():.1f} μm

Analysis Time Point: {final_time/60:.1f} hours
"""
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
               verticalalignment='top', fontfamily='monospace',
               bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        
        # Save files
        png_file = self.figures_dir / "morphological_analysis.png"
        pdf_file = self.pdfs_dir / "morphological_analysis.pdf"
        
        plt.savefig(str(png_file), dpi=300, bbox_inches='tight')
        plt.savefig(str(pdf_file), format='pdf', bbox_inches='tight')
        plt.close()
        
        print(f"Morphological analysis saved: {png_file}")
        print(f"PDF version: {pdf_file}")
    
    def run_complete_advanced_visualization(self):
        """Run complete advanced visualization suite"""
        print("Starting Advanced Visualization Suite")
        print("=" * 50)
        
        try:
            # Create all visualizations
            self.create_3d_heatmaps()
            self.create_2d_cell_dynamics_gif()
            self.create_3d_cell_dynamics_gif() 
            self.create_treatment_efficacy_analysis()
            self.create_morphological_analysis()
            
            print("\n" + "=" * 50)
            print("Advanced Visualization Suite Completed!")
            print("=" * 50)
            
            print(f"\nGenerated files:")
            print(f"📁 Figures: {self.figures_dir}")
            print(f"📁 PDFs: {self.pdfs_dir}")
            print(f"📁 Animations: {self.animations_dir}")
            print(f"📁 Interactive: {self.interactive_dir}")
            
            print(f"\nFile counts:")
            print(f"   PNG files: {len(list(self.figures_dir.glob('*.png')))}")
            print(f"   PDF files: {len(list(self.pdfs_dir.glob('*.pdf')))}")
            print(f"   GIF files: {len(list(self.animations_dir.glob('*.gif')))}")
            print(f"   HTML files: {len(list(self.interactive_dir.glob('*.html')))}")
            
            return True
            
        except Exception as e:
            print(f"ERROR: Advanced visualization failed: {e}")
            return False

if __name__ == "__main__":
    # Run advanced visualization suite
    visualizer = AdvancedVisualization()
    visualizer.run_complete_advanced_visualization()
