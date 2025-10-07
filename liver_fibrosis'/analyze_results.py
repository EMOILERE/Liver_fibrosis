#!/usr/bin/env python3
"""
iGEM Liver Fibrosis Simulation Results Analysis Script
Analyzes PhysiCell simulation output and generates synergy reports
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
import xml.etree.ElementTree as ET

# Set font for plots
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False

class LiverFibrosisAnalyzer:
 def __init__(self, base_dir="."):
 """Initialize analyzer"""
 self.base_dir = Path(base_dir)
 self.experimental_groups = [
 'control', 'positive_model', 'natural_exosomes', 'NC_mimic',
 'miR_455_3p', 'miR_148a_5p', 'dual_miRNA_1_1', 'dual_miRNA_1_2', 'dual_miRNA_2_1'
 ]
 self.results = {}
 
 def load_simulation_data(self):
 """Load simulation data for all experimental groups"""
 print("Loading simulation data...")
 
 for group in self.experimental_groups:
 output_dir = self.base_dir / f"output_{group}"
 if not output_dir.exists():
 print(f"Warning: {output_dir} not found, skipping {group}")
 continue
 
 self.results[group] = self._parse_group_data(output_dir)
 print(f"Loaded data for {group}")
 
 def _parse_group_data(self, output_dir):
 """Parse data for a single experimental group"""
 group_data = {
 'time_points': [],
 'total_cells': [],
 'activated_cells': [],
 'activation_fraction': [],
 'alpha_SMA_levels': [],
 'collagen_I_levels': [],
 'miR_455_3p_levels': [],
 'miR_148a_5p_levels': []
 }
 
 # Find all output files
 output_files = sorted(glob.glob(str(output_dir / "output*.xml")))
 
 for output_file in output_files:
 try:
 time_point, cell_data = self._parse_output_file(output_file)
 if time_point is not None:
 group_data['time_points'].append(time_point)
 group_data['total_cells'].append(cell_data['total_cells'])
 group_data['activated_cells'].append(cell_data['activated_cells'])
 group_data['activation_fraction'].append(
 cell_data['activated_cells'] / max(cell_data['total_cells'], 1)
 )
 group_data['alpha_SMA_levels'].append(cell_data['alpha_SMA_avg'])
 group_data['collagen_I_levels'].append(cell_data['collagen_I_avg'])
 group_data['miR_455_3p_levels'].append(cell_data['miR_455_3p_avg'])
 group_data['miR_148a_5p_levels'].append(cell_data['miR_148a_5p_avg'])
 except Exception as e:
 print(f"Error parsing {output_file}: {e}")
 continue
 
 return group_data
 
 def _parse_output_file(self, output_file):
 """Parse a single output file"""
 try:
 tree = ET.parse(output_file)
 root = tree.getroot()
 
 # Get time point
 time_element = root.find('.//current_time')
 if time_element is None:
 return None, None
 time_point = float(time_element.text)
 
 # Initialize cell data
 cell_data = {
 'total_cells': 0,
 'activated_cells': 0,
 'alpha_SMA_avg': 0.0,
 'collagen_I_avg': 0.0,
 'miR_455_3p_avg': 0.0,
 'miR_148a_5p_avg': 0.0
 }
 
 # Parse cell data
 cells = root.findall('.//cell')
 alpha_SMA_values = []
 collagen_I_values = []
 miR_455_values = []
 miR_148_values = []
 
 for cell in cells:
 cell_data['total_cells'] += 1
 
 # Check cell type
 cell_type = cell.find('.//cell_type')
 if cell_type is not None and cell_type.text == '1': # LX2_activated
 cell_data['activated_cells'] += 1
 
 # Get custom data
 custom_vars = cell.findall('.//custom_variables/custom_variable')
 for var in custom_vars:
 name = var.get('name')
 value = float(var.text) if var.text else 0.0
 
 if name == 'alpha_SMA_production':
 alpha_SMA_values.append(value)
 elif name == 'collagen_I_production':
 collagen_I_values.append(value)
 elif name == 'miR_455_3p_level':
 miR_455_values.append(value)
 elif name == 'miR_148a_5p_level':
 miR_148_values.append(value)
 
 # Calculate averages
 cell_data['alpha_SMA_avg'] = np.mean(alpha_SMA_values) if alpha_SMA_values else 0.0
 cell_data['collagen_I_avg'] = np.mean(collagen_I_values) if collagen_I_values else 0.0
 cell_data['miR_455_3p_avg'] = np.mean(miR_455_values) if miR_455_values else 0.0
 cell_data['miR_148a_5p_avg'] = np.mean(miR_148_values) if miR_148_values else 0.0
 
 return time_point, cell_data
 
 except Exception as e:
 print(f"Error parsing file {output_file}: {e}")
 return None, None
 
 def calculate_synergy_effects(self):
 """Calculate synergistic enhancement effects"""
 print("Calculating synergy effects...")
 
 synergy_analysis = {}
 
 # Get endpoint data (last time point)
 endpoint_data = {}
 for group in self.results:
 if self.results[group]['time_points']:
 last_idx = -1
 endpoint_data[group] = {
 'activation_fraction': self.results[group]['activation_fraction'][last_idx],
 'alpha_SMA': self.results[group]['alpha_SMA_levels'][last_idx],
 'collagen_I': self.results[group]['collagen_I_levels'][last_idx]
 }
 
 # Calculate inhibition efficiency
 if 'positive_model' in endpoint_data:
 baseline = endpoint_data['positive_model']
 print(f"Baseline (positive_model) data: {baseline}")
 
 for group in endpoint_data:
 if group == 'positive_model':
 continue
 
 data = endpoint_data[group]
 
 # Safe division calculation to avoid division by zero
 def safe_inhibition_calc(treatment_value, baseline_value):
 if baseline_value == 0:
 if treatment_value == 0:
 return 0 # Both are 0, no change
 else:
 return 0 # Baseline is 0 but treatment has value, cannot calculate inhibition rate
 else:
 return max(0, 1 - treatment_value / baseline_value)
 
 synergy_analysis[group] = {
 'activation_inhibition': safe_inhibition_calc(data['activation_fraction'], baseline['activation_fraction']),
 'alpha_SMA_inhibition': safe_inhibition_calc(data['alpha_SMA'], baseline['alpha_SMA']),
 'collagen_I_inhibition': safe_inhibition_calc(data['collagen_I'], baseline['collagen_I'])
 }
 
 # Analyze synergistic effects
 single_miR_groups = ['miR_455_3p', 'miR_148a_5p']
 dual_miR_groups = ['dual_miRNA_1_1', 'dual_miRNA_1_2', 'dual_miRNA_2_1']
 
 for dual_group in dual_miR_groups:
 if dual_group in synergy_analysis:
 # Calculate expected additive effects
 expected_effect = {}
 if all(g in synergy_analysis for g in single_miR_groups):
 for metric in ['activation_inhibition', 'alpha_SMA_inhibition', 'collagen_I_inhibition']:
 expected = (synergy_analysis['miR_455_3p'][metric] + 
 synergy_analysis['miR_148a_5p'][metric])
 observed = synergy_analysis[dual_group][metric]
 
 synergy_analysis[dual_group][f'{metric}_synergy'] = observed - expected
 synergy_analysis[dual_group][f'{metric}_synergy_fold'] = observed / max(expected, 0.001)
 
 self.synergy_results = synergy_analysis
 return synergy_analysis
 
 def plot_time_series(self):
 """Plot time series graphs"""
 print("Generating time series plots...")
 
 fig, axes = plt.subplots(2, 2, figsize=(15, 12))
 fig.suptitle('Liver Fibrosis Simulation Time Series Analysis', fontsize=16, fontweight='bold')
 
 # Define color scheme
 colors = {
 'control': 'gray',
 'positive_model': 'red',
 'natural_exosomes': 'lightcoral',
 'NC_mimic': 'pink',
 'miR_455_3p': 'lightblue',
 'miR_148a_5p': 'lightgreen',
 'dual_miRNA_1_1': 'blue',
 'dual_miRNA_1_2': 'green',
 'dual_miRNA_2_1': 'purple'
 }
 
 # Cell activation rate
 ax1 = axes[0, 0]
 for group in self.results:
 data = self.results[group]
 if data['time_points']:
 ax1.plot(data['time_points'], data['activation_fraction'], 
 label=group, color=colors.get(group, 'black'), linewidth=2)
 ax1.set_title('Cell Activation Rate Over Time')
 ax1.set_xlabel('Time (min)')
 ax1.set_ylabel('Activation Rate')
 ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
 ax1.grid(True, alpha=0.3)
 
 # α-SMA levels
 ax2 = axes[0, 1]
 for group in self.results:
 data = self.results[group]
 if data['time_points']:
 ax2.plot(data['time_points'], data['alpha_SMA_levels'], 
 label=group, color=colors.get(group, 'black'), linewidth=2)
 ax2.set_title('α-SMA Production Level')
 ax2.set_xlabel('Time (min)')
 ax2.set_ylabel('α-SMA Production Rate')
 ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
 ax2.grid(True, alpha=0.3)
 
 # Collagen I levels
 ax3 = axes[1, 0]
 for group in self.results:
 data = self.results[group]
 if data['time_points']:
 ax3.plot(data['time_points'], data['collagen_I_levels'], 
 label=group, color=colors.get(group, 'black'), linewidth=2)
 ax3.set_title('Collagen I Production Level')
 ax3.set_xlabel('Time (min)')
 ax3.set_ylabel('Collagen I Production Rate')
 ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
 ax3.grid(True, alpha=0.3)
 
 # miRNA levels
 ax4 = axes[1, 1]
 for group in ['miR_455_3p', 'miR_148a_5p', 'dual_miRNA_1_1']:
 if group in self.results:
 data = self.results[group]
 if data['time_points']:
 ax4.plot(data['time_points'], data['miR_455_3p_levels'], 
 label=f'{group}_miR455', color=colors.get(group, 'black'), 
 linewidth=2, linestyle='-')
 ax4.plot(data['time_points'], data['miR_148a_5p_levels'], 
 label=f'{group}_miR148', color=colors.get(group, 'black'), 
 linewidth=2, linestyle='--')
 ax4.set_title('Intracellular miRNA Levels')
 ax4.set_xlabel('Time (min)')
 ax4.set_ylabel('miRNA Concentration')
 ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
 ax4.grid(True, alpha=0.3)
 
 plt.tight_layout()
 plt.savefig('time_series_analysis.png', dpi=300, bbox_inches='tight')
 plt.show()
 
 def plot_synergy_analysis(self):
 """Plot synergy analysis graphs"""
 if not hasattr(self, 'synergy_results'):
 self.calculate_synergy_effects()
 
 print("Generating synergy analysis plots...")
 
 fig, axes = plt.subplots(2, 2, figsize=(15, 12))
 fig.suptitle('miRNA Synergistic Enhancement Analysis', fontsize=16, fontweight='bold')
 
 # Prepare data
 groups = []
 activation_inhibition = []
 alpha_SMA_inhibition = []
 collagen_I_inhibition = []
 synergy_scores = []
 
 for group in ['miR_455_3p', 'miR_148a_5p', 'dual_miRNA_1_1', 'dual_miRNA_1_2', 'dual_miRNA_2_1']:
 if group in self.synergy_results:
 groups.append(group)
 data = self.synergy_results[group]
 activation_inhibition.append(data.get('activation_inhibition', 0))
 alpha_SMA_inhibition.append(data.get('alpha_SMA_inhibition', 0))
 collagen_I_inhibition.append(data.get('collagen_I_inhibition', 0))
 
 # Calculate overall synergy score
 synergy_score = (data.get('activation_inhibition_synergy', 0) + 
 data.get('alpha_SMA_inhibition_synergy', 0) + 
 data.get('collagen_I_inhibition_synergy', 0)) / 3
 synergy_scores.append(synergy_score)
 
 # Inhibition efficiency comparison
 ax1 = axes[0, 0]
 x = np.arange(len(groups))
 width = 0.25
 ax1.bar(x - width, activation_inhibition, width, label='Cell Activation Inhibition', alpha=0.8)
 ax1.bar(x, alpha_SMA_inhibition, width, label='α-SMA Inhibition', alpha=0.8)
 ax1.bar(x + width, collagen_I_inhibition, width, label='Collagen I Inhibition', alpha=0.8)
 ax1.set_title('Inhibition Efficiency Comparison')
 ax1.set_xlabel('Experimental Groups')
 ax1.set_ylabel('Inhibition Efficiency')
 ax1.set_xticks(x)
 ax1.set_xticklabels(groups, rotation=45)
 ax1.legend()
 ax1.grid(True, alpha=0.3)
 
 # Synergy score
 ax2 = axes[0, 1]
 colors = ['lightblue' if 'miR' in g and 'dual' not in g else 'orange' for g in groups]
 ax2.bar(groups, synergy_scores, color=colors, alpha=0.8)
 ax2.set_title('Synergistic Enhancement Score')
 ax2.set_xlabel('Experimental Groups')
 ax2.set_ylabel('Synergy Score')
 ax2.set_xticklabels(groups, rotation=45)
 ax2.grid(True, alpha=0.3)
 ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5)
 
 # Heatmap: effects of different ratios
 ax3 = axes[1, 0]
 ratio_data = []
 ratio_labels = []
 if 'dual_miRNA_1_1' in self.synergy_results:
 ratio_data.append([
 self.synergy_results['dual_miRNA_1_1'].get('activation_inhibition', 0),
 self.synergy_results['dual_miRNA_1_1'].get('alpha_SMA_inhibition', 0),
 self.synergy_results['dual_miRNA_1_1'].get('collagen_I_inhibition', 0)
 ])
 ratio_labels.append('1:1')
 if 'dual_miRNA_1_2' in self.synergy_results:
 ratio_data.append([
 self.synergy_results['dual_miRNA_1_2'].get('activation_inhibition', 0),
 self.synergy_results['dual_miRNA_1_2'].get('alpha_SMA_inhibition', 0),
 self.synergy_results['dual_miRNA_1_2'].get('collagen_I_inhibition', 0)
 ])
 ratio_labels.append('1:2')
 if 'dual_miRNA_2_1' in self.synergy_results:
 ratio_data.append([
 self.synergy_results['dual_miRNA_2_1'].get('activation_inhibition', 0),
 self.synergy_results['dual_miRNA_2_1'].get('alpha_SMA_inhibition', 0),
 self.synergy_results['dual_miRNA_2_1'].get('collagen_I_inhibition', 0)
 ])
 ratio_labels.append('2:1')
 
 if ratio_data:
 sns.heatmap(ratio_data, 
 xticklabels=['Cell Activation', 'α-SMA', 'Collagen I'],
 yticklabels=ratio_labels,
 annot=True, fmt='.3f', cmap='RdYlBu_r', ax=ax3)
 ax3.set_title('Inhibition Effect Heatmap by miRNA Ratio')
 
 # Synergy coefficient comparison
 ax4 = axes[1, 1]
 dual_groups = [g for g in groups if 'dual' in g]
 if dual_groups:
 synergy_metrics = ['activation_inhibition_synergy', 'alpha_SMA_inhibition_synergy', 'collagen_I_inhibition_synergy']
 synergy_matrix = []
 for group in dual_groups:
 row = []
 for metric in synergy_metrics:
 row.append(self.synergy_results[group].get(metric, 0))
 synergy_matrix.append(row)
 
 sns.heatmap(synergy_matrix,
 xticklabels=['Cell Activation', 'α-SMA', 'Collagen I'],
 yticklabels=dual_groups,
 annot=True, fmt='.3f', cmap='RdBu_r', center=0, ax=ax4)
 ax4.set_title('Synergistic Enhancement Coefficient Heatmap')
 
 plt.tight_layout()
 plt.savefig('synergy_analysis.png', dpi=300, bbox_inches='tight')
 plt.show()
 
 def generate_report(self):
 """Generate analysis report"""
 print("Generating analysis report...")
 
 with open('analysis_report.md', 'w', encoding='utf-8') as f:
 f.write("# iGEM Liver Fibrosis Simulation Analysis Report\n\n")
 f.write(f"Analysis Time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
 
 f.write("## Experimental Groups Overview\n\n")
 f.write("| Group | Cell Count | Time Points | Data Status |\n")
 f.write("|-------|------------|-------------|-------------|\n")
 
 for group in self.experimental_groups:
 if group in self.results:
 data = self.results[group]
 cell_count = data['total_cells'][-1] if data['total_cells'] else 0
 time_points = len(data['time_points'])
 status = " Complete" if time_points > 0 else " Missing"
 else:
 cell_count = 0
 time_points = 0
 status = " Not Found"
 
 f.write(f"| {group} | {cell_count} | {time_points} | {status} |\n")
 
 if hasattr(self, 'synergy_results'):
 f.write("\n## Synergistic Enhancement Analysis\n\n")
 
 f.write("### Key Findings\n\n")
 
 # Find the best synergistic combination
 best_synergy = None
 best_score = -999
 for group in self.synergy_results:
 if 'dual' in group:
 score = (self.synergy_results[group].get('activation_inhibition_synergy', 0) + 
 self.synergy_results[group].get('alpha_SMA_inhibition_synergy', 0) + 
 self.synergy_results[group].get('collagen_I_inhibition_synergy', 0)) / 3
 if score > best_score:
 best_score = score
 best_synergy = group
 
 if best_synergy:
 f.write(f"- **Best Synergistic Combination**: {best_synergy} (Synergy Score: {best_score:.3f})\n")
 
 f.write("\n### Detailed Results\n\n")
 f.write("| Group | Cell Activation Inhibition | α-SMA Inhibition | Collagen I Inhibition | Synergy Effect |\n")
 f.write("|-------|----------------------------|------------------|----------------------|----------------|\n")
 
 for group in ['miR_455_3p', 'miR_148a_5p', 'dual_miRNA_1_1', 'dual_miRNA_1_2', 'dual_miRNA_2_1']:
 if group in self.synergy_results:
 data = self.synergy_results[group]
 synergy_score = (data.get('activation_inhibition_synergy', 0) + 
 data.get('alpha_SMA_inhibition_synergy', 0) + 
 data.get('collagen_I_inhibition_synergy', 0)) / 3
 
 f.write(f"| {group} | {data.get('activation_inhibition', 0):.3f} | "
 f"{data.get('alpha_SMA_inhibition', 0):.3f} | "
 f"{data.get('collagen_I_inhibition', 0):.3f} | "
 f"{synergy_score:.3f} |\n")
 
 f.write("\n## \n\n")
 f.write("1. checkexperimentgroupdata\n")
 f.write("2. analysistime\n")
 f.write("3. verifysynergistic effect\n")
 f.write("4. experimentdataverify\n")
 
 print("Analysis report generated: analysis_report.md")
 
 def run_full_analysis(self):
 """runanalysis"""
 print("Starting full analysis pipeline...")
 
 self.load_simulation_data()
 self.calculate_synergy_effects()
 self.plot_time_series()
 self.plot_synergy_analysis()
 self.generate_report()
 
 print("\nAnalysis completed!")
 print("Generated files:")
 print("- time_series_analysis.png")
 print("- synergy_analysis.png") 
 print("- analysis_report.md")

if __name__ == "__main__":
 # runanalysis
 analyzer = LiverFibrosisAnalyzer()
 analyzer.run_full_analysis()
