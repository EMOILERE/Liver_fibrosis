#!/usr/bin/env python3
"""
Quick Configuration Script
For rapid modification of key parameters in simulation_parameters.xml
"""

import xml.etree.ElementTree as ET
import argparse
import sys
import os

def load_config(filename):
    """Load configuration file"""
    try:
        tree = ET.parse(filename)
        return tree, tree.getroot()
    except Exception as e:
        print(f"Error loading config file: {e}")
        return None, None

def save_config(tree, filename):
    """Save configuration file"""
    try:
        tree.write(filename, encoding='utf-8', xml_declaration=True)
        print(f"Configuration saved to: {filename}")
        return True
    except Exception as e:
        print(f"Error saving config file: {e}")
        return False

def set_experiment_group(root, group_name):
    """Set experiment group - supports multiple configuration formats"""
    updated = False
    
    # 1. Try simulation_parameters.xml format
    experiment_node = root.find('.//experimental_conditions/current_experiment')
    if experiment_node is not None:
        experiment_node.text = group_name
        updated = True
    
    # 2. Try PhysiCell_settings.xml format - string elements in user_parameters
    user_params = root.find('.//user_parameters')
    if user_params is not None:
        # Look for existing experimental_group parameter
        for param in user_params.findall('string'):
            if param.get('name') == 'experimental_group':
                param.text = group_name
                updated = True
                break
        else:
            # If not found, create new experimental_group parameter
            new_param = ET.SubElement(user_params, 'string')
            new_param.set('name', 'experimental_group')
            new_param.set('units', 'dimensionless')
            new_param.text = group_name
            updated = True
    
    # 3. Try PhysiCell_settings.xml format - standalone experimental_group element
    exp_group_node = root.find('.//experimental_group')
    if exp_group_node is not None:
        exp_group_node.text = group_name
        updated = True
    
    if updated:
        print(f"Set experiment group to: {group_name}")
        return True
    else:
        print("Error: Could not find experimental group configuration in any format")
        return False

def set_synergy_factor(root, factor):
    """Set synergistic enhancement factor"""
    synergy_node = root.find('.//synergy_mechanism/max_synergy_factor')
    if synergy_node is not None:
        synergy_node.text = str(factor)
        print(f"Set max synergy factor to: {factor}")
        return True
    else:
        print("Error: Could not find max_synergy_factor node")
        return False

def set_tgf_concentration(root, concentration):
    """Set TGF-β1 concentration"""
    tgf_node = root.find('.//microenvironment/TGF_beta1/concentration')
    if tgf_node is not None:
        tgf_node.text = str(concentration)
        print(f"Set TGF-β1 concentration to: {concentration} ng/mL")
        return True
    else:
        print("Error: Could not find TGF_beta1 concentration node")
        return False

def set_cell_count(root, count):
    """Set initial cell count"""
    cell_node = root.find('.//cell_population/initial_cell_count')
    if cell_node is not None:
        cell_node.text = str(count)
        print(f"Set initial cell count to: {count}")
        return True
    else:
        print("Error: Could not find initial_cell_count node")
        return False

def set_simulation_time(root, time_hours):
    """Set simulation time (in hours)"""
    time_node = root.find('.//basic_settings/simulation_time')
    if time_node is not None:
        time_minutes = time_hours * 60
        time_node.text = str(time_minutes)
        print(f"Set simulation time to: {time_hours} hours ({time_minutes} minutes)")
        return True
    else:
        print("Error: Could not find simulation_time node")
        return False

def set_activation_threshold(root, threshold):
    """Set cell activation threshold"""
    threshold_node = root.find('.//cell_activation/TGF_beta1_response/activation_threshold')
    if threshold_node is not None:
        threshold_node.text = str(threshold)
        print(f"Set activation threshold to: {threshold} ng/mL")
        return True
    else:
        print("Error: Could not find activation_threshold node")
        return False

def list_available_experiments(root):
    """List all available experiment groups"""
    print("\nAvailable experiment groups:")
    groups_node = root.find('.//experimental_conditions/experiment_groups')
    if groups_node is not None:
        for i, group in enumerate(groups_node, 1):
            description = group.get('description', 'No description')
            print(f"  {i:2d}. {group.tag:<20} - {description}")
    print()

def set_simulation_duration(root, duration_minutes):
    """Set extended simulation duration in minutes"""
    updated = False
    
    # Try PhysiCell_settings.xml format first
    user_params = root.find('.//user_parameters')
    if user_params is not None:
        # Look for existing parameter
        for param in user_params.findall('double'):
            if param.get('name') == 'simulation_duration_hours':
                param.text = str(duration_minutes / 60.0)
                updated = True
                break
        else:
            # Create new parameter
            new_param = ET.SubElement(user_params, 'double')
            new_param.set('name', 'simulation_duration_hours')
            new_param.text = str(duration_minutes / 60.0)
            updated = True
    
    # Also update max_time in overall_parameters
    overall_params = root.find('.//overall')
    if overall_params is not None:
        max_time = overall_params.find('max_time')
        if max_time is not None:
            max_time.text = str(duration_minutes)
            updated = True
    
    if updated:
        print(f"Set simulation duration to {duration_minutes} minutes ({duration_minutes/60:.1f} hours)")
    
    return updated

def set_3d_simulation(root, enable_3d):
    """Enable or disable 3D simulation"""
    updated = False
    
    user_params = root.find('.//user_parameters')
    if user_params is not None:
        # Look for existing parameter
        for param in user_params.findall('bool'):
            if param.get('name') == 'enable_3D_simulation':
                param.text = 'true' if enable_3d else 'false'
                updated = True
                break
        else:
            # Create new parameter
            new_param = ET.SubElement(user_params, 'bool')
            new_param.set('name', 'enable_3D_simulation')
            new_param.text = 'true' if enable_3d else 'false'
            updated = True
    
    if updated:
        mode = "3D" if enable_3d else "2D"
        print(f"Set simulation mode to {mode}")
    
    return updated

def set_noise_reduction(root, reduce_noise):
    """Enable or disable noise reduction for ideal data"""
    updated = False
    
    user_params = root.find('.//user_parameters')
    if user_params is not None:
        # Look for existing parameter
        for param in user_params.findall('bool'):
            if param.get('name') == 'reduce_noise':
                param.text = 'true' if reduce_noise else 'false'
                updated = True
                break
        else:
            # Create new parameter
            new_param = ET.SubElement(user_params, 'bool')
            new_param.set('name', 'reduce_noise')
            new_param.text = 'true' if reduce_noise else 'false'
            updated = True
    
    if updated:
        status = "enabled" if reduce_noise else "disabled"
        print(f"Noise reduction {status}")
    
    return updated

def reset_to_defaults(root):
    """Reset all parameters to default values"""
    defaults = {
        'TGF_beta1_concentration': 10.0,
        'exosome_concentration': 50.0,
        'experimental_group': 'positive_model',
        'miR_455_3p_ratio': 0.5,
        'miR_148a_5p_ratio': 0.5,
        'activation_threshold': 3.0,
        'miRNA_degradation_rate': 0.005,
        'max_synergy_factor': 3.5,
        'initial_LX2_count': 120,
        'simulation_duration_hours': 120.0,
        'enable_3D_simulation': False,
        'reduce_noise': True
    }
    
    updated = False
    user_params = root.find('.//user_parameters')
    if user_params is not None:
        for param_name, default_value in defaults.items():
            param_type = 'double' if isinstance(default_value, float) else 'int' if isinstance(default_value, int) else 'bool' if isinstance(default_value, bool) else 'string'
            
            # Find existing parameter
            found = False
            for param in user_params.findall(param_type):
                if param.get('name') == param_name:
                    if param_type == 'bool':
                        param.text = 'true' if default_value else 'false'
                    else:
                        param.text = str(default_value)
                    found = True
                    updated = True
                    break
            
            # Create if not found
            if not found:
                new_param = ET.SubElement(user_params, param_type)
                new_param.set('name', param_name)
                if param_type == 'bool':
                    new_param.text = 'true' if default_value else 'false'
                else:
                    new_param.text = str(default_value)
                updated = True
    
    if updated:
        print("Reset all parameters to enhanced defaults")
    
    return updated

def create_preset_configs():
    """Create preset configuration files"""
    presets = {
        'fast_test': {
            'simulation_time': 720,  # 12 hours
            'initial_cell_count': 50,
            'save_interval': 30,
            'description': 'Fast test configuration'
        },
        'high_synergy': {
            'max_synergy_factor': 3.5,
            'synergy_threshold': 0.05,
            'description': 'High synergy effect configuration'
        },
        'low_synergy': {
            'max_synergy_factor': 1.5,
            'synergy_threshold': 0.2,
            'description': 'Low synergy effect configuration'
        },
        'high_dose': {
            'TGF_beta1_concentration': 20.0,
            'exosome_concentration': 100.0,
            'description': 'High dose configuration'
        }
    }
    
    # Try multiple possible paths
    possible_paths = [
        'simulation_parameters.xml',
        'config/simulation_parameters.xml',
        '../config/simulation_parameters.xml',
        '../../config/simulation_parameters.xml',
        '../../../config/simulation_parameters.xml'
    ]
    
    base_file = None
    for path in possible_paths:
        if os.path.exists(path):
            base_file = path
            break
    
    if base_file is None:
        print(f"Error: Base config file simulation_parameters.xml not found in any of: {possible_paths}")
        return
    
    tree, root = load_config(base_file)
    if tree is None:
        return
    
    for preset_name, params in presets.items():
        # Copy base configuration
        preset_tree = ET.ElementTree(ET.fromstring(ET.tostring(root)))
        preset_root = preset_tree.getroot()
        
        # Apply preset parameters
        for param, value in params.items():
            if param == 'description':
                continue
            elif param == 'simulation_time':
                set_simulation_time(preset_root, value / 60)  # Convert to hours
            elif param == 'initial_cell_count':
                set_cell_count(preset_root, value)
            elif param == 'max_synergy_factor':
                set_synergy_factor(preset_root, value)
            elif param == 'synergy_threshold':
                node = preset_root.find('.//synergy_mechanism/synergy_threshold')
                if node is not None:
                    node.text = str(value)
            elif param == 'TGF_beta1_concentration':
                set_tgf_concentration(preset_root, value)
            elif param == 'exosome_concentration':
                node = preset_root.find('.//microenvironment/exosomes/concentration')
                if node is not None:
                    node.text = str(value)
        
        # Save preset configuration
        preset_filename = f'simulation_parameters_{preset_name}.xml'
        save_config(preset_tree, preset_filename)
        print(f"Created preset: {preset_filename} - {params.get('description', '')}")

def main():
    parser = argparse.ArgumentParser(description='Quick configuration for liver fibrosis simulation parameters')
    parser.add_argument('--config', '-c', default='config/PhysiCell_settings.xml', 
                       help='Configuration file path (default: config/PhysiCell_settings.xml)')
    parser.add_argument('--output', '-o', help='Output file path (default: overwrite input file)')
    
    # Parameter setting options
    parser.add_argument('--experiment', '-e', help='Set experiment group')
    parser.add_argument('--synergy', '-s', type=float, help='Set synergistic enhancement factor')
    parser.add_argument('--tgf', '-t', type=float, help='Set TGF-β1 concentration (ng/mL)')
    parser.add_argument('--cells', type=int, help='Set initial cell count')
    parser.add_argument('--time', type=float, help='Set simulation time (hours)')
    parser.add_argument('--threshold', type=float, help='Set activation threshold (ng/mL)')
    
    # Enhanced parameters
    parser.add_argument('--duration', type=float, help='Set simulation duration in minutes (default: 2880 = 48h)')
    parser.add_argument('--3d', action='store_true', help='Enable 3D simulation (default: 2D)')
    parser.add_argument('--2d', action='store_true', help='Force 2D simulation')
    parser.add_argument('--reduce-noise', action='store_true', help='Enable noise reduction for ideal data')
    parser.add_argument('--reset-defaults', action='store_true', help='Reset all parameters to defaults')
    
    # Other options
    parser.add_argument('--list', '-l', action='store_true', help='List all available experiment groups')
    parser.add_argument('--create-presets', action='store_true', help='Create preset configuration files')
    
    args = parser.parse_args()
    
    # Create preset configurations
    if args.create_presets:
        create_preset_configs()
        return
    
    # Check if configuration file exists, use flexible path search
    config_file = args.config
    if not os.path.exists(config_file):
        # If specified path doesn't exist, try common locations
        possible_paths = [
            args.config,
            f'config/{args.config}',
            f'../{args.config}',
            'config/PhysiCell_settings.xml',
            '../config/PhysiCell_settings.xml',
            'PhysiCell_settings.xml',
            'config/simulation_parameters.xml',
            '../config/simulation_parameters.xml'
        ]
        
        config_file = None
        for path in possible_paths:
            if os.path.exists(path):
                config_file = path
                break
        
        if config_file is None:
            print(f"Error: Config file not found in any of: {possible_paths}")
            return
    
    # Load configuration
    tree, root = load_config(config_file)
    if tree is None:
        return
    
    # List experiment groups
    if args.list:
        list_available_experiments(root)
        return
    
    # Apply parameter modifications
    modified = False
    
    if args.experiment:
        if set_experiment_group(root, args.experiment):
            modified = True
    
    if args.synergy:
        if set_synergy_factor(root, args.synergy):
            modified = True
    
    if args.tgf:
        if set_tgf_concentration(root, args.tgf):
            modified = True
    
    if args.cells:
        if set_cell_count(root, args.cells):
            modified = True
    
    if args.time:
        if set_simulation_time(root, args.time):
            modified = True
    
    if args.threshold:
        if set_activation_threshold(root, args.threshold):
            modified = True
    
    # Handle new enhanced parameters
    if args.duration:
        if set_simulation_duration(root, args.duration):
            modified = True
    
    if getattr(args, '3d', False):
        if set_3d_simulation(root, True):
            modified = True
    elif getattr(args, '2d', False):
        if set_3d_simulation(root, False):
            modified = True
    
    if args.reduce_noise:
        if set_noise_reduction(root, True):
            modified = True
    
    if args.reset_defaults:
        if reset_to_defaults(root):
            modified = True
    
    # Save configuration
    if modified:
        output_file = args.output if args.output else args.config
        save_config(tree, output_file)
    else:
        print("No modifications made. Use --help for usage information.")

if __name__ == '__main__':
    main()
