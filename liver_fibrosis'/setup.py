#!/usr/bin/env python3
"""
Setup script for Liver Fibrosis Simulation Package
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))

# Package metadata
setup(
    name="liver-fibrosis-simulation",
    version="1.0.0",
    author="iGEM Liver Fibrosis Team",
    author_email="contact@liverfibrosis.simulation",
    description="Advanced PhysiCell-based liver fibrosis simulation with miRNA therapeutic analysis",
    long_description="A comprehensive simulation platform for studying hepatic stellate cell activation and miRNA-based therapeutic interventions in liver fibrosis using PhysiCell framework",
    long_description_content_type="text/markdown",
    url="https://github.com/igem/liver-fibrosis-simulation",
    
    # Package configuration
    packages=find_packages(),
    package_data={
        'liver_fibrosis_sim': [
            'config/*.xml',
            'config/*.csv',
            'example_configs/*.xml',
            'custom_modules/*.cpp',
            'custom_modules/*.h',
            'Makefile',
            '*.cpp'
        ]
    },
    include_package_data=True,
    
    # Dependencies
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.3.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "scipy>=1.7.0",
        "plotly>=5.0.0",
        "scikit-learn>=1.0.0",
        "pillow>=8.0.0",
        "kaleido>=0.2.1"  # For plotly image export
    ],
    
    # Optional dependencies
    extras_require={
        'dev': [
            'pytest>=6.0.0',
            'pytest-cov>=2.0.0',
            'black>=21.0.0',
            'flake8>=3.8.0'
        ],
        'full': [
            'networkx>=2.5.0',
            'dash>=2.0.0',
            'jupyterlab>=3.0.0'
        ]
    },
    
    # Entry points for command-line tools
    entry_points={
        'console_scripts': [
            'liver-fibrosis-demo=liver_fibrosis_sim.demo:main',
            'liver-fibrosis-analyze=liver_fibrosis_sim.analysis:main',
            'liver-fibrosis-visualize=liver_fibrosis_sim.visualization:main',
        ],
    },
    
    # Classification
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
    
    # Python version requirement
    python_requires='>=3.7',
    
    # Keywords for searching
    keywords="liver fibrosis, miRNA, PhysiCell, simulation, systems biology, therapeutic analysis",
    
    # Project URLs
    project_urls={
        'Bug Reports': 'https://github.com/igem/liver-fibrosis-simulation/issues',
        'Source': 'https://github.com/igem/liver-fibrosis-simulation',
        'Documentation': 'https://liver-fibrosis-simulation.readthedocs.io/',
    }
)
