# Agent-Based Model Simulation Platform for miRNA-Based Liver Fibrosis Therapy

## Multi-Scale Computational Biology Modeling and Systems Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![PhysiCell](https://img.shields.io/badge/PhysiCell-1.10+-green.svg)](http://physicell.org/)

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [System Architecture](#system-architecture)
- [Experimental Design](#experimental-design)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Advanced Usage](#advanced-usage)
- [Mathematical Framework](#mathematical-framework)
- [Biological Models](#biological-models)
- [Visualization Suite](#visualization-suite)
- [Performance Optimization](#performance-optimization)
- [License](#license)

---

## Overview

This project presents a comprehensive **Agent-Based Model (ABM) simulation platform** for investigating miRNA-based therapeutic interventions in liver fibrosis. The platform integrates multi-scale computational biology modeling, from molecular interactions to tissue-level dynamics, providing a powerful tool for understanding hepatic stellate cell (HSC) activation mechanisms and evaluating novel therapeutic strategies.

### Research Focus

The platform specifically models the therapeutic potential of **miR-455-3p** and **miR-148a-5p** delivered via exosomes for treating liver fibrosis. Through systematic computational experiments, we investigate:

- Individual and synergistic effects of dual miRNA therapies
- Dose-response relationships and optimal treatment protocols
- Mechanistic insights into TGF-β signaling pathway modulation
- Spatial-temporal dynamics of fibrotic tissue remodeling

### Scientific Impact

This simulation platform bridges the gap between molecular mechanisms and clinical applications, offering:

- **Predictive modeling** for therapeutic efficacy assessment
- **Cost-effective screening** of treatment combinations
- **Mechanistic insights** into complex biological processes
- **Experimental design optimization** for wet lab validation

---

## Key Features

### 🧬 Multi-Scale Biological Modeling

- **Molecular Level**: TGF-β signaling networks, miRNA regulation, protein interactions
- **Cellular Level**: HSC activation, proliferation, apoptosis, migration dynamics
- **Tissue Level**: Extracellular matrix remodeling, mechanical properties, spatial organization

### 🔬 Comprehensive Experimental Design

- **9 Experimental Configurations**: Control, positive/negative controls, single and dual miRNA treatments
- **Systematic Parameter Variation**: Dose-response curves, ratio optimization, temporal dynamics
- **Statistical Rigor**: Multiple replicates, uncertainty quantification, power analysis

### 📊 Advanced Visualization Suite

- **2D/3D Spatial Visualization**: Cell distributions, molecular gradients, tissue architecture
- **Time-Lapse Animations**: Dynamic cellular behavior, treatment progression
- **Statistical Analysis**: Correlation matrices, dose-response curves, comparative analysis
- **Interactive Dashboards**: Real-time parameter exploration, sensitivity analysis

### ⚡ High-Performance Computing

- **Parallel Processing**: Multi-threading, distributed computing, GPU acceleration
- **Scalable Architecture**: Modular design, efficient memory management
- **Adaptive Algorithms**: Dynamic load balancing, intelligent task scheduling

### 🎯 Therapeutic Applications

- **Drug Discovery**: Virtual screening, combination therapy optimization
- **Clinical Translation**: Biomarker identification, patient stratification
- **Personalized Medicine**: Individual response prediction, treatment customization

---

## System Architecture

The platform employs a modular, hierarchical architecture designed for scalability and extensibility:

```
Liver_fibrosis/liver_fibrosis/
├── mathematical_core/          # Core mathematical algorithms
│   ├── multiscale_solvers.py   # ADI, RKF, operator splitting
│   ├── stochastic_processes.py # MCMC, Gillespie, SDE solvers
│   └── optimization_algorithms.py # Bayesian optimization, NSGA-II
├── biological_models/          # Biological mechanism models
│   ├── cellular_dynamics.py    # Cell cycle, metabolism, migration
│   ├── molecular_networks.py   # TGF-β pathway, miRNA regulation
│   └── membrane_transport.py   # Endocytosis, exocytosis, trafficking
├── multiscale_integration/     # Multi-scale coupling systems
│   ├── scale_coupling.py       # Cross-scale information transfer
│   ├── reaction_diffusion.py   # Molecular transport and reactions
│   └── tissue_mechanics.py     # Continuum mechanics, stress fields
├── advanced_algorithms/        # High-performance computing
│   ├── parallel_computing.py   # Load balancing, distributed memory
│   ├── data_analysis.py        # t-SNE, PCA, statistical analysis
│   └── parameter_calibration.py # Bayesian inference, MCMC
├── visualization/              # Comprehensive visualization suite
│   ├── professional_visualization.py # Publication-quality figures
│   ├── 3D_visualization.py     # Interactive 3D rendering
│   └── advanced_visualization.py # Statistical and network plots
└── config/                     # Configuration and parameters
    ├── experiment_config.py    # Experimental design parameters
    └── PhysiCell_settings.xml  # PhysiCell simulation settings
```

### Core Components

#### 1. Mathematical Core
- **Numerical Solvers**: Advanced PDE/ODE solvers with adaptive time stepping
- **Stochastic Algorithms**: Monte Carlo methods, random process simulation
- **Optimization Framework**: Multi-objective optimization, parameter estimation

#### 2. Biological Models
- **Cellular Dynamics**: 33 state variables, 12 behavior modules
- **Molecular Networks**: TGF-β signaling, miRNA regulation, synergy effects
- **Transport Systems**: Membrane trafficking, vesicle dynamics, autophagy

#### 3. Multi-Scale Integration
- **Scale Coupling**: Molecular ↔ Cellular ↔ Tissue information transfer
- **Reaction-Diffusion**: 12 molecular species with anisotropic diffusion
- **Tissue Mechanics**: Elastic-viscous behavior, mechanical signaling

#### 4. Advanced Algorithms
- **Parallel Computing**: OpenMP/MPI hybrid parallelization
- **Data Analysis**: Dimensionality reduction, clustering, correlation analysis
- **Parameter Calibration**: Bayesian inference, uncertainty quantification

---

## Experimental Design

### Nine Experimental Configurations

The platform implements a systematic experimental design with nine distinct configurations:

| Configuration | TGF-β1 (ng/mL) | Exosomes (μg/mL) | miR-455-3p | miR-148a-5p | Expected Activation |
|---------------|-----------------|------------------|------------|-------------|-------------------|
| **Control** | 2.0 | 0.0 | 0.0 | 0.0 | 8% |
| **Positive Control** | 15.0 | 0.0 | 0.0 | 0.0 | 75% |
| **Natural Exosomes** | 15.0 | 25.0 | 0.05 | 0.05 | 65% |
| **NC Mimic** | 15.0 | 50.0 | 0.02* | 0.02* | 70% |
| **miR-455-3p Only** | 15.0 | 50.0 | 0.4 | 0.0 | 45% |
| **miR-148a-5p Only** | 15.0 | 50.0 | 0.0 | 0.4 | 50% |
| **Dual miRNA (1:1)** | 15.0 | 75.0 | 0.25 | 0.25 | 28% |
| **Dual miRNA (2:1)** | 15.0 | 75.0 | 0.35 | 0.20 | 30% |
| **Dual miRNA (1:2)** | 15.0 | 75.0 | 0.20 | 0.35 | 25% |

*Scrambled sequence controls

### Mathematical Modeling Framework

#### Dose-Response Relationships
miRNA therapeutic effects follow Hill equation kinetics:

$$E_{miRNA} = E_{max} \cdot \frac{D^{n_H}}{EC_{50}^{n_H} + D^{n_H}}$$

where $E_{max}$ is maximum effect, $EC_{50}$ is half-effect concentration, and $n_H$ is the Hill coefficient.

#### Synergy Quantification
Dual miRNA synergy is evaluated using Bliss independence:

$$E_{synergy} = E_{obs} - E_{exp}$$

where $E_{exp} = E_A + E_B - E_A \cdot E_B$ and synergy index $CI = E_{obs}/E_{exp}$.

#### Temporal Dynamics
All experiments use standardized temporal parameters:
- **Total simulation time**: 2880 minutes (48 hours)
- **Sampling interval**: 60 minutes
- **Key time points**: 24h (early), 48h (steady-state), 72h (long-term)

---

## Installation

### Prerequisites

- **Python 3.8+** with scientific computing libraries
- **PhysiCell 1.10+** for agent-based modeling
- **C++ compiler** (GCC 7+ or MSVC 2019+)
- **Git** for version control

### System Requirements

- **Memory**: 8GB RAM minimum, 16GB recommended
- **Storage**: 5GB free space for installation and results
- **CPU**: Multi-core processor recommended for parallel computing
- **OS**: Windows 10+, macOS 10.15+, or Linux (Ubuntu 18.04+)

### Installation Steps

1. **Clone the Repository**
```bash
   git clone https://github.com/your-org/liver-fibrosis-simulation.git
   cd liver-fibrosis-simulation/PhysiCell-master/liver_fibrosis/liver_fibrosis
```

2. **Install Python Dependencies**
```bash
    pip install -r requirements.txt
```

3. **Compile PhysiCell**
```bash
   make clean
   make -j4
```


### Docker Installation (Recommended)

```bash
docker pull physicell/liver-fibrosis:latest
docker run -it -v $(pwd):/workspace physicell/liver-fibrosis:latest
```

---

## Quick Start

### Basic Simulation

Run a complete simulation with default parameters:

```python
from liver_fibrosis_sim import LiverFibrosisSimulation

# Initialize simulation
sim = LiverFibrosisSimulation()

# Run single experiment
results = sim.run_experiment('dual_miRNA_1_1', duration_hours=48)

# Generate visualizations
sim.create_visualizations(results, output_dir='results/')
```

### Multi-Experiment Analysis

Execute all nine experimental configurations:

```python
from experiment_config import experiment_manager
from professional_visualization import ProfessionalBioVisualization

# Run all experiments
manager = experiment_manager
results = manager.run_all_experiments()

# Create comprehensive analysis
viz = ProfessionalBioVisualization()
viz.run_professional_visualization_suite()
```

### 3D Visualization

Generate interactive 3D visualizations:

```python
from advanced_visualization import Advanced3DVisualization

# Create 3D analyzer
viz_3d = Advanced3DVisualization()

# Generate molecular gradients
viz_3d.create_3D_molecular_gradients('dual_miRNA_1_1')

# Create time-lapse animation
viz_3d.create_3D_time_lapse_animation('dual_miRNA_1_1')
```

---

## Advanced Usage

### Custom Experimental Design

Define custom experimental parameters:

```python
from experiment_config import ExperimentConfig

# Create custom configuration
custom_config = ExperimentConfig(
    name='custom_treatment',
    miR_455_3p_dose=0.6,
    miR_148a_5p_dose=0.3,
    TGF_beta1_concentration=12.0,
    exosome_concentration=60.0
)

# Run custom experiment
results = sim.run_experiment(custom_config)
```

### Parameter Sensitivity Analysis

Perform global sensitivity analysis:

```python
from advanced_algorithms.parameter_calibration import SensitivityAnalyzer

# Initialize analyzer
analyzer = SensitivityAnalyzer()

# Define parameter bounds
bounds = [
    (0.0, 1.0),  # miR-455-3p dose
    (0.0, 1.0),  # miR-148a-5p dose
    (5.0, 20.0), # TGF-β1 concentration
]

# Calculate Sobol indices
sobol_results = analyzer.sobol_indices(
    model_function=sim.run_model,
    parameter_bounds=bounds,
    n_samples=1000
)
```

### Parallel Computing

Enable high-performance parallel execution:

```python
from advanced_algorithms.parallel_computing import ParallelComputingFramework

# Initialize parallel framework
parallel = ParallelComputingFramework(
    strategy='hybrid',
    num_workers=8
)

# Submit parallel tasks
tasks = [create_simulation_task(config) for config in all_configs]
results = parallel.execute_parallel(tasks)
```

### Bayesian Parameter Calibration

Calibrate model parameters using experimental data:

```python
from advanced_algorithms.parameter_calibration import BayesianParameterInference
from advanced_algorithms.parameter_calibration import MarkovChainMonteCarlo

# Define priors
priors = {
    'activation_rate': PriorDistribution(PriorType.LOGNORMAL, {'mean': 0.0, 'sigma': 0.5}),
    'degradation_rate': PriorDistribution(PriorType.GAMMA, {'shape': 2.0, 'scale': 0.1})
}

# Initialize Bayesian inference
inference = BayesianParameterInference(
    model_function=sim.run_model,
    parameter_names=['activation_rate', 'degradation_rate'],
    priors=priors
)

# Run MCMC sampling
mcmc = MarkovChainMonteCarlo(inference)
calibration_results = mcmc.metropolis_hastings(
    observed_data=experimental_data,
    n_samples=10000,
    burn_in=2000
)
```

---

## Mathematical Framework

### Multi-Scale Coupling Theory

The platform implements rigorous multi-scale coupling based on scale separation principles:

#### Molecular Scale
Reaction-diffusion equations govern molecular transport:

$$\frac{\partial C_i}{\partial t} = \nabla \cdot (D_i \nabla C_i) + R_i(\mathbf{C}) + S_i(\mathbf{x}, t)$$

where $C_i$ is concentration of species $i$, $D_i$ is the diffusion tensor, $R_i$ represents reactions, and $S_i$ are source terms.

#### Cellular Scale
Agent-based dynamics with stochastic state transitions:

$$\frac{d\mathbf{X}_j}{dt} = \mathbf{v}_j(\mathbf{X}_j, \mathbf{C}, t) + \boldsymbol{\eta}_j(t)$$

where $\mathbf{X}_j$ is the state of cell $j$, $\mathbf{v}_j$ is the deterministic velocity, and $\boldsymbol{\eta}_j$ represents stochastic forces.

#### Tissue Scale
Continuum mechanics with elastic-viscous behavior:

$$\rho \frac{\partial^2 \mathbf{u}}{\partial t^2} = \nabla \cdot \boldsymbol{\sigma} + \mathbf{f}$$

where $\mathbf{u}$ is displacement, $\boldsymbol{\sigma}$ is the stress tensor, and $\mathbf{f}$ represents body forces.

### Numerical Methods

#### Advanced Solvers
- **ADI (Alternating Direction Implicit)**: Efficient 3D diffusion solving
- **Runge-Kutta-Fehlberg**: Adaptive time stepping for stiff ODEs
- **Operator Splitting**: Decoupling of reaction and diffusion processes
- **Multiple Time Stepping**: Multi-scale temporal integration

#### Stochastic Algorithms
- **Gillespie Algorithm**: Exact stochastic chemical kinetics
- **Metropolis-Hastings MCMC**: Bayesian parameter inference
- **Ornstein-Uhlenbeck Processes**: Correlated noise modeling

#### Optimization Techniques
- **Bayesian Optimization**: Efficient parameter space exploration
- **NSGA-II**: Multi-objective evolutionary optimization
- **Sobol Sensitivity Analysis**: Global parameter sensitivity

---

## Biological Models

### Cellular Behavior Modules

The platform models 12 distinct cellular behavior modules:

#### 1. Cell Cycle Dynamics
- **G0/G1/S/G2/M Phase Transitions**: Checkpoint-controlled progression
- **Growth Factor Dependence**: External signal integration
- **DNA Damage Response**: p53-mediated cell cycle arrest

#### 2. Metabolic Networks
- **ATP Dynamics**: Glycolysis and oxidative phosphorylation
- **Metabolic Switching**: Adaptation to environmental conditions
- **Energy-Dependent Processes**: Transport, synthesis, maintenance

#### 3. Stress Response Systems
- **Oxidative Stress**: ROS generation and antioxidant responses
- **ER Stress**: Unfolded protein response and adaptation
- **Mechanical Stress**: Mechanosensitive channel activation

#### 4. Migration and Motility
- **Chemotaxis**: Gradient-sensing and directional migration
- **Random Motility**: Brownian motion and persistence
- **Contact Inhibition**: Density-dependent migration suppression

### Molecular Interaction Networks

#### TGF-β Signaling Pathway
Comprehensive model including:
- **Smad-Dependent Pathway**: Canonical TGF-β signaling
- **Non-Smad Pathways**: TAK1, p38, JNK activation
- **Feedback Regulation**: Smad7-mediated inhibition
- **miRNA Modulation**: Post-transcriptional regulation

#### miRNA Regulatory Networks
- **Processing Pathway**: pri-miRNA → pre-miRNA → mature miRNA
- **Target Recognition**: Seed sequence matching and binding kinetics
- **Synergistic Effects**: Cooperative target regulation
- **Degradation Dynamics**: miRNA stability and turnover

### Membrane Transport Systems

#### Endocytic Pathways
- **Clathrin-Mediated Endocytosis**: Receptor-mediated uptake
- **Caveolin-Mediated Endocytosis**: Lipid raft-dependent internalization
- **Macropinocytosis**: Bulk fluid uptake mechanisms

#### Intracellular Trafficking
- **Endosome Maturation**: Early to late endosome progression
- **Lysosomal Fusion**: Cargo degradation pathways
- **Recycling Pathways**: Receptor and membrane recycling

---

## Visualization Suite

### Professional 2D Visualizations

#### Morphological Analysis
- **Cell Size Distributions**: Population heterogeneity analysis
- **Elongation Metrics**: Shape factor quantification
- **Stress Fiber Density**: Cytoskeletal organization assessment
- **Proliferation/Apoptosis Rates**: Cell fate analysis

#### Molecular Heatmaps
- **Concentration Fields**: Spatial distribution of key molecules
- **Gradient Analysis**: Directional transport visualization
- **Treatment Efficacy Maps**: Therapeutic effect quantification

#### Statistical Analysis
- **Dose-Response Curves**: Hill equation fitting and EC50 determination
- **Correlation Matrices**: Multi-parameter relationship analysis
- **Time Series Analysis**: Temporal dynamics and trend detection

### Interactive 3D Visualizations

#### Spatial Cell Distributions
- **3D Cell Positioning**: Real-time cell tracking and visualization
- **Activation State Mapping**: Color-coded cellular phenotypes
- **Migration Trajectories**: Path analysis and velocity fields

#### Molecular Gradient Fields
- **Isosurface Rendering**: 3D concentration visualization
- **Vector Field Display**: Flux and gradient directions
- **Time-Lapse Animation**: Dynamic molecular transport

#### Tissue Architecture
- **Extracellular Matrix**: Fiber orientation and density
- **Mechanical Properties**: Stress and strain field visualization
- **Vascular Networks**: Blood vessel architecture modeling

### Advanced Analytics

#### Dimensionality Reduction
- **t-SNE Analysis**: High-dimensional data visualization
- **PCA Decomposition**: Principal component identification
- **Clustering Analysis**: Cell population segmentation

#### Network Analysis
- **Pathway Visualization**: Molecular interaction networks
- **Graph Metrics**: Centrality and connectivity analysis
- **Dynamic Networks**: Time-evolving interaction patterns

---

## Performance Optimization

### Parallel Computing Architecture

#### Multi-Threading Support
- **OpenMP Integration**: Shared memory parallelization
- **Thread Pool Management**: Efficient task distribution
- **Load Balancing**: Dynamic work allocation

#### Distributed Computing
- **MPI Implementation**: Message passing for cluster computing
- **Hybrid Parallelization**: Combined OpenMP/MPI approach
- **Scalable Architecture**: Linear scaling to 100+ cores

#### GPU Acceleration
- **CUDA Integration**: GPU-accelerated computations
- **Memory Optimization**: Efficient GPU memory management
- **Kernel Optimization**: High-performance computing kernels

### Memory Management

#### Efficient Data Structures
- **Sparse Matrices**: Memory-efficient linear algebra
- **Compressed Storage**: Reduced memory footprint
- **Cache Optimization**: CPU cache-friendly algorithms

#### Adaptive Mesh Refinement
- **Dynamic Grids**: Automatic spatial resolution adjustment
- **Error Estimation**: Local refinement criteria
- **Memory Pooling**: Efficient allocation strategies

### Algorithmic Optimizations

#### Adaptive Time Stepping
- **Error Control**: Automatic step size adjustment
- **Stability Analysis**: CFL condition monitoring
- **Multi-Rate Integration**: Different time scales handling

#### Preconditioning
- **Matrix Preconditioning**: Accelerated linear solver convergence
- **Multigrid Methods**: Hierarchical solution strategies
- **Domain Decomposition**: Parallel solver optimization

---

## Contributing

We welcome contributions from the computational biology, systems biology, and biomedical engineering communities.

### Development Guidelines

#### Code Standards
- **Python Style**: Follow PEP 8 guidelines
- **Documentation**: Comprehensive docstrings and comments
- **Testing**: Unit tests with >90% coverage
- **Version Control**: Git workflow with feature branches

#### Contribution Process
1. **Fork Repository**: Create personal fork on GitHub
2. **Feature Branch**: Develop in dedicated feature branch
3. **Testing**: Ensure all tests pass and add new tests
4. **Documentation**: Update documentation and examples
5. **Pull Request**: Submit PR with detailed description

#### Areas for Contribution
- **New Biological Models**: Additional cellular behaviors and pathways
- **Numerical Methods**: Advanced solvers and optimization algorithms
- **Visualization Tools**: Novel analysis and presentation methods
- **Performance Optimization**: Parallelization and GPU acceleration
- **Experimental Validation**: Comparison with wet lab data


---



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Third-Party Libraries

- **PhysiCell**: BSD 3-Clause License
- **NumPy/SciPy**: BSD License
- **Matplotlib**: PSF License
- **scikit-learn**: BSD License
- **Plotly**: MIT License

---

