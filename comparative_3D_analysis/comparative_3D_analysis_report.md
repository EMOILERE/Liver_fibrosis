# Comparative 3D Analysis Report

Generated: 2025-10-07 13:47:24

## Experiment Overview

### Control
- **Configuration**: control
- **miR-455-3p dose**: 0.0
- **miR-148a-5p dose**: 0.0
- **Exosome delivery**: No

### Dual miRNA (1:1)
- **Configuration**: dual_miRNA_1_1
- **miR-455-3p dose**: 1.0
- **miR-148a-5p dose**: 1.0
- **Exosome delivery**: Yes

### Dual miRNA (2:1)
- **Configuration**: dual_miRNA_2_1
- **miR-455-3p dose**: 2.0
- **miR-148a-5p dose**: 1.0
- **Exosome delivery**: Yes

### Dual miRNA (1:2)
- **Configuration**: dual_miRNA_1_2
- **miR-455-3p dose**: 1.0
- **miR-148a-5p dose**: 2.0
- **Exosome delivery**: Yes

### miR-455-3p Only
- **Configuration**: miR_455_only
- **miR-455-3p dose**: 1.0
- **miR-148a-5p dose**: 0.0
- **Exosome delivery**: Yes

### miR-148a-5p Only
- **Configuration**: miR_148_only
- **miR-455-3p dose**: 0.0
- **miR-148a-5p dose**: 1.0
- **Exosome delivery**: Yes

### Natural Exosomes
- **Configuration**: natural_exosomes
- **miR-455-3p dose**: 0.0
- **miR-148a-5p dose**: 0.0
- **Exosome delivery**: Yes

### NC Mimic
- **Configuration**: negative_control
- **miR-455-3p dose**: 0.0
- **miR-148a-5p dose**: 0.0
- **Exosome delivery**: Yes

### Positive Control
- **Configuration**: positive_model
- **miR-455-3p dose**: 0.0
- **miR-148a-5p dose**: 0.0
- **Exosome delivery**: No

## 3D Analysis Components

### Generated for Each Experiment
1. **Interactive 3D Cell Visualization**
   - Real-time 3D scatter plots of cells colored by biological state
   - Interactive rotation, zoom, and hover information
   - Culture well boundary visualization

2. **3D Molecular Gradient Visualization**
   - Six key molecules: Oxygen, TGF-β1, α-SMA, Glucose, Lactate, Collagen-I
   - 3D spatial distribution with optimized colorbars
   - Realistic gradient patterns based on biology

3. **3D Spatial Pattern Analysis**
   - Vertical (Z-axis) cell distribution
   - Radial distribution from culture center
   - DBSCAN clustering analysis
   - Activation vs position correlations
   - miRNA co-expression patterns
   - Metabolic activity gradients

4. **Dynamic Time-lapse Animation**
   - Interactive 48-hour time course
   - Play/pause controls and time slider
   - Real-time statistics overlay
   - Cell state evolution tracking

## Key 3D Insights

- **Spatial Heterogeneity**: 3D culture captures realistic tissue-level complexity
- **Gradient Formation**: Physiologically relevant molecular gradients established
- **Cell-Cell Interactions**: Enhanced 3D cell communication patterns
- **Treatment Penetration**: miRNA delivery shows realistic 3D distribution
- **Temporal Dynamics**: Time-course reveals treatment kinetics

## File Organization

```
3D_visualization_outputs/
├── experiment_control/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_dual_miRNA_1_1/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_dual_miRNA_2_1/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_dual_miRNA_1_2/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_miR_455_only/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_miR_148_only/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_natural_exosomes/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_negative_control/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
├── experiment_positive_model/
│   ├── interactive_html/    # 3D visualizations
│   ├── images/              # Static analyses
│   ├── animations/          # Time-lapse files
│   └── reports/             # Analysis reports
└── comparative_3D_analysis/    # Cross-experiment comparisons
```

---

**Advanced 3D Visualization System**
Enhanced platform for 3D tissue culture simulation analysis
