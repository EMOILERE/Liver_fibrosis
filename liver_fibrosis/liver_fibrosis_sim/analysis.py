"""
Analysis Module for Liver Fibrosis Simulation
Statistical analysis and data processing tools
"""

from pathlib import Path
import sys
from typing import Optional, Union, Dict, Any

# Add parent directory to path to import analysis modules
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

class EnhancedAnalyzer:
    """Wrapper for enhanced analysis functionality"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize enhanced analyzer"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
        
        # Import and initialize the actual analyzer class
        try:
            from enhanced_analysis import EnhancedAnalyzer as EAnalyzer
            self.analyzer = EAnalyzer(str(self.base_dir))
        except ImportError as e:
            print(f"Warning: Could not import enhanced analyzer: {e}")
            self.analyzer = None
    
    def load_data(self):
        """Load simulation data"""
        if self.analyzer:
            return self.analyzer.load_data()
        return False
    
    def analyze_time_series(self):
        """Analyze time series data"""
        if self.analyzer:
            return self.analyzer.analyze_time_series()
        return {}
    
    def analyze_cellular_states(self):
        """Analyze cellular state distributions"""
        if self.analyzer:
            return self.analyzer.analyze_cellular_states()
        return {}
    
    def analyze_treatment_efficacy(self):
        """Analyze treatment efficacy"""
        if self.analyzer:
            return self.analyzer.analyze_treatment_efficacy()
        return {}
    
    def calculate_synergy_index(self):
        """Calculate synergy index for combination treatments"""
        if self.analyzer:
            return self.analyzer.calculate_synergy_index()
        return {}
    
    def generate_statistical_summary(self):
        """Generate comprehensive statistical summary"""
        if self.analyzer:
            return self.analyzer.generate_statistical_summary()
        return {}
    
    def create_analysis_plots(self):
        """Create analysis plots"""
        if self.analyzer:
            return self.analyzer.create_enhanced_time_series_analysis()
        return False
    
    def run_complete_analysis(self):
        """Run complete analysis pipeline"""
        if self.analyzer:
            return self.analyzer.run_complete_enhanced_analysis()
        return False

class StatisticalAnalyzer:
    """Statistical analysis utilities"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize statistical analyzer"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
    
    def perform_ttest(self, group1_data, group2_data):
        """Perform t-test between two groups"""
        try:
            from scipy import stats
            statistic, p_value = stats.ttest_ind(group1_data, group2_data)
            return {
                'statistic': statistic,
                'p_value': p_value,
                'significant': p_value < 0.05
            }
        except ImportError:
            print("scipy not available for statistical tests")
            return None
    
    def perform_anova(self, *groups):
        """Perform one-way ANOVA"""
        try:
            from scipy import stats
            statistic, p_value = stats.f_oneway(*groups)
            return {
                'statistic': statistic,
                'p_value': p_value,
                'significant': p_value < 0.05
            }
        except ImportError:
            print("scipy not available for statistical tests")
            return None
    
    def calculate_effect_size(self, group1_data, group2_data):
        """Calculate Cohen's d effect size"""
        try:
            import numpy as np
            
            mean1, mean2 = np.mean(group1_data), np.mean(group2_data)
            std1, std2 = np.std(group1_data), np.std(group2_data)
            n1, n2 = len(group1_data), len(group2_data)
            
            # Pooled standard deviation
            pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))
            
            # Cohen's d
            cohens_d = (mean1 - mean2) / pooled_std
            
            return {
                'cohens_d': cohens_d,
                'effect_size': 'small' if abs(cohens_d) < 0.5 else 'medium' if abs(cohens_d) < 0.8 else 'large'
            }
        except ImportError:
            print("numpy not available for effect size calculation")
            return None

class DataProcessor:
    """Data processing utilities"""
    
    def __init__(self, base_dir: Optional[Union[str, Path]] = None):
        """Initialize data processor"""
        if base_dir is None:
            base_dir = Path(__file__).parent.parent
        self.base_dir = Path(base_dir)
    
    def load_csv_data(self, filepath: Union[str, Path]):
        """Load CSV data"""
        try:
            import pandas as pd
            return pd.read_csv(filepath)
        except ImportError:
            print("pandas not available for CSV loading")
            return None
    
    def load_xml_data(self, filepath: Union[str, Path]):
        """Load PhysiCell XML data"""
        try:
            import xml.etree.ElementTree as ET
            tree = ET.parse(filepath)
            root = tree.getroot()
            
            # Extract basic information
            data = {
                'current_time': root.find('.//current_time'),
                'cell_count': len(root.findall('.//cell')),
                'cells': []
            }
            
            # Extract cell data
            for cell in root.findall('.//cell'):
                cell_data = {
                    'id': cell.get('ID'),
                    'type': cell.find('.//cell_type'),
                    'position': {
                        'x': cell.find('.//position/x'),
                        'y': cell.find('.//position/y'),
                        'z': cell.find('.//position/z')
                    }
                }
                data['cells'].append(cell_data)
            
            return data
        except ImportError:
            print("xml.etree not available for XML loading")
            return None
    
    def clean_data(self, data):
        """Clean and preprocess data"""
        try:
            import pandas as pd
            import numpy as np
            
            if isinstance(data, pd.DataFrame):
                # Remove NaN values
                data = data.dropna()
                
                # Remove outliers (beyond 3 standard deviations)
                numeric_columns = data.select_dtypes(include=[np.number]).columns
                for col in numeric_columns:
                    mean = data[col].mean()
                    std = data[col].std()
                    data = data[np.abs(data[col] - mean) <= 3 * std]
                
                return data
            
            return data
        except ImportError:
            print("pandas/numpy not available for data cleaning")
            return data
    
    def normalize_data(self, data, method='minmax'):
        """Normalize data"""
        try:
            import pandas as pd
            import numpy as np
            from sklearn.preprocessing import MinMaxScaler, StandardScaler
            
            if isinstance(data, pd.DataFrame):
                numeric_columns = data.select_dtypes(include=[np.number]).columns
                
                if method == 'minmax':
                    scaler = MinMaxScaler()
                elif method == 'standard':
                    scaler = StandardScaler()
                else:
                    return data
                
                data[numeric_columns] = scaler.fit_transform(data[numeric_columns])
                return data
            
            return data
        except ImportError:
            print("sklearn not available for data normalization")
            return data

def main():
    """Command-line interface for analysis"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Liver Fibrosis Analysis')
    parser.add_argument('--enhanced', action='store_true', help='Run enhanced analysis')
    parser.add_argument('--statistical', action='store_true', help='Run statistical analysis')
    parser.add_argument('--all', action='store_true', help='Run all analysis modules')
    
    args = parser.parse_args()
    
    success = True
    
    if args.enhanced or args.all:
        print("Running enhanced analysis...")
        analyzer = EnhancedAnalyzer()
        success &= analyzer.run_complete_analysis()
    
    if args.statistical or args.all:
        print("Running statistical analysis...")
        # Statistical analysis would be implemented here
        pass
    
    if not any([args.enhanced, args.statistical, args.all]):
        print("No analysis specified. Use --help for options.")
        success = False
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
