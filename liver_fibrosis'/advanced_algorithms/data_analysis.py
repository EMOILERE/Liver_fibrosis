"""
Data Analysis Module

Implements advanced data analysis techniques for high-dimensional biological
data including dimensionality reduction, clustering, correlation analysis,
and statistical testing.

Key Features:
- t-SNE and PCA for dimensionality reduction
- Gaussian mixture modeling for clustering
- Spatial correlation analysis with Moran's I
- Time series analysis and forecasting
- Comprehensive statistical test suite
- Information theory metrics
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Callable, Any, Union
from dataclasses import dataclass
from enum import Enum
import warnings
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score


class DistanceMetric(Enum):
    """Distance metrics for analysis."""
    EUCLIDEAN = "euclidean"
    MANHATTAN = "manhattan"
    COSINE = "cosine"
    CORRELATION = "correlation"
    JACCARD = "jaccard"


class ClusteringMethod(Enum):
    """Clustering methods."""
    KMEANS = "kmeans"
    GAUSSIAN_MIXTURE = "gaussian_mixture"
    HIERARCHICAL = "hierarchical"
    DBSCAN = "dbscan"


@dataclass
class AnalysisResult:
    """Container for analysis results."""
    method: str
    parameters: Dict[str, Any]
    results: Dict[str, Any]
    metrics: Dict[str, float]
    timestamp: float
    
    def __post_init__(self):
        if self.timestamp == 0:
            import time
            self.timestamp = time.time()


class DimensionalityReducer:
    """
    Advanced dimensionality reduction techniques for biological data.
    """
    
    def __init__(self):
        """Initialize dimensionality reducer."""
        self.fitted_models = {}
        self.scalers = {}
        
    def fit_pca(self, data: np.ndarray, n_components: Optional[int] = None,
                explained_variance_threshold: float = 0.95) -> Dict[str, Any]:
        """
        Fit PCA model to data.
        
        Args:
            data: Input data matrix (samples x features)
            n_components: Number of components (auto if None)
            explained_variance_threshold: Variance threshold for auto selection
            
        Returns:
            PCA results and model
        """
        # Standardize data
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)
        self.scalers['pca'] = scaler
        
        # Determine number of components
        if n_components is None:
            # Find components needed for variance threshold
            pca_temp = PCA()
            pca_temp.fit(data_scaled)
            
            cumsum_var = np.cumsum(pca_temp.explained_variance_ratio_)
            n_components = np.argmax(cumsum_var >= explained_variance_threshold) + 1
        
        # Fit final PCA model
        pca = PCA(n_components=n_components)
        transformed_data = pca.fit_transform(data_scaled)
        
        self.fitted_models['pca'] = pca
        
        return {
            'transformed_data': transformed_data,
            'explained_variance_ratio': pca.explained_variance_ratio_,
            'cumulative_variance': np.cumsum(pca.explained_variance_ratio_),
            'components': pca.components_,
            'n_components': n_components,
            'total_variance_explained': np.sum(pca.explained_variance_ratio_)
        }
    
    def fit_tsne(self, data: np.ndarray, n_components: int = 2,
                perplexity: float = 30.0, learning_rate: float = 200.0,
                n_iter: int = 1000) -> Dict[str, Any]:
        """
        Fit t-SNE model to data.
        
        Args:
            data: Input data matrix
            n_components: Number of embedding dimensions
            perplexity: t-SNE perplexity parameter
            learning_rate: Learning rate for optimization
            n_iter: Number of iterations
            
        Returns:
            t-SNE results
        """
        # Standardize data
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)
        self.scalers['tsne'] = scaler
        
        # Fit t-SNE
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            n_iter=n_iter,
            random_state=42
        )
        
        transformed_data = tsne.fit_transform(data_scaled)
        
        return {
            'transformed_data': transformed_data,
            'kl_divergence': tsne.kl_divergence_,
            'n_iter_final': tsne.n_iter_,
            'perplexity': perplexity,
            'learning_rate': learning_rate
        }
    
    def calculate_intrinsic_dimensionality(self, data: np.ndarray,
                                         method: str = "mle") -> float:
        """
        Estimate intrinsic dimensionality of data.
        
        Args:
            data: Input data matrix
            method: Estimation method ("mle", "correlation", "pca")
            
        Returns:
            Estimated intrinsic dimensionality
        """
        if method == "pca":
            # Use PCA eigenvalue decay
            pca = PCA()
            pca.fit(data)
            
            # Find elbow in eigenvalue spectrum
            eigenvalues = pca.explained_variance_
            ratios = eigenvalues[:-1] / eigenvalues[1:]
            
            # Simple elbow detection
            intrinsic_dim = np.argmax(ratios) + 1
            
        elif method == "correlation":
            # Use correlation dimension
            distances = pdist(data)
            
            # Count pairs within different distance thresholds
            thresholds = np.logspace(-2, 0, 20)
            counts = []
            
            for threshold in thresholds:
                count = np.sum(distances < threshold)
                counts.append(count)
            
            # Fit power law: log(count) ~ d * log(threshold)
            log_thresholds = np.log(thresholds[1:])  # Avoid log(0)
            log_counts = np.log(np.array(counts[1:]) + 1)  # Avoid log(0)
            
            if len(log_thresholds) > 1:
                slope, _, _, _, _ = stats.linregress(log_thresholds, log_counts)
                intrinsic_dim = slope
            else:
                intrinsic_dim = data.shape[1]
                
        else:  # MLE method (simplified)
            # Maximum likelihood estimation using nearest neighbor distances
            from sklearn.neighbors import NearestNeighbors
            
            k = min(10, data.shape[0] - 1)
            nbrs = NearestNeighbors(n_neighbors=k+1).fit(data)
            distances, _ = nbrs.kneighbors(data)
            
            # Use k-th nearest neighbor distances
            kth_distances = distances[:, k]
            
            # Simple MLE estimate
            if np.std(kth_distances) > 0:
                intrinsic_dim = np.mean(kth_distances) / np.std(kth_distances)
            else:
                intrinsic_dim = 1.0
        
        return max(1.0, min(float(data.shape[1]), intrinsic_dim))


class TSNEAnalyzer:
    """
    Specialized t-SNE analyzer with advanced features.
    """
    
    def __init__(self):
        """Initialize t-SNE analyzer."""
        self.embeddings = {}
        self.perplexity_analysis = {}
        
    def perplexity_analysis(self, data: np.ndarray,
                           perplexity_range: Tuple[float, float] = (5, 50),
                           n_trials: int = 10) -> Dict[str, Any]:
        """
        Analyze optimal perplexity for t-SNE.
        
        Args:
            data: Input data matrix
            perplexity_range: Range of perplexity values to test
            n_trials: Number of trials per perplexity
            
        Returns:
            Perplexity analysis results
        """
        perplexities = np.logspace(
            np.log10(perplexity_range[0]),
            np.log10(perplexity_range[1]),
            n_trials
        )
        
        results = {
            'perplexities': perplexities,
            'kl_divergences': [],
            'silhouette_scores': [],
            'embeddings': []
        }
        
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)
        
        for perplexity in perplexities:
            tsne = TSNE(
                n_components=2,
                perplexity=perplexity,
                random_state=42,
                n_iter=1000
            )
            
            embedding = tsne.fit_transform(data_scaled)
            
            # Calculate quality metrics
            kl_div = tsne.kl_divergence_
            
            # Silhouette score (requires clustering)
            from sklearn.cluster import KMeans
            n_clusters = min(8, max(2, int(np.sqrt(data.shape[0] / 2))))
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            labels = kmeans.fit_predict(embedding)
            
            if len(np.unique(labels)) > 1:
                sil_score = silhouette_score(embedding, labels)
            else:
                sil_score = 0.0
            
            results['kl_divergences'].append(kl_div)
            results['silhouette_scores'].append(sil_score)
            results['embeddings'].append(embedding)
        
        # Find optimal perplexity
        # Normalize metrics and combine
        kl_norm = 1.0 - np.array(results['kl_divergences']) / np.max(results['kl_divergences'])
        sil_norm = np.array(results['silhouette_scores']) / np.max(results['silhouette_scores'])
        
        combined_score = 0.6 * sil_norm + 0.4 * kl_norm
        optimal_idx = np.argmax(combined_score)
        
        results['optimal_perplexity'] = perplexities[optimal_idx]
        results['optimal_embedding'] = results['embeddings'][optimal_idx]
        
        return results
    
    def guided_tsne(self, data: np.ndarray, labels: np.ndarray,
                   guidance_weight: float = 0.1) -> np.ndarray:
        """
        Perform guided t-SNE using known labels.
        
        Args:
            data: Input data matrix
            labels: Known labels for guidance
            guidance_weight: Weight for guidance term
            
        Returns:
            Guided t-SNE embedding
        """
        # This is a simplified implementation
        # Full guided t-SNE would require custom loss function
        
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)
        
        # Standard t-SNE
        tsne = TSNE(n_components=2, random_state=42)
        embedding = tsne.fit_transform(data_scaled)
        
        # Post-process embedding to respect labels (simplified)
        unique_labels = np.unique(labels)
        
        for label in unique_labels:
            mask = labels == label
            if np.sum(mask) > 1:
                # Move points with same label closer together
                centroid = np.mean(embedding[mask], axis=0)
                embedding[mask] += guidance_weight * (centroid - embedding[mask])
        
        return embedding


class PCAAnalyzer:
    """
    Advanced PCA analyzer with feature interpretation.
    """
    
    def __init__(self):
        """Initialize PCA analyzer."""
        self.loadings_analysis = {}
        
    def analyze_loadings(self, pca_model, feature_names: List[str],
                        n_components: int = 5) -> Dict[str, Any]:
        """
        Analyze PCA loadings for feature interpretation.
        
        Args:
            pca_model: Fitted PCA model
            feature_names: Names of original features
            n_components: Number of components to analyze
            
        Returns:
            Loadings analysis results
        """
        components = pca_model.components_[:n_components]
        
        results = {
            'component_contributions': {},
            'feature_importance': {},
            'top_features_per_component': {}
        }
        
        for i, component in enumerate(components):
            # Feature contributions to this component
            contributions = dict(zip(feature_names, component))
            results['component_contributions'][f'PC{i+1}'] = contributions
            
            # Top contributing features
            sorted_features = sorted(contributions.items(), 
                                   key=lambda x: abs(x[1]), reverse=True)
            results['top_features_per_component'][f'PC{i+1}'] = sorted_features[:10]
        
        # Overall feature importance across components
        feature_importance = {}
        for feature_name in feature_names:
            importance = 0.0
            for i in range(n_components):
                loading = pca_model.components_[i, feature_names.index(feature_name)]
                variance_weight = pca_model.explained_variance_ratio_[i]
                importance += abs(loading) * variance_weight
            
            feature_importance[feature_name] = importance
        
        results['feature_importance'] = dict(
            sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)
        )
        
        return results
    
    def biplot_data(self, pca_model, transformed_data: np.ndarray,
                   feature_names: List[str]) -> Dict[str, np.ndarray]:
        """
        Generate data for PCA biplot visualization.
        
        Args:
            pca_model: Fitted PCA model
            transformed_data: PCA-transformed data
            feature_names: Original feature names
            
        Returns:
            Biplot data for visualization
        """
        # Sample scores (first two components)
        scores = transformed_data[:, :2]
        
        # Feature loadings (first two components)
        loadings = pca_model.components_[:2].T
        
        # Scale loadings for visualization
        loading_scale = np.max(np.abs(scores)) / np.max(np.abs(loadings))
        loadings_scaled = loadings * loading_scale * 0.8
        
        return {
            'scores': scores,
            'loadings': loadings_scaled,
            'feature_names': feature_names,
            'explained_variance': pca_model.explained_variance_ratio_[:2]
        }


class GaussianMixtureAnalyzer:
    """
    Gaussian mixture model analyzer for clustering and density estimation.
    """
    
    def __init__(self):
        """Initialize Gaussian mixture analyzer."""
        self.fitted_models = {}
        
    def fit_optimal_gmm(self, data: np.ndarray, max_components: int = 10,
                       covariance_types: List[str] = None) -> Dict[str, Any]:
        """
        Fit optimal Gaussian mixture model using model selection.
        
        Args:
            data: Input data matrix
            max_components: Maximum number of components to test
            covariance_types: Covariance types to test
            
        Returns:
            Optimal GMM results
        """
        if covariance_types is None:
            covariance_types = ['full', 'tied', 'diag', 'spherical']
        
        # Standardize data
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)
        
        best_model = None
        best_bic = np.inf
        best_params = {}
        
        results = {
            'model_comparison': [],
            'bic_scores': {},
            'aic_scores': {}
        }
        
        for cov_type in covariance_types:
            bic_scores = []
            aic_scores = []
            
            for n_components in range(1, max_components + 1):
                try:
                    gmm = GaussianMixture(
                        n_components=n_components,
                        covariance_type=cov_type,
                        random_state=42,
                        max_iter=200
                    )
                    
                    gmm.fit(data_scaled)
                    
                    bic = gmm.bic(data_scaled)
                    aic = gmm.aic(data_scaled)
                    
                    bic_scores.append(bic)
                    aic_scores.append(aic)
                    
                    results['model_comparison'].append({
                        'n_components': n_components,
                        'covariance_type': cov_type,
                        'bic': bic,
                        'aic': aic,
                        'converged': gmm.converged_
                    })
                    
                    if bic < best_bic and gmm.converged_:
                        best_bic = bic
                        best_model = gmm
                        best_params = {
                            'n_components': n_components,
                            'covariance_type': cov_type
                        }
                
                except Exception as e:
                    warnings.warn(f"GMM fitting failed for {n_components} components, "
                                f"{cov_type} covariance: {str(e)}")
            
            results['bic_scores'][cov_type] = bic_scores
            results['aic_scores'][cov_type] = aic_scores
        
        if best_model is not None:
            # Analyze best model
            labels = best_model.predict(data_scaled)
            probabilities = best_model.predict_proba(data_scaled)
            
            results.update({
                'best_model': best_model,
                'best_params': best_params,
                'labels': labels,
                'probabilities': probabilities,
                'means': best_model.means_,
                'covariances': best_model.covariances_,
                'weights': best_model.weights_,
                'scaler': scaler
            })
            
            # Calculate clustering metrics
            if len(np.unique(labels)) > 1:
                results['silhouette_score'] = silhouette_score(data_scaled, labels)
            
            self.fitted_models['optimal'] = best_model
        
        return results
    
    def anomaly_detection(self, data: np.ndarray, contamination: float = 0.1) -> Dict[str, Any]:
        """
        Perform anomaly detection using GMM.
        
        Args:
            data: Input data matrix
            contamination: Expected fraction of anomalies
            
        Returns:
            Anomaly detection results
        """
        # Fit GMM
        gmm_results = self.fit_optimal_gmm(data)
        
        if 'best_model' not in gmm_results:
            return {'error': 'Failed to fit GMM model'}
        
        model = gmm_results['best_model']
        scaler = gmm_results['scaler']
        data_scaled = scaler.transform(data)
        
        # Calculate log-likelihood scores
        log_likelihood = model.score_samples(data_scaled)
        
        # Determine threshold for anomalies
        threshold = np.percentile(log_likelihood, contamination * 100)
        
        # Identify anomalies
        is_anomaly = log_likelihood < threshold
        
        return {
            'log_likelihood': log_likelihood,
            'threshold': threshold,
            'is_anomaly': is_anomaly,
            'anomaly_indices': np.where(is_anomaly)[0],
            'n_anomalies': np.sum(is_anomaly)
        }


class SpatialCorrelationAnalyzer:
    """
    Analyzer for spatial correlation patterns in biological data.
    """
    
    def __init__(self):
        """Initialize spatial correlation analyzer."""
        pass
    
    def morans_i(self, values: np.ndarray, coordinates: np.ndarray,
                distance_threshold: Optional[float] = None) -> Dict[str, float]:
        """
        Calculate Moran's I spatial autocorrelation statistic.
        
        Args:
            values: Variable values at each location
            coordinates: Spatial coordinates
            distance_threshold: Maximum distance for neighbors
            
        Returns:
            Moran's I statistics
        """
        n = len(values)
        
        # Calculate distance matrix
        distances = squareform(pdist(coordinates))
        
        # Create spatial weights matrix
        if distance_threshold is None:
            # Use average nearest neighbor distance
            min_distances = np.partition(distances, 1, axis=1)[:, 1]
            distance_threshold = np.mean(min_distances) * 2
        
        weights = (distances <= distance_threshold) & (distances > 0)
        weights = weights.astype(float)
        
        # Row-normalize weights
        row_sums = np.sum(weights, axis=1)
        weights[row_sums > 0] = weights[row_sums > 0] / row_sums[row_sums > 0, np.newaxis]
        
        # Calculate Moran's I
        mean_value = np.mean(values)
        deviations = values - mean_value
        
        numerator = 0.0
        denominator = np.sum(deviations**2)
        
        for i in range(n):
            for j in range(n):
                numerator += weights[i, j] * deviations[i] * deviations[j]
        
        if denominator > 0:
            morans_i = (n / np.sum(weights)) * (numerator / denominator)
        else:
            morans_i = 0.0
        
        # Calculate expected value and variance (under null hypothesis)
        expected_i = -1.0 / (n - 1)
        
        # Simplified variance calculation
        variance_i = (n**2 - 3*n + 3) / ((n - 1) * (n - 2) * (n - 3)) - expected_i**2
        
        # Z-score and p-value
        if variance_i > 0:
            z_score = (morans_i - expected_i) / np.sqrt(variance_i)
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        else:
            z_score = 0.0
            p_value = 1.0
        
        return {
            'morans_i': morans_i,
            'expected_i': expected_i,
            'variance_i': variance_i,
            'z_score': z_score,
            'p_value': p_value,
            'distance_threshold': distance_threshold
        }
    
    def local_morans_i(self, values: np.ndarray, coordinates: np.ndarray,
                      distance_threshold: Optional[float] = None) -> Dict[str, np.ndarray]:
        """
        Calculate local Moran's I (LISA) statistics.
        
        Args:
            values: Variable values at each location
            coordinates: Spatial coordinates
            distance_threshold: Maximum distance for neighbors
            
        Returns:
            Local Moran's I statistics for each location
        """
        n = len(values)
        
        # Calculate distance matrix and weights
        distances = squareform(pdist(coordinates))
        
        if distance_threshold is None:
            min_distances = np.partition(distances, 1, axis=1)[:, 1]
            distance_threshold = np.mean(min_distances) * 2
        
        weights = (distances <= distance_threshold) & (distances > 0)
        weights = weights.astype(float)
        
        # Calculate local Moran's I for each location
        mean_value = np.mean(values)
        deviations = values - mean_value
        variance = np.var(values)
        
        local_i = np.zeros(n)
        
        for i in range(n):
            if np.sum(weights[i]) > 0:
                weighted_sum = np.sum(weights[i] * deviations)
                local_i[i] = (deviations[i] / variance) * weighted_sum
        
        return {
            'local_morans_i': local_i,
            'deviations': deviations,
            'weights_matrix': weights,
            'distance_threshold': distance_threshold
        }
    
    def spatial_clustering_analysis(self, coordinates: np.ndarray,
                                  cluster_labels: np.ndarray) -> Dict[str, Any]:
        """
        Analyze spatial clustering patterns.
        
        Args:
            coordinates: Spatial coordinates
            cluster_labels: Cluster assignments
            
        Returns:
            Spatial clustering analysis results
        """
        unique_labels = np.unique(cluster_labels)
        n_clusters = len(unique_labels)
        
        results = {
            'cluster_centroids': {},
            'cluster_dispersions': {},
            'inter_cluster_distances': np.zeros((n_clusters, n_clusters)),
            'spatial_compactness': {}
        }
        
        # Calculate cluster centroids and dispersions
        for i, label in enumerate(unique_labels):
            mask = cluster_labels == label
            cluster_coords = coordinates[mask]
            
            if len(cluster_coords) > 0:
                centroid = np.mean(cluster_coords, axis=0)
                results['cluster_centroids'][label] = centroid
                
                # Calculate dispersion (average distance from centroid)
                distances_to_centroid = np.linalg.norm(cluster_coords - centroid, axis=1)
                dispersion = np.mean(distances_to_centroid)
                results['cluster_dispersions'][label] = dispersion
                
                # Spatial compactness (ratio of area to perimeter squared)
                if len(cluster_coords) > 2:
                    from scipy.spatial import ConvexHull
                    try:
                        hull = ConvexHull(cluster_coords)
                        area = hull.volume  # In 2D, volume is area
                        perimeter = hull.area  # In 2D, area is perimeter
                        compactness = 4 * np.pi * area / (perimeter**2) if perimeter > 0 else 0
                        results['spatial_compactness'][label] = compactness
                    except:
                        results['spatial_compactness'][label] = 0.0
        
        # Calculate inter-cluster distances
        centroids = list(results['cluster_centroids'].values())
        if len(centroids) > 1:
            centroid_distances = squareform(pdist(centroids))
            results['inter_cluster_distances'] = centroid_distances
        
        return results


class TimeSeriesAnalyzer:
    """
    Time series analysis for temporal biological data.
    """
    
    def __init__(self):
        """Initialize time series analyzer."""
        pass
    
    def trend_analysis(self, time_series: np.ndarray, time_points: np.ndarray) -> Dict[str, Any]:
        """
        Analyze trends in time series data.
        
        Args:
            time_series: Time series values
            time_points: Time points
            
        Returns:
            Trend analysis results
        """
        # Linear trend
        slope, intercept, r_value, p_value, std_err = stats.linregress(time_points, time_series)
        
        # Polynomial trends
        poly_degrees = [2, 3]
        polynomial_fits = {}
        
        for degree in poly_degrees:
            coeffs = np.polyfit(time_points, time_series, degree)
            poly_fit = np.polyval(coeffs, time_points)
            
            # R-squared for polynomial fit
            ss_res = np.sum((time_series - poly_fit)**2)
            ss_tot = np.sum((time_series - np.mean(time_series))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            polynomial_fits[f'degree_{degree}'] = {
                'coefficients': coeffs,
                'fitted_values': poly_fit,
                'r_squared': r_squared
            }
        
        # Seasonal decomposition (simplified)
        # Detect periodicity using autocorrelation
        autocorr = np.correlate(time_series, time_series, mode='full')
        autocorr = autocorr[autocorr.size // 2:]
        autocorr = autocorr / autocorr[0]  # Normalize
        
        # Find peaks in autocorrelation (potential periods)
        from scipy.signal import find_peaks
        peaks, _ = find_peaks(autocorr[1:], height=0.1)
        
        potential_periods = peaks + 1  # Add 1 because we started from index 1
        
        return {
            'linear_trend': {
                'slope': slope,
                'intercept': intercept,
                'r_value': r_value,
                'p_value': p_value,
                'std_err': std_err
            },
            'polynomial_fits': polynomial_fits,
            'autocorrelation': autocorr,
            'potential_periods': potential_periods,
            'is_stationary': self._test_stationarity(time_series)
        }
    
    def _test_stationarity(self, time_series: np.ndarray) -> Dict[str, Any]:
        """Test for stationarity using simple statistics."""
        # Split series into two halves
        n = len(time_series)
        first_half = time_series[:n//2]
        second_half = time_series[n//2:]
        
        # Compare means and variances
        mean_diff = abs(np.mean(first_half) - np.mean(second_half))
        var_ratio = np.var(first_half) / (np.var(second_half) + 1e-12)
        
        # Simple stationarity test
        is_stationary = (mean_diff < 0.1 * np.std(time_series)) and (0.5 < var_ratio < 2.0)
        
        return {
            'is_stationary': is_stationary,
            'mean_difference': mean_diff,
            'variance_ratio': var_ratio
        }
    
    def changepoint_detection(self, time_series: np.ndarray,
                            min_segment_length: int = 10) -> Dict[str, Any]:
        """
        Detect changepoints in time series.
        
        Args:
            time_series: Time series values
            min_segment_length: Minimum length of segments
            
        Returns:
            Changepoint detection results
        """
        n = len(time_series)
        changepoints = []
        
        # Simple variance-based changepoint detection
        for i in range(min_segment_length, n - min_segment_length):
            left_segment = time_series[:i]
            right_segment = time_series[i:]
            
            # Calculate variance change
            left_var = np.var(left_segment)
            right_var = np.var(right_segment)
            total_var = np.var(time_series)
            
            # Variance reduction score
            weighted_var = (len(left_segment) * left_var + len(right_segment) * right_var) / n
            variance_reduction = total_var - weighted_var
            
            if variance_reduction > 0.1 * total_var:  # Threshold
                changepoints.append({
                    'index': i,
                    'variance_reduction': variance_reduction,
                    'left_mean': np.mean(left_segment),
                    'right_mean': np.mean(right_segment)
                })
        
        # Sort by variance reduction
        changepoints.sort(key=lambda x: x['variance_reduction'], reverse=True)
        
        return {
            'changepoints': changepoints[:5],  # Top 5 changepoints
            'n_changepoints': len(changepoints)
        }


class StatisticalTestSuite:
    """
    Comprehensive statistical testing suite for biological data.
    """
    
    def __init__(self):
        """Initialize statistical test suite."""
        pass
    
    def normality_tests(self, data: np.ndarray) -> Dict[str, Any]:
        """
        Perform normality tests on data.
        
        Args:
            data: Data to test
            
        Returns:
            Normality test results
        """
        results = {}
        
        # Shapiro-Wilk test
        if len(data) <= 5000:  # Shapiro-Wilk has sample size limitations
            shapiro_stat, shapiro_p = stats.shapiro(data)
            results['shapiro_wilk'] = {
                'statistic': shapiro_stat,
                'p_value': shapiro_p,
                'is_normal': shapiro_p > 0.05
            }
        
        # Kolmogorov-Smirnov test
        ks_stat, ks_p = stats.kstest(data, 'norm', args=(np.mean(data), np.std(data)))
        results['kolmogorov_smirnov'] = {
            'statistic': ks_stat,
            'p_value': ks_p,
            'is_normal': ks_p > 0.05
        }
        
        # Anderson-Darling test
        ad_stat, ad_critical, ad_significance = stats.anderson(data, dist='norm')
        results['anderson_darling'] = {
            'statistic': ad_stat,
            'critical_values': ad_critical,
            'significance_levels': ad_significance,
            'is_normal': ad_stat < ad_critical[2]  # 5% significance level
        }
        
        return results
    
    def correlation_tests(self, x: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
        """
        Perform correlation tests between two variables.
        
        Args:
            x: First variable
            y: Second variable
            
        Returns:
            Correlation test results
        """
        results = {}
        
        # Pearson correlation
        pearson_r, pearson_p = stats.pearsonr(x, y)
        results['pearson'] = {
            'correlation': pearson_r,
            'p_value': pearson_p,
            'is_significant': pearson_p < 0.05
        }
        
        # Spearman correlation
        spearman_r, spearman_p = stats.spearmanr(x, y)
        results['spearman'] = {
            'correlation': spearman_r,
            'p_value': spearman_p,
            'is_significant': spearman_p < 0.05
        }
        
        # Kendall's tau
        kendall_tau, kendall_p = stats.kendalltau(x, y)
        results['kendall'] = {
            'correlation': kendall_tau,
            'p_value': kendall_p,
            'is_significant': kendall_p < 0.05
        }
        
        return results
    
    def group_comparison_tests(self, groups: List[np.ndarray]) -> Dict[str, Any]:
        """
        Perform group comparison tests.
        
        Args:
            groups: List of group data arrays
            
        Returns:
            Group comparison test results
        """
        results = {}
        
        if len(groups) == 2:
            # Two-group tests
            group1, group2 = groups
            
            # Independent t-test
            t_stat, t_p = stats.ttest_ind(group1, group2)
            results['t_test'] = {
                'statistic': t_stat,
                'p_value': t_p,
                'is_significant': t_p < 0.05
            }
            
            # Mann-Whitney U test
            u_stat, u_p = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            results['mann_whitney'] = {
                'statistic': u_stat,
                'p_value': u_p,
                'is_significant': u_p < 0.05
            }
            
        elif len(groups) > 2:
            # Multi-group tests
            
            # One-way ANOVA
            f_stat, f_p = stats.f_oneway(*groups)
            results['anova'] = {
                'statistic': f_stat,
                'p_value': f_p,
                'is_significant': f_p < 0.05
            }
            
            # Kruskal-Wallis test
            h_stat, h_p = stats.kruskal(*groups)
            results['kruskal_wallis'] = {
                'statistic': h_stat,
                'p_value': h_p,
                'is_significant': h_p < 0.05
            }
        
        return results
    
    def multiple_testing_correction(self, p_values: np.ndarray,
                                  method: str = "fdr_bh") -> Dict[str, np.ndarray]:
        """
        Apply multiple testing correction.
        
        Args:
            p_values: Array of p-values
            method: Correction method ("bonferroni", "fdr_bh", "fdr_by")
            
        Returns:
            Corrected p-values and significance indicators
        """
        from statsmodels.stats.multitest import multipletests
        
        rejected, p_corrected, alpha_sidak, alpha_bonf = multipletests(
            p_values, method=method, alpha=0.05
        )
        
        return {
            'p_values_corrected': p_corrected,
            'rejected': rejected,
            'alpha_sidak': alpha_sidak,
            'alpha_bonferroni': alpha_bonf,
            'method': method
        }
