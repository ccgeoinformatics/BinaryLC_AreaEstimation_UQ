"""
Script: Multiple Binary Land Cover Map Error Analysis, Accuracy Assessment, Area Estimation, and Uncertainty Quantification

Description:
    This script performs error analysis, accuracy assessment, unbiased area estimation and
	uncertainty quantification for multiple binary land cover classification maps based on an error
	matrix and mapped pixel counts. It calculates several key metrics 
    to evaluate classification accuracy, including user’s accuracy, producer’s accuracy, 
    overall accuracy, error-adjusted area estimates, and confidence intervals.

    The methodology follows Olofsson et al. (2013):
        Olofsson, P., Foody, G. M., Stehman, S. V., & Woodcock, C. E. (2013). 
        Making better use of accuracy data in land change studies: 
        Estimating area and change accuracy. Remote Sensing of Environment, 129, 122-131.
		
	The formulas for the variance and standard errors of user's and producer's accuracy
	and overall accuracy are based on Olofsson et al. (2014):
		Olofsson, P., Foody, G. M., Herold, M., Stehman, S. V., Woodcock, C. E.,
		& Wulder, M. A. (2014). Good practices for estimating area and assessing accuracy
		of land change. Remote sensing of Environment, 148, 42-57.
		https://doi.org/10.1016/j.rse.2014.02.015

    Each land cover map’s data is read from a CSV file containing fields for `lcmap`, error matrix values,
    mapped pixel counts, pixel size, and unit.
	
	CSV Structure:
        lcmap, n11, n12, n21, n22, mapped_pixels_class_1, mapped_pixels_class_0, pixel_size, unit

    Example Usage:
        python multiple_binary_lc_accuracy_areaEstimation_uq.py

Author:
    Jojene R. Santillan
    Institute of Photogrammetry and GeoInformation (IPI), Leibniz University Hannover, Germany 
    & Caraga Center for Geo-Informatics & Department of Geodetic Engineering, College of Engineering and Geosciences, 
    Caraga State University, Butuan City, Philippines
    santillan@ipi.uni-hannover.de, jrsantillan@carsu.edu.ph
    4 November 2024
"""

import numpy as np
import pandas as pd
import scipy.stats as st

def calculate_weights(mapped_pixels):
    total_pixels = np.sum(mapped_pixels, dtype=np.float64)
    weights = mapped_pixels / total_pixels
    return weights.astype(np.float64)

def convert_to_area_proportion(error_matrix, weights):
    area_proportion_matrix = np.zeros_like(error_matrix, dtype=np.float64)
    for i in range(error_matrix.shape[0]):
        row_total = error_matrix[i, :].sum()
        area_proportion_matrix[i, :] = weights[i] * (error_matrix[i, :] / row_total) if row_total != 0 else np.nan
    return area_proportion_matrix

def calculate_accuracy_metrics(area_proportion_matrix, error_matrix, weights, confidence_level=0.95):
    q = area_proportion_matrix.shape[0]
    user_accuracy = []
    user_accuracy_se = []
    user_accuracy_ci_value = []
    
    producer_accuracy = []
    producer_accuracy_se = []
    producer_accuracy_ci_value = []

    z_score = st.norm.ppf(1 - (1 - confidence_level) / 2)

    # User's accuracy calculations
    for i in range(q):
        p_ii = area_proportion_matrix[i, i]
        p_i_dot = area_proportion_matrix[i, :].sum()
        U_i = p_ii / p_i_dot if p_i_dot != 0 else np.nan
        user_accuracy.append(U_i)

        # Variance, SE, and CI for user's accuracy
        n_i_dot = error_matrix[i, :].sum()
        if n_i_dot > 1:
            variance_U_i = U_i * (1 - U_i) / (n_i_dot - 1)
            se_U_i = np.sqrt(variance_U_i)
            user_accuracy_se.append(se_U_i)
            user_accuracy_ci_value.append(z_score * se_U_i)
        else:
            user_accuracy_se.append(np.nan)
            user_accuracy_ci_value.append(np.nan)

    # Producer's accuracy calculations
    for j in range(q):
        p_jj = area_proportion_matrix[j, j]
        p_dot_j = area_proportion_matrix[:, j].sum()
        P_j = p_jj / p_dot_j if p_dot_j != 0 else np.nan
        producer_accuracy.append(P_j)

        # Variance, SE, and CI for producer's accuracy
        n_dot_j = error_matrix[:, j].sum()
        if n_dot_j > 1:
            variance_P_j = P_j * (1 - P_j) / (n_dot_j - 1)
            se_P_j = np.sqrt(variance_P_j)
            producer_accuracy_se.append(se_P_j)
            producer_accuracy_ci_value.append(z_score * se_P_j)
        else:
            producer_accuracy_se.append(np.nan)
            producer_accuracy_ci_value.append(np.nan)

    # Overall accuracy calculations
    overall_accuracy = np.trace(area_proportion_matrix)
    overall_accuracy_variance = sum(
        weights[i] ** 2 * user_accuracy[i] * (1 - user_accuracy[i]) / (error_matrix[i, :].sum() - 1)
        for i in range(q) if error_matrix[i, :].sum() > 1
    )
    overall_accuracy_se = np.sqrt(overall_accuracy_variance)
    overall_accuracy_ci_value = z_score * overall_accuracy_se

    return {
        "user_accuracy": user_accuracy,
        "user_accuracy_se": user_accuracy_se,
        "user_accuracy_ci_value": user_accuracy_ci_value,
        "producer_accuracy": producer_accuracy,
        "producer_accuracy_se": producer_accuracy_se,
        "producer_accuracy_ci_value": producer_accuracy_ci_value,
        "overall_accuracy": overall_accuracy,
        "overall_accuracy_se": overall_accuracy_se,
        "overall_accuracy_ci_value": overall_accuracy_ci_value
    }

def calculate_error_adjusted_area(area_proportion_matrix, total_area):
    adjusted_areas = total_area * area_proportion_matrix.sum(axis=0)
    return adjusted_areas

def calculate_standard_error_and_ci(error_matrix, weights, area_adjusted_estimates, total_area, confidence_level=0.95):
    q = error_matrix.shape[0]
    standard_errors = []

    for j in range(q):
        se_j_squared = 0
        for i in range(q):
            n_ij = error_matrix[i, j]
            n_i_dot = error_matrix[i, :].sum()
            if n_i_dot > 1:
                se_j_squared += weights[i] ** 2 * (n_ij / n_i_dot) * (1 - n_ij / n_i_dot) / (n_i_dot - 1)
        standard_errors.append(np.sqrt(se_j_squared) * total_area)

    z_score = st.norm.ppf(1 - (1 - confidence_level) / 2)
    ci_values = [z_score * se for se in standard_errors]

    return {
        "standard_errors": standard_errors,
        "ci_values": ci_values
    }

def read_error_matrix_from_csv(filepath):
    return pd.read_csv(filepath)

def process_land_cover_map(row):
    # Read values from the row
    unit = row['unit']  # Read the unit from the CSV

    # Define error matrix and mapped pixels from row data
    error_matrix = np.array([
        [row['n11'], row['n12']],
        [row['n21'], row['n22']]
    ], dtype=np.float64)
    mapped_pixels = np.array([row['mapped_pixels_class_1'], row['mapped_pixels_class_0']], dtype=np.float64)
    pixel_size = row['pixel_size']

    # Calculate weights, area proportions, and other metrics
    weights = calculate_weights(mapped_pixels)
    area_proportion_matrix = convert_to_area_proportion(error_matrix, weights)
    pixel_area = pixel_size ** 2
    total_area = np.sum(mapped_pixels) * pixel_area

    accuracy_metrics = calculate_accuracy_metrics(area_proportion_matrix, error_matrix, weights)
    error_adjusted_area = calculate_error_adjusted_area(area_proportion_matrix, total_area)
    standard_error_and_ci = calculate_standard_error_and_ci(error_matrix, weights, error_adjusted_area, total_area)

    return {
        "lcmap": row['lcmap'],
        "accuracy_metrics": accuracy_metrics,
        "error_adjusted_area": error_adjusted_area,
        "standard_error_and_ci": standard_error_and_ci,
        "unit": unit  # Pass the unit to be used in column naming
    }

def main():
    filepath = "input.csv"  # Update with the actual CSV file path
    df = read_error_matrix_from_csv(filepath)
    
    all_results = []
    for _, row in df.iterrows():
        result = process_land_cover_map(row)
        all_results.append(result)

    for result in all_results:
        # Retrieve the unit for each map
        unit = result["unit"]

        # Create DataFrame for results with units dynamically added to column names
        df_result = pd.DataFrame({
            "Class": ["Class_1", "Class_0"],
            "User_Accuracy": result["accuracy_metrics"]["user_accuracy"],
            "User_Accuracy_SE": result["accuracy_metrics"]["user_accuracy_se"],
            "User_Accuracy_CI_Value": result["accuracy_metrics"]["user_accuracy_ci_value"],
            "Producer_Accuracy": result["accuracy_metrics"]["producer_accuracy"],
            "Producer_Accuracy_SE": result["accuracy_metrics"]["producer_accuracy_se"],
            "Producer_Accuracy_CI_Value": result["accuracy_metrics"]["producer_accuracy_ci_value"],
            "Overall_Accuracy": [result["accuracy_metrics"]["overall_accuracy"]] * 2,
            "Overall_Accuracy_SE": [result["accuracy_metrics"]["overall_accuracy_se"]] * 2,
            "Overall_Accuracy_CI_Value": [result["accuracy_metrics"]["overall_accuracy_ci_value"]] * 2,
            f"Error_Adjusted_Area [{unit}]": result["error_adjusted_area"],
            f"Standard_Error [{unit}]": result["standard_error_and_ci"]["standard_errors"],
            f"CI_Value [{unit}]": result["standard_error_and_ci"]["ci_values"]
        })
        
        # Save to CSV
        df_result.to_csv(f"{result['lcmap']}_results.csv", index=False)

# Run the main function
if __name__ == "__main__":
    main()
