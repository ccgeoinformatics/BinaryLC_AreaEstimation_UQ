"""
Script: Binary Land Cover Map Error Analysis, Accuracy Assessment, Area Estimation, and Uncertainty Quantification

Description:
    This script performs error analysis, accuracy assessment, unbiased area estimation and
	uncertainty quantification for a binary land cover classification map based on an error
	matrix and mapped pixel counts. It calculates several key metrics 
    to evaluate classification accuracy, including user’s accuracy, producer’s accuracy, 
    overall accuracy, error-adjusted area estimates, and confidence intervals.
    
    The methodology for calculating error-adjusted area estimates, accuracy metrics, and 
    confidence intervals follows the approach outlined in Olofsson et al. (2013):
        Olofsson, P., Foody, G. M., Stehman, S. V., & Woodcock, C. E. (2013). 
        Making better use of accuracy data in land change studies: 
        Estimating area and change accuracy. Remote Sensing of Environment, 129, 122-131.
        https://doi.org/10.1016/j.rse.2012.10.031
    
	The formulas for the variance and standard errors of user's and producer's accuracy
	and overall accuracy are based on Olofsson et al. (2014):
		Olofsson, P., Foody, G. M., Herold, M., Stehman, S. V., Woodcock, C. E.,
		& Wulder, M. A. (2014). Good practices for estimating area and assessing accuracy
		of land change. Remote sensing of Environment, 148, 42-57.
		https://doi.org/10.1016/j.rse.2014.02.015
	
    The script accepts an error matrix representing sample counts for two classes and 
    the number of pixels mapped for each class, then uses these to compute:
        - Class weights based on mapped pixels
        - Area proportion matrix
        - User’s and producer’s accuracy for each class
        - Overall accuracy
		- Standard errors and confidence intervals for the above accuracy metrics
        - Error-adjusted area estimates for each class
        - Standard errors and confidence intervals for the area estimates
        
    ****Note: For this script, the classes are labelled as '1' (class #1) and '0' (class #2) ****.

Key Functions:
    - calculate_weights: Computes weights (W_i) based on the number of pixels mapped to each class.
    - convert_to_area_proportion: Converts an error matrix of sample counts to area proportions using weights.
    - calculate_accuracy_metrics: Calculates user’s accuracy, producer’s accuracy, and overall accuracy, including their standard errors and confidence intervals.
    - calculate_error_adjusted_area: Computes error-adjusted area estimates for each class.
    - calculate_standard_error_and_ci: Calculates standard errors and confidence intervals for area estimates.

Input Parameters:
    - error_matrix: 2D numpy array representing the sample counts for each class.
    - mapped_pixels: Array containing the number of pixels mapped to each class.
    - pixel_size: Side length of each pixel in meters, used to calculate pixel area.
    - confidence_level: Confidence level for calculating confidence intervals (default is 0.95).

Outputs:
    - A CSV file with a summary of calculated metrics, including user’s accuracy, producer’s accuracy, 
      overall accuracy, error-adjusted area estimates, standard errors, and confidence intervals.
    - CSV files containing the original error matrix with totals and the area proportion matrix with totals.

Requirements:
    - Python 3.x
    - numpy
    - pandas
    - scipy

Example Usage:
    Run the script with the specified error matrix and mapped pixel counts:
    
        python single_binary_lc_accuracy_areaEstimation_uq.py

Author:
    Jojene R. Santillan
    Institute of Photogrammetry and GeoInformation (IPI), Leibniz University Hannover, Germany 
    & Caraga Center for Geo-Informatics & Department of Geodetic Engineering, College of Engineering and Geosciences, 
    Caraga State University, Butuan City, Philippines
    santillan@ipi.uni-hannover.de, jrsantillan@carsu.edu.ph
    28 October 2024
"""

import numpy as np
import pandas as pd
import scipy.stats as st  # Importing scipy.stats for confidence interval calculation

def calculate_weights(mapped_pixels):
    """
    Calculate weights (W_i) based on the number of pixels mapped to each class,
    with maximum precision.
    """
    total_pixels = np.sum(mapped_pixels, dtype=np.float64)
    weights = mapped_pixels / total_pixels
    return weights.astype(np.float64)

def convert_to_area_proportion(error_matrix, weights):
    """
    Convert an error matrix of sample counts to area proportions using calculated weights.
    """
    area_proportion_matrix = np.zeros_like(error_matrix, dtype=np.float64)
    for i in range(error_matrix.shape[0]):
        row_total = error_matrix[i, :].sum()
        area_proportion_matrix[i, :] = weights[i] * (error_matrix[i, :] / row_total) if row_total != 0 else np.nan
    return area_proportion_matrix

def calculate_accuracy_metrics(area_proportion_matrix, error_matrix, weights, confidence_level=0.95):
    """
    Calculate user’s accuracy, producer’s accuracy, and overall accuracy from area proportions.
    Also calculates the standard error, confidence intervals, and CI value for both user’s, producer’s, and overall accuracies.
    """
    q = area_proportion_matrix.shape[0]
    user_accuracy = []
    user_accuracy_se = []
    user_accuracy_ci_lower = []
    user_accuracy_ci_upper = []
    user_accuracy_ci_value = []

    producer_accuracy = []
    producer_accuracy_se = []
    producer_accuracy_ci_lower = []
    producer_accuracy_ci_upper = []
    producer_accuracy_ci_value = []

    z_score = st.norm.ppf(1 - (1 - confidence_level) / 2)  # z-score for the specified confidence level

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
            
            ci_lower = U_i - z_score * se_U_i
            ci_upper = U_i + z_score * se_U_i
            user_accuracy_ci_lower.append(ci_lower)
            user_accuracy_ci_upper.append(ci_upper)
            ci_value =  z_score * se_U_i
            user_accuracy_ci_value.append(ci_value)
        else:
            user_accuracy_se.append(np.nan)
            user_accuracy_ci_lower.append(np.nan)
            user_accuracy_ci_upper.append(np.nan)
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
            variance_P_j = (1 / n_dot_j ** 2) * (
                (n_dot_j ** 2 * (1 - P_j) ** 2 * U_i * (1 - U_i)) / (n_i_dot - 1)
                + P_j ** 2 * sum(
                    (error_matrix[i, j] / n_i_dot) * (1 - (error_matrix[i, j] / n_i_dot)) / (n_i_dot - 1)
                    for i in range(q) if i != j
                )
            )
            se_P_j = np.sqrt(variance_P_j)
            producer_accuracy_se.append(se_P_j)
            
            ci_lower = P_j - z_score * se_P_j
            ci_upper = P_j + z_score * se_P_j
            producer_accuracy_ci_lower.append(ci_lower)
            producer_accuracy_ci_upper.append(ci_upper)
            ci_value =  z_score * se_P_j
            producer_accuracy_ci_value.append(ci_value)
        else:
            producer_accuracy_se.append(np.nan)
            producer_accuracy_ci_lower.append(np.nan)
            producer_accuracy_ci_upper.append(np.nan)
            producer_accuracy_ci_value.append(np.nan)

    # Overall accuracy calculations
    overall_accuracy = np.trace(area_proportion_matrix)
    overall_accuracy_variance = sum(
        weights[i] ** 2 * user_accuracy[i] * (1 - user_accuracy[i]) / (error_matrix[i, :].sum() - 1)
        for i in range(q) if error_matrix[i, :].sum() > 1
    )
    overall_accuracy_se = np.sqrt(overall_accuracy_variance)
    
    # Confidence interval for overall accuracy
    overall_accuracy_ci_lower = overall_accuracy - z_score * overall_accuracy_se
    overall_accuracy_ci_upper = overall_accuracy + z_score * overall_accuracy_se
    overall_accuracy_ci_value =  z_score * overall_accuracy_se

    return {
        "user_accuracy": [float(ua) for ua in user_accuracy],
        "user_accuracy_se": user_accuracy_se,
        "user_accuracy_ci_lower": user_accuracy_ci_lower,
        "user_accuracy_ci_upper": user_accuracy_ci_upper,
        "user_accuracy_ci_value": user_accuracy_ci_value,
        "producer_accuracy": [float(pa) for pa in producer_accuracy],
        "producer_accuracy_se": producer_accuracy_se,
        "producer_accuracy_ci_lower": producer_accuracy_ci_lower,
        "producer_accuracy_ci_upper": producer_accuracy_ci_upper,
        "producer_accuracy_ci_value": producer_accuracy_ci_value,
        "overall_accuracy": float(overall_accuracy),
        "overall_accuracy_se": overall_accuracy_se,
        "overall_accuracy_ci_lower": overall_accuracy_ci_lower,
        "overall_accuracy_ci_upper": overall_accuracy_ci_upper,
        "overall_accuracy_ci_value": overall_accuracy_ci_value
    }


def calculate_error_adjusted_area(area_proportion_matrix, total_area):
    """
    Calculate error-adjusted area estimates for each category without including units.
    """
    adjusted_areas = total_area * area_proportion_matrix.sum(axis=0)
    adjusted_areas_no_units = [f"{area:.4f}" for area in adjusted_areas]
    return adjusted_areas_no_units


def calculate_standard_error_and_ci(error_matrix, weights, area_adjusted_estimates, total_area, confidence_level=0.95):
    """
    Calculate the standard error and confidence interval for error-adjusted area estimates.
    
    Parameters:
    - error_matrix: 2D numpy array, original error matrix with sample counts.
    - weights: List of mapped area proportions for each class.
    - area_adjusted_estimates: List of error-adjusted area estimates for each class.
    - total_area: Total area of the map.
    - confidence_level: Confidence level for the interval (default is 0.95).
    
    Returns:
    - A dictionary containing standard errors and confidence intervals for each class.
    """
    q = error_matrix.shape[0]
    standard_errors = []

    # Calculate standard errors for each class using original error matrix
    for j in range(q):
        se_j_squared = 0
        for i in range(q):
            n_ij = error_matrix[i, j]
            n_i_dot = error_matrix[i, :].sum()
            if n_i_dot > 1:
                se_j_squared += float(weights[i]) ** 2 * (n_ij / n_i_dot) * (1 - n_ij / n_i_dot) / (n_i_dot - 1)
        standard_errors.append(np.sqrt(se_j_squared))

    # Convert to area standard errors by scaling with total_area
    area_standard_errors = [total_area * se_j for se_j in standard_errors]
    
    # Confidence interval multiplier
    z_score = st.norm.ppf(1 - (1 - confidence_level) / 2)
    
    # Calculate confidence intervals based on error-adjusted area estimates
    confidence_intervals = [
        (
            area_adjusted_estimates[j] - z_score * area_standard_errors[j],
            area_adjusted_estimates[j] + z_score * area_standard_errors[j]
        )
        for j in range(q)
    ]

    # Format results without units
    standard_errors_formatted = [f"{se:.4f}" for se in area_standard_errors]
    confidence_intervals_formatted = [(f"{ci[0]:.4f}", f"{ci[1]:.4f}") for ci in confidence_intervals]

    return {
        "standard_errors": standard_errors_formatted,
        "confidence_intervals": confidence_intervals_formatted
    }



# Input the error matrix values here and mapped number of pixels per class
error_matrix = np.array([[354, 46], [33, 67]], dtype=np.float64) # Error matrix of sample counts.
mapped_pixels = np.array([42640, 4005], dtype=np.float64)  # Number of pixels mapped to each class
pixel_size = 30  # Pixel size --> ground length of one side of the pixel; specify units below for the equivalent ground area of a pixel)
pixel_area = (pixel_size ** 2)
unit = "m²" # Unit of area; change accordingly

# Calculate total area based on pixel size and mapped pixels
total_area = np.sum(mapped_pixels) * pixel_area

# Calculate weights based on mapped pixels
weights = calculate_weights(mapped_pixels)

# Step 1: Display the original error matrix with totals, mapped pixels, and weights (with W_i total added)
mapped_pixels_total = np.sum(mapped_pixels)
mapped_pixels_column = np.append(mapped_pixels, mapped_pixels_total)

weights_total = np.sum(weights)
weights_column = np.append(weights, weights_total)

error_matrix_with_totals = np.hstack([error_matrix, error_matrix.sum(axis=1, keepdims=True)])
totals_row = np.append(error_matrix.sum(axis=0), error_matrix.sum())

error_matrix_df = pd.DataFrame(
    np.vstack([error_matrix_with_totals, totals_row]),
    columns=["Reference_1", "Reference_0", "Total"],
    index=["Mapped_1", "Mapped_0", "Total"]
)
error_matrix_df["Mapped_Pixels"] = mapped_pixels_column
error_matrix_df["W_i"] = weights_column

# Step 2: Convert error matrix to area proportions
area_proportion_matrix = convert_to_area_proportion(error_matrix, weights)

# Step 3: Add row and column totals to the area proportion matrix
area_proportion_with_totals = np.hstack([area_proportion_matrix, area_proportion_matrix.sum(axis=1, keepdims=True)])
area_totals_row = np.append(area_proportion_matrix.sum(axis=0), area_proportion_matrix.sum())

area_proportion_df = pd.DataFrame(
    np.vstack([area_proportion_with_totals, area_totals_row]),
    columns=["Reference_1", "Reference_0", "Total"],
    index=["Mapped_1", "Mapped_0", "Total"]
)

# Step 4: Calculate overall accuracy as part of the accuracy metrics and include it in the output
accuracy_metrics = calculate_accuracy_metrics(area_proportion_matrix, error_matrix, weights)

# Step 5: Calculate error-adjusted area estimates with units
error_adjusted_area = calculate_error_adjusted_area(area_proportion_matrix, total_area)

# Step 6: Calculate standard errors and confidence intervals, using error-adjusted area estimates
# Ensure 'error_adjusted_area' is converted to numeric values before passing
area_adjusted_estimates = [float(area.split(" ")[0]) for area in error_adjusted_area]
standard_error_and_ci = calculate_standard_error_and_ci(error_matrix, weights, area_adjusted_estimates, total_area)

# Save Results as CSV
error_matrix_df.to_csv("error_matrix_with_totals.csv", index=True)
area_proportion_df.to_csv("area_proportion_matrix_with_totals.csv", index=True)

# Separate confidence intervals into lower and upper bounds without units
confidence_interval_lower = [ci[0].split(" ")[0] for ci in standard_error_and_ci["confidence_intervals"]]
confidence_interval_upper = [ci[1].split(" ")[0] for ci in standard_error_and_ci["confidence_intervals"]]

# Define confidence level
confidence_level = 0.95  # or any other level you prefer

# Calculate the CI values for area estimates by removing units and converting to float
area_standard_errors_numeric = [float(se.split(" ")[0]) for se in standard_error_and_ci["standard_errors"]]
z_score = st.norm.ppf(1 - (1 - confidence_level) / 2)
area_ci_values = [f"{z_score * se:.4f}" for se in area_standard_errors_numeric]  # CI Value for Area Estimates

# Save calculated metrics as a summary CSV file with units specified in the column names
metrics_summary = pd.DataFrame({
    "Class": ["Class_1", "Class_0", "Overall Accuracy"],
    
    # User Accuracy
    "User_Accuracy": accuracy_metrics["user_accuracy"] + [None],
    "User_Accuracy_SE": accuracy_metrics["user_accuracy_se"] + [None],
    "User_Accuracy_CI_Value": accuracy_metrics["user_accuracy_ci_value"] + [None],  # CI Value for User Accuracy
    
    # Producer Accuracy
    "Producer_Accuracy": accuracy_metrics["producer_accuracy"] + [None],
    "Producer_Accuracy_SE": accuracy_metrics["producer_accuracy_se"] + [None],
    "Producer_Accuracy_CI_Value": accuracy_metrics["producer_accuracy_ci_value"] + [None],  # CI Value for Producer Accuracy
    
    # Overall Accuracy
    "Overall_Accuracy": [None, None, accuracy_metrics["overall_accuracy"]],
    "Overall_Accuracy_SE": [None, None, accuracy_metrics["overall_accuracy_se"]],
    "Overall_Accuracy_CI_Value": [None, None, accuracy_metrics["overall_accuracy_ci_value"]],  # CI Value for Overall Accuracy
    
    # Error-Adjusted Area Estimates
    "Error_Adjusted_Area [m²]": area_adjusted_estimates + [None],
    
    # Standard Errors and Confidence Intervals for Area Estimates
    "Area_Standard_Error [m²]": standard_error_and_ci["standard_errors"] + [None],
    "Area_CI_Value [m²]": area_ci_values + [None]  # CI Value for Area Estimates
})

# Save to CSV with UTF-8 encoding and no extraneous characters
metrics_summary.to_csv("calculated_metrics_summary.csv", index=False, encoding="utf-8-sig")

# Display a message indicating the file has been saved successfully
print("Metrics summary has been saved to 'calculated_metrics_summary.csv'.")


# Print the Original Error Matrix with Totals, Mapped Pixels, and Weights
print("Original Error Matrix with Totals, Mapped Pixels, and Weights:")
print(error_matrix_df.to_string(index=True))
print("\n")

# Print the Area Proportion Matrix with Totals
print("Area Proportion Matrix with Totals:")
print(area_proportion_df.to_string(index=True))
print("\n")

# Print the Accuracy Metrics
print("Accuracy Metrics:")
print(f"User's Accuracy: {accuracy_metrics['user_accuracy']}")
print(f"User's Accuracy SE: {np.round(accuracy_metrics['user_accuracy_se'], 4)}")
print(f"User's Accuracy CI Lower: {np.round(accuracy_metrics['user_accuracy_ci_lower'], 4)}")
print(f"User's Accuracy CI Upper: {np.round(accuracy_metrics['user_accuracy_ci_upper'], 4)}")
print(f"User's Accuracy CI Value: {np.round(accuracy_metrics['user_accuracy_ci_value'], 4)}\n")

print(f"Producer's Accuracy: {accuracy_metrics['producer_accuracy']}")
print(f"Producer's Accuracy SE: {np.round(accuracy_metrics['producer_accuracy_se'], 4)}")
print(f"Producer's Accuracy CI Lower: {np.round(accuracy_metrics['producer_accuracy_ci_lower'], 4)}")
print(f"Producer's Accuracy CI Upper: {np.round(accuracy_metrics['producer_accuracy_ci_upper'], 4)}")
print(f"Producer's Accuracy CI Value: {np.round(accuracy_metrics['producer_accuracy_ci_value'], 4)}\n")

print(f"Overall Accuracy: {accuracy_metrics['overall_accuracy']:.4f}")
print(f"Overall Accuracy SE: {accuracy_metrics['overall_accuracy_se']:.4f}")
print(f"Overall Accuracy CI Lower: {accuracy_metrics['overall_accuracy_ci_lower']:.4f}")
print(f"Overall Accuracy CI Upper: {accuracy_metrics['overall_accuracy_ci_upper']:.4f}")
print(f"Overall Accuracy CI Value: {accuracy_metrics['overall_accuracy_ci_value']:.4f}\n")

# Print Error-Adjusted Area Estimates
print("Error-Adjusted Area Estimates:", error_adjusted_area)

# Print Standard Errors and Confidence Intervals for Area Estimates
print("Standard Errors:", [f"{se}" for se in standard_error_and_ci["standard_errors"]])
print("Confidence Intervals Lower Bound:", [f"{ci}" for ci in confidence_interval_lower])
print("Confidence Intervals Upper Bound:", [f"{ci}" for ci in confidence_interval_upper])
