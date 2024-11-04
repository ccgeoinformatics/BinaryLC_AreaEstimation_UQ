# BinaryLC_AreaEstimation_UQ
This repository provides Python scripts for performing error analysis, accuracy assessment, unbiased area estimation, and uncertainty quantification for binary land cover classification maps. The methodology is based on established practices outlined in Olofsson et al. (2013, 2014).
Scripts
1. single_binary_lc_accuracy_areaEstimation_uq.py
This script processes a single binary land cover map using an error matrix and mapped pixel counts as inputs. It calculates key accuracy metrics, including user’s accuracy, producer’s accuracy, and overall accuracy. Additionally, it performs error-adjusted area estimation and calculates standard errors and confidence intervals for the accuracy metrics and area estimates.

Key Calculations:

User's, producer's, and overall accuracy, with standard errors and confidence intervals.
Error-adjusted area estimates.
Standard errors and confidence intervals for area estimates.
References:

Olofsson, P., et al. (2013). Making better use of accuracy data in land change studies: Estimating area and change accuracy. Remote Sensing of Environment, 129, 122-131.
Olofsson, P., et al. (2014). Good practices for estimating area and assessing accuracy of land change. Remote Sensing of Environment, 148, 42-57.
2. multiple_binary_lc_accuracy_areaEstimation_uq.py
This script allows processing of multiple binary land cover maps in an iterative fashion. It reads data from a CSV file where each row represents a different land cover map, including columns for error matrix values, mapped pixel counts, pixel size, and unit. The script calculates accuracy metrics, error-adjusted area estimates, and uncertainty measures for each map.

CSV Input Structure:

lcmap, n11, n12, n21, n22, mapped_pixels_class_1, mapped_pixels_class_0, pixel_size, unit
Example row: GISD30_2000,354,46,33,67,42640,4005,30,m2
Key Calculations:

For each map in the CSV file:
User's, producer's, and overall accuracy, with standard errors and confidence intervals.
Error-adjusted area estimates, with standard errors and confidence intervals.
Usage
Requirements
The scripts require Python 3 and the following libraries:

numpy
pandas
scipy
Running the Scripts
Clone this repository and navigate to the folder containing the scripts.

Prepare the necessary input data:

For single_binary_lc_accuracy_areaEstimation_uq.py, define the error matrix, mapped pixels, pixel size, and unit in the script.
For multiple_binary_lc_accuracy_areaEstimation_uq.py, provide a CSV file (input.csv by default) with the structure specified above.
Run the script with the command:

bash
Copy code
python single_binary_lc_accuracy_areaEstimation_uq.py
or

bash
Copy code
python multiple_binary_lc_accuracy_areaEstimation_uq.py
Outputs
Each script will produce:

CSV files containing the results of error matrix analysis, area proportion matrix, accuracy metrics, error-adjusted area estimates, and confidence intervals.
Printed output summarizing the results.
Example
After running multiple_binary_lc_accuracy_areaEstimation_uq.py, results for each map will be saved as {lcmap}_results.csv, where {lcmap} is the name of each land cover map specified in the CSV file.

License
This project is licensed under the MIT License. See the LICENSE file for details.

Contact
For questions or further assistance, contact Jojene R. Santillan at:

Email: santillan@ipi.uni-hannover.de, jrsantillan@carsu.edu.ph
Institute of Photogrammetry and GeoInformation (IPI), Leibniz University Hannover, Germany
Caraga State University, Butuan City, Philippines
