
# Binary/Two-class Land Cover Map Accuracy Assessment, Area Estimation, and Uncertainty Quantification

This repository provides tools for assessing the accuracy of binary or two-class land cover classification maps, including unbiased class area estimation and uncertainty quantifications following the methodologies presented in Olofsson et al. (2013, 2014). It includes scripts for single and multiple binary land cover maps, allowing for iterative accuracy assessment and uncertainty quantification.

## Features

- **Accuracy Assessment**: Calculate user's, producer's, and overall accuracy metrics.
- **Area Estimation**: Provide unbiased area estimates with associated uncertainty metrics.
- **Uncertainty Quantification**: Calculate confidence intervals and standard errors for accuracy and area estimates.
- **Iterative Analysis**: Support for multiple land cover maps in one iteration.

## Methodology

The methodology follows these references:
1. **Olofsson et al. (2013)**: Making better use of accuracy data in land change studies: Estimating accuracy and area and quantifying uncertainty using stratified estimation. *Remote Sensing of Environment*, 129, 122-131.
   [https://doi.org/10.1016/j.rse.2012.10.031](https://doi.org/10.1016/j.rse.2012.10.031)
2. **Olofsson et al. (2014)**: Good practices for estimating area and assessing accuracy of land change. *Remote Sensing of Environment*, 148, 42-57.
   [https://doi.org/10.1016/j.rse.2014.02.015](https://doi.org/10.1016/j.rse.2014.02.015)
   
Users of the scripts are advised to read the above references for the theoretical concepts behind the methods, including the assumptions (for example, the error matrix shall be generated from a stratified random sampling design).

## Scripts

1. **single_binary_lc_accuracy_areaEstimation_uq.py**: Processes a single binary land cover map, performing accuracy assessment, area estimation, and uncertainty quantification. The entries of the error (confusion) matrix shall be manually entered or provided within the script, including the number of pixels mapped/classified for each class.
2. **multiple_binary_lc_accuracy_areaEstimation_uq.py**: Processes multiple binary land cover maps from a CSV file.

### CSV File Structure for Multiple Land Cover Map Script

The input CSV file for the multiple land cover map script should have the following columns:

- `lcmap`: Name of the land cover map
- `n11, n12, n21, n22`: Error matrix values
- `mapped_pixels_class_1, mapped_pixels_class_0`: Mapped pixel counts for each class
- `pixel_size`: Pixel size in ground units
- `unit`: Unit of area (e.g., m²)

Example contents of input.csv file:
```csv
lcmap,n11,n12,n21,n22,mapped_pixels_class_1,mapped_pixels_class_0,pixel_size,unit
GISD30_2000,354,46,33,67,42640,4005,30,m²
GISD30_2005,363,37,33,67,43091,3554,30,m²
```

## Two-Class Error Matrix

Both scripts assume that the error/confusion matrix is structured as follows ('1' and '0' are the two classes; but this can be modified in the scripts, e.g, '1' and '2', or 'Builtup' and'Non-Builtup', etc.)

|               | Reference Class 1 | Reference Class 0 | Total  |
|---------------|-------------------|-------------------|--------|
| **Mapped 1**  | n11               | n12               | n1.    |
| **Mapped 0**  | n21               | n22               | n2.    |
| **Total**     | n.1               | n.2               | N      |

- **n11**: The number of pixels correctly classified as **Class 1** (true positives for Class 1).
- **n12**: The number of pixels misclassified as **Class 1** but are actually **Class 0** (false positives for Class 1).
- **n21**: The number of pixels misclassified as **Class 0** but are actually **Class 1** (false negatives for Class 1).
- **n22**: The number of pixels correctly classified as **Class 0** (true negatives for Class 1).

### Row Totals:
- **n1.**: Total pixels mapped as **Class 1** (sum of n11 and n12).
- **n2.**: Total pixels mapped as **Class 0** (sum of n21 and n22).

### Column Totals:
- **n.1**: Total reference pixels for **Class 1** (sum of n11 and n21).
- **n.2**: Total reference pixels for **Class 0** (sum of n12 and n22).

### Overall Total:
- **N**: Total number of pixels in the error matrix (sum of all elements).



## Installation

### Option 1: Clone the Repository

Download the repository using Git:

```bash
git clone https://github.com/ccgeoinformatics/BinaryLC_AreaEstimation_UQ.git
cd BinaryLC_AreaEstimation_UQ
```
### Option 2: Download as a ZIP File
1. Go to the [repository on GitHub](https://github.com/ccgeoinformatics/BinaryLC_AreaEstimation_UQ).
2. Click on the "Code" button.
3. Select "Download ZIP" and extract the downloaded file.

### Install Required Packages
After obtaining the files via either option, install the dependencies using the requirements.txt file:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Run the single land cover map script:
   ```bash
   python single_binary_lc_accuracy_areaEstimation_uq.py
   ```
2. Run the multiple land cover map script:
   ```bash
   python multiple_binary_lc_accuracy_areaEstimation_uq.py
   ```

## Contact

For any questions or issues, please reach out to:

**Jojene R. Santillan**  
Institute of Photogrammetry and GeoInformation (IPI), Leibniz University Hannover, Germany  
& Caraga Center for Geo-Informatics & Department of Geodetic Engineering, College of Engineering and Geosciences, Caraga State University, Butuan City, Philippines  
santillan@ipi.uni-hannover.de, jrsantillan@carsu.edu.ph
