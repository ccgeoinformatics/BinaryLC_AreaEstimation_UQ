
# Binary/Two-class Land Cover Map Accuracy Assessment, Area Estimation, and Uncertainty Quantification

This repository provides tools for assessing the accuracy of binary or two-class land cover classification maps, including unbiased class area estimation and uncertainty quantifications following the methodologies presented in Olofsson et al. (2013, 2014). It includes scripts for single and multiple binary land cover maps, allowing for iterative accuracy assessment and uncertainty quantification.

## Features

- **Accuracy Assessment**: Calculate user's, producer's, and overall accuracy metrics.
- **Area Estimation**: Provide unbiased area estimates with associated uncertainty metrics.
- **Uncertainty Quantification**: Calculate confidence intervals and standard errors for accuracy and area estimates.
- **Iterative Analysis**: Support for multiple land cover maps in one iteration.

## Methodology

The methodology follows these references:
1. **Olofsson et al. (2013)**: Making better use of accuracy data in land change studies: Estimating area and change accuracy. *Remote Sensing of Environment*, 129, 122-131.
   [https://doi.org/10.1016/j.rse.2012.10.031](https://doi.org/10.1016/j.rse.2012.10.031)
2. **Olofsson et al. (2014)**: Good practices for estimating area and assessing accuracy of land change. *Remote Sensing of Environment*, 148, 42-57.
   [https://doi.org/10.1016/j.rse.2014.02.015](https://doi.org/10.1016/j.rse.2014.02.015)

## Scripts

1. **single_binary_lc_accuracy_areaEstimation_uq.py**: Processes a single binary land cover map, performing accuracy assessment, area estimation, and uncertainty quantification.
2. **multiple_binary_lc_accuracy_areaEstimation_uq.py**: Processes multiple binary land cover maps from a CSV file.

### CSV File Structure for Multiple Land Cover Map Script

The input CSV file for the multiple land cover map script should have the following columns:

- `lcmap`: Name of the land cover map
- `n11, n12, n21, n22`: Error matrix values
- `mapped_pixels_class_1, mapped_pixels_class_0`: Mapped pixel counts for each class
- `pixel_size`: Pixel size in ground units
- `unit`: Unit of area (e.g., m²)

Example:
```csv
lcmap,n11,n12,n21,n22,mapped_pixels_class_1,mapped_pixels_class_0,pixel_size,unit
GISD30_2000,354,46,33,67,42640,4005,30,m²
GISD30_2005,363,37,33,67,43091,3554,30,m²
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/ccgeoinformatics/BinaryLC_AreaEstimation_UQ.git
   ```
2. Install required Python packages:
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
