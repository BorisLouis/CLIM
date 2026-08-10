# Correlation Clustering Imaging (CLIM)

Correlation Clustering Imaging (CLIM) is an image-analysis method that uses spatiotemporal correlations in intensity fluctuations to identify functional domains in luminescent materials.

CLIM analyzes time-dependent wide-field fluorescence microscopy images and groups pixels displaying correlated temporal intensity dynamics.

The method was introduced and validated in:

> **B. Louis, S. Seth, Q. An, R. Ji, Y. Vaynzof, J. Hofkens, and I. G. Scheblykin**  
> *In Operando Locally-Resolved Photophysics in Perovskite Solar Cells by Correlation Clustering Imaging*  
> **Advanced Materials** 37, 2413126 (2025)  
> https://doi.org/10.1002/adma.202413126

The CLIM algorithm was developed by Boris Louis, with additional input from Sudipta Seth.

---

## License and permitted use

CLIM is provided for **non-commercial academic research and teaching only**.

Commercial use is **not permitted** under the academic license. This includes, among other things, use of CLIM by or on behalf of a commercial company, integration into commercial software or services, and commercial research and development.

Commercial organizations interested in using CLIM should contact:

**Boris Louis**  
boris.louis@kuleuven.be

for information regarding separate permission or commercial licensing.

If CLIM contributes to scientific work, publications, presentations, reports, or other publicly disseminated research, please cite the original CLIM paper:

> B. Louis et al., *Advanced Materials* 37, 2413126 (2025).  
> https://doi.org/10.1002/adma.202413126

See the `LICENSE` file for the complete terms of use.

> **Note:** The licensing wording in this repository should be reviewed and approved by the applicable institutional rights holder(s), particularly KU Leuven and/or Lund University.

---

## Disclaimer

This software is research code provided to support reproduction and further academic use of the CLIM method described in the publication above.

The software is provided **as is**, without warranty that it will function for a particular dataset, material system, microscope, operating system, or application.

For questions regarding the code or potential applications of CLIM, contact:

**boris.louis@kuleuven.be**

---

# Using CLIM

The main script for running CLIM is:

```matlab
mainCorrMovie.m
```

CLIM expects a time-dependent image sequence. The current implementation supports the following input formats:

- `.mat`
- `.tif`
- `.ome.tif`
- `.spe`

If your data are stored in another format, they can for example be converted to TIFF using Fiji/ImageJ.

---

## Input data

The provided test dataset is a reduced version of the data used to demonstrate CLIM. It contains a smaller field of view and fewer frames to allow the analysis pipeline to be tested quickly.

`file.path` must point to a folder containing the image sequence to be analyzed, and `file.ext` specifies its extension.

Example:

```matlab
%% User input

file.path = 'path\to\your\data';
file.ext  = '.tif';
```

---

# Main parameters

## Run mode

```matlab
info.runMethod = 'load';
```

Two options are available:

- `'load'` — CLIM attempts to reuse results from a previous analysis when compatible processed data are found.
- `'run'` — forces the analysis to be recalculated even when previous results are present.

---

## Drift correction

```matlab
info.driftCorr = true;
```

- `true` — perform drift correction.
- `false` — do not perform drift correction.

Drift correction is recommended when spatial drift is present during the recording, as CLIM relies on temporal correlations between spatially fixed pixels.

---

## Global correlated signal correction

```matlab
deconvolve = true;
```

When enabled, CLIM attempts to remove intensity dynamics that are correlated over the entire image. This can help separate local fluctuations from global changes in excitation, photobleaching, or other common-mode intensity variations.

---

## Background removal

```matlab
backgroundThresh = 0.1;
```

This parameter is intended to control segmentation of signal-containing pixels from the background.

> **Implementation note:** In the current CLIM v1.0 implementation, the background threshold inside `CorrClusterMovie.deconvolve()` is fixed to `0.1`. The user-defined `backgroundThresh` value will only control this step after the implementation is updated to use the passed parameter.

---

# Correlation threshold

CLIM currently provides three clustering modes:

```matlab
info.thresholdMode = 'auto';   % 'auto', 'fixed', or 'None'
```

## Automatic threshold

```matlab
info.thresholdMode = 'auto';

minCorr  = 0.4;
stepCorr = 0.05;
maxCorr  = 0.9;
```

In automatic mode, CLIM analyzes a smaller region of the image using a range of correlation thresholds:

```text
minCorr : stepCorr : maxCorr
```

For each threshold, a clustering solution is generated and evaluated.

The current implementation selects the threshold using a metric combining:

1. the mean silhouette value of the clustering; and
2. the fraction of the image assigned to clusters.

This reduces user bias in threshold selection while penalizing thresholds that produce apparently good clusters for only a very small fraction of the image.

---

## Fixed threshold

```matlab
info.thresholdMode = 'fixed';
threshold = 0.5;
```

In this mode, CLIM uses the specified correlation threshold directly.

This can be useful when comparing datasets for which the same clustering criterion should be maintained.

---

## Alternative clustering mode

```matlab
info.thresholdMode = 'None';
```

This activates the alternative `corrClusteringNoThresh` implementation.

This mode resolves competing cluster assignments by comparing pixels with candidate clusters instead of relying exclusively on the sequential threshold-based growth used by the standard clustering algorithm.

Note that the current implementation still uses `minCorr` as a minimum local correlation criterion.

This mode should presently be considered experimental.

---

# Region used for automatic threshold optimization

```matlab
testROIRadius = 64;
```

In automatic threshold mode, CLIM initially evaluates thresholds within a smaller region centered in the dataset.

`testROIRadius` determines the radius, in pixels, of this region.

---

# Analyze a specific region of interest

```matlab
info.ROI = false;

% ROI format:
% [x y width height]

ROI = [];

% Example:
% ROI = [5 71 230 120];
```

Set:

```matlab
info.ROI = true;
```

to restrict the complete CLIM analysis to the specified region.

---

# Frames used for correlation analysis

```matlab
frame2Process = 1:2000;
```

This variable specifies which movie frames are used to calculate temporal correlations.

The optimal number of frames depends on the temporal dynamics of the system, signal-to-noise ratio, acquisition rate, and total recording duration.

In the original CLIM paper, the experimental PL movies consisted of 6000 frames recorded with a 50 ms exposure time. Users should ensure that the selected frames sufficiently sample the relevant temporal dynamics of their system.

---

# Cluster intensity extraction

```matlab
method = 'Mean';
```

The default method extracts the mean time-dependent intensity of the pixels assigned to each CLIM cluster.

These traces can subsequently be used to investigate the temporal dynamics of individual functional domains.

---

# CLIM outputs

CLIM converts a time-dependent fluorescence movie into several spatial maps describing the local temporal dynamics of the sample.

## 1. Correlation map

The **correlation map** displays the average temporal correlation of each pixel with its directly neighboring pixels.

High values indicate areas in which neighboring pixels undergo similar temporal intensity dynamics.

Low-correlation regions can reveal boundaries between independently behaving functional domains even when these boundaries are not visible in a time-averaged fluorescence image.

---

## 2. Cluster map

The **cluster map** groups pixels into spatial regions displaying highly correlated temporal behavior.

Each integer label corresponds to one CLIM cluster. When displayed as an RGB image, colors are used only to distinguish different clusters; identical or similar colors do not imply similar photophysical behavior.

The clusters should generally be interpreted as **functional domains** rather than automatically as structural domains.

In the MAPI thin films investigated in the original CLIM publication, the clusters closely corresponded to individual morphological grains. This correspondence should not be assumed for other material systems.

---

## 3. Internal correlation map

The **internal correlation map** evaluates the average correlation of each pixel with other pixels in the cluster to which it has been assigned.

It therefore provides a measure of the internal consistency and quality of the clustering.

---

## 4. Silhouette map

The **silhouette map** evaluates how much better a pixel correlates with its assigned cluster compared with nearby alternative clusters.

For a pixel, CLIM compares:

- its average correlation with pixels belonging to the same cluster (**internal correlation**); and
- its strongest correlation with pixels belonging to nearby clusters (**external correlation**).

The silhouette score is calculated from the difference between these values divided by the larger of them.

A high silhouette value indicates that the pixel is strongly associated with its own cluster and comparatively weakly correlated with neighboring clusters.

A low silhouette value indicates that the pixel also displays dynamics similar to another cluster.

The silhouette map therefore provides both:

- a measure of clustering quality; and
- information about correlations or possible communication between neighboring functional domains.

In the perovskite systems studied in the original publication, low silhouette values between grains were interpreted as evidence of correlated dynamics between otherwise distinct domains.

---

## 5. Cleaned cluster mask

CLIM additionally generates a cleaned cluster mask in which pixels with low silhouette values are removed.

In the current implementation:

```matlab
silhouette < 0.2
```

is excluded from the cleaned mask.

---

## 6. Cluster intensity traces

For every identified cluster, CLIM can extract its time-dependent intensity trace.

These traces enable further analysis of the dynamics responsible for the correlated behavior detected by CLIM.

---

# Saved results

The final structure generated by:

```matlab
corrOutput = myMovie.generateResults;
```

contains the main CLIM results, including:

```text
corrOutput.Image
corrOutput.corrMap
corrOutput.corrMask
corrOutput.cleanMask
corrOutput.silMap
corrOutput.corrClustMap
corrOutput.threshold
corrOutput.results
```

### `corrOutput.corrMap`

Pixel-to-neighbor correlation map.

### `corrOutput.corrMask`

Raw CLIM cluster-label map.

### `corrOutput.cleanMask`

Cluster map after removing pixels with low silhouette values.

### `corrOutput.silMap`

Pixel-resolved silhouette map.

### `corrOutput.corrClustMap`

Map describing the correlation of pixels within their assigned clusters.

### `corrOutput.threshold`

Correlation threshold used for the analysis.

### `corrOutput.results`

Per-cluster quantitative information, including metrics such as:

- number of pixels;
- mean intra-cluster correlation;
- minimum correlation;
- correlation variability;
- mean intensity;
- fluctuation amplitude;
- mean inter-cluster correlation;
- maximum inter-cluster correlation; and
- mean silhouette value.

---

# Scientific interpretation

CLIM identifies regions exhibiting correlated temporal intensity dynamics.

The original publication demonstrated CLIM using photoluminescence movies of MAPbI3 perovskite thin films and operating solar cells.

In perovskite thin films, CLIM clusters closely reproduced the morphological grain structure observed by scanning electron microscopy.

In functioning solar cells, however, correlated CLIM domains could extend across several grains, demonstrating that CLIM is sensitive to functional dynamics rather than morphology alone.

CLIM can therefore potentially be applied to other luminescent materials and devices provided that measurable local temporal intensity dynamics are present.

For detailed information concerning the CLIM algorithm, validation, resolution, noise sensitivity, and physical interpretation, please consult the original publication.

---

# Citation

If CLIM contributes to your work, please cite:

B. Louis, S. Seth, Q. An, R. Ji, Y. Vaynzof, J. Hofkens, and I. G. Scheblykin,  
**In Operando Locally-Resolved Photophysics in Perovskite Solar Cells by Correlation Clustering Imaging**,  
*Advanced Materials* **37**, 2413126 (2025).  
https://doi.org/10.1002/adma.202413126

---

# Contact

For scientific questions, technical support, collaborations, or licensing:

**Boris Louis**  
KU Leuven  
boris.louis@kuleuven.be
