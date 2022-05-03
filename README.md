
# stJointcount: Quantifying spatial dependencies between clusters for spatial transcriptomic data

Spatial dependency helps determine whether there is a match between locational and attribute similarity. The measure looks at whether an attribute of a variable observed at one location is independent of values observed at neighboring locations. If these values are similar and statistically significant, positive spatial autocorrelation exists in the spatial distribution. In the case of neighboring locations exhibiting differences (dissimilarities), there is little or no spatial autocorrelation in the spatial distribution. Join count analysis is a quantitative method that determines the degree of clustering or dispersion among a set of spatially adjacent polygons. 

**stJoincount**, the application of join count analysis to the spatial transcriptomics dataset. This tool converts the spatial map into a raster object (a two-dimensional image as a rectangular matrix or grid of square pixels), where clusters labelled spots are converted to adjacent pixels with a calculated resolution. A neighbors' list was created based on the rasterized sample, which identifies adjacent and diagonal neighbors for each pixel. After adding binary spatial weights to the neighbors' list, a multi-categorical join count analysis is then performed, allowing all possible combinations of cluster pairings to be tabulated. The function returns the observed join counts, the expected count under conditions of spatial randomness, and the variance of observed to expected calculated under non-free sampling. The z-score is then calculated as the difference between observed and expected counts, divided by the square root of the variance. 

## Package User Guide

Detail introduction of each function in this package can be found in `/vignettes/stJointcount-vignette.Rmd`

Following codes show main functions of rasterization and joint count analysis. In this demo, we are going to use an human breast cancer sample that could be downloaded from 10x Genomics website [here](https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0). User can download the related folders and files directly from the website, or they can follow first two blocks of codes in the vignette, saving them in the R temp folder.

### Data input

The input is a directory containing the H5 file, the image data (in a subdirectory ‘spatial’), and analysis results (in a subdirectory ‘analysis’).

```r
sample <- spatialDataPrep(outDir)
```

Note: Users can skip this step if they want to use the data with customized clusters attached. The data format required by this package is the column name of customized cluster labels must be named as "Cluster" (the current version only accepts cluster labels as integer "1,2,3..."). It can be checked by `head(sample@meta.data)` and can be changed as needed.

### Resolution calculation

```r
resolutionList <- resolutionCalc(sample)
```

### Rasterization

```r
mosaicIntegration <- rasterizeEachCluster(sample)
mosaicIntPlot(mosaicIntegration)
```

<p align="center"><img src="https://github.com/Nina-Song/stJoincount/blob/master/inst/extdata/rasterization.png" height="400"></p>

### Joint count analysis

```r
joincount.result <- joincountAnalysis(mosaicIntegration)
```

<p align="center"><img src="https://github.com/Nina-Song/stJoincount/blob/master/inst/extdata/joincountResult.png" height="400"></p>

### Z-score Heatmap

```r
heatmapMatrix <- zscoreMatrix(sample, joincount.result)

zscorePlot(heatmapMatrix)
```

<p align="center"><img src="https://github.com/Nina-Song/stJoincount/blob/master/inst/extdata/zscoreHeatmap.png" height="400"></p>


### Installation
```r
library(devtools)
install_github("Nina-Song/stJoincount")
```  

### Contact
songjiar@usc.edu
