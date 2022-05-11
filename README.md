
# stJointcount: Quantifying spatial dependencies in labeled and mapped spatial transcriptomic data

Spatial dependency is the relationship between locational and attribute similarity. The measure reflects whether an attribute of a variable observed at one location is independent of values observed at neighboring locations. Positive spatial dependency exists when neighboring attributes are more similar than what could be explained by chance. Likewise, a negative spatial dependency is reflected by a dissimilarity of neighboring attributes. Join count analysis allows for quantification of the spatial dependencies of nominal data in an arrangement of spatially adjacent polygons. 

**stJoincount** facilitates the application of join count analysis to spatial transcriptomic data generated from the 10x Genomics Visium platform. This tool first converts a labeled spatial tissue map into a raster object, in which each spatial feature is represented by a pixel coded by label assignment. This process includes automatic calculation of optimal raster resolution and extent for the sample. A neighbors list is then created from the rasterized sample, in which adjacent and diagonal neighbors for each pixel are identified. After adding binary spatial weights to the neighbors list, a multi-categorical join count analysis is performed to tabulate "joins" between all possible combinations of label pairs. The function returns the observed join counts, the expected count under conditions of spatial randomness, and the variance calculated under non-free sampling. The z-score is then calculated as the difference between observed and expected counts, divided by the square root of the variance. 

## Package User Guide

Detailed explanations of each function in this package can be found in `/vignettes/stJointcount-vignette.Rmd`

The following code demonstrates the process of rasterization and join count analysis. In this demo, we use a human breast cancer sample that can be downloaded from 10x Genomics website [here](https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0). Users can download the related folders and files directly from the website, or they can follow first two blocks of codes in the vignette, saving them in the R temp folder.

### Data input, normalization, and clustering

The input is a directory containing the H5 file, the image data (in a subdirectory ‘spatial’), and analysis results (in a subdirectory ‘analysis’). This directory is produced by the 10x Genomics Spaceranger pipelines.

```r
sample <- spatialDataPrep(outDir)
```

This command will load a Seurat object and add the results of Spaceranger's graph-based clustering analysis to the object's metadata. 

**Note:** Users can skip this step if they want to use the data with custom labels attached. In this case, users should load their object and add label assignments to the object metadata under the title "Cluster". 

```r
#customlabels is a vector of label assignments for each spatial feature in the object
sample$Cluster <- customlabels
```

Currently, the package only accepts integer labels (i.e. 1, 2, 3, ...). We therefore suggest users with categorical labels create a corresponding integer key. Labels can be checked by `head(sample@meta.data)` and can be changed as needed.

### Resolution calculation

```r
resolutionList <- resolutionCalc(sample)
```

This function calculates the optimal resolution parameters for rasterization of the sample. 

### Rasterization

```r
mosaicIntegration <- rasterizeEachCluster(sample)
mosaicIntPlot(mosaicIntegration)
```
This function completes rasterization and label coding of the sample.

<p align="center"><img src="https://github.com/Nina-Song/stJoincount/blob/master/inst/extdata/rasterization.png" height="400"></p>

### Joint count analysis

```r
joincount.result <- joincountAnalysis(mosaicIntegration)
```

This command performes multi-categorical join count analysis of the rasterized sample. The image below shows a portion of the results.

<p align="center"><img src="https://github.com/Nina-Song/stJoincount/blob/master/inst/extdata/joincountResult.png" height="400"></p>

### Results visualization

```r
heatmapMatrix <- zscoreMatrix(sample, joincount.result)

zscorePlot(heatmapMatrix)
```

This command provides a heatmap of z-scores resulting from the join count analysis for all possible label pairs.

<p align="center"><img src="https://github.com/Nina-Song/stJoincount/blob/master/inst/extdata/zscoreHeatmap.png" height="400"></p>


### Installation
```r
library(devtools)
install_github("Nina-Song/stJoincount")
```  

### Contact
songjiar@usc.edu
