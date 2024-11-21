# stJointcount: Quantifying spatial dependencies in labeled and mapped spatial transcriptomic data

Spatial dependency is the relationship between locational and attribute similarity. The measure reflects whether an attribute of a variable observed at one location is independent of values observed at neighboring locations. Positive spatial dependency exists when neighboring attributes are more similar than what could be explained by chance. Likewise, a negative spatial dependency is reflected by a dissimilarity of neighboring attributes. Join count analysis allows for quantification of the spatial dependencies of nominal data in an arrangement of spatially adjacent polygons. 

**stJoincount** facilitates the application of join count analysis to spatial transcriptomic data generated from the 10x Genomics Visium platform. This tool first converts a labeled spatial tissue map into a raster object, in which each spatial feature is represented by a pixel coded by label assignment. This process includes automatic calculation of optimal raster resolution and extent for the sample. A neighbors list is then created from the rasterized sample, in which adjacent and diagonal neighbors for each pixel are identified. After adding binary spatial weights to the neighbors list, a multi-categorical join count analysis is performed to tabulate "joins" between all possible combinations of label pairs. The function returns the observed join counts, the expected count under conditions of spatial randomness, and the variance calculated under non-free sampling. The z-score is then calculated as the difference between observed and expected counts, divided by the square root of the variance. 

## Package User Guide

Detailed explanations of each function in this package can be found in `/vignettes/stJointcount-vignette.Rmd`

Users can install `stJoincount` with:

```{r "install", eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
  }
BiocManager::install("stJoincount")
```

Examples of how to run this tool are below:  

## Preprocessing

In this vignette, we are going to use an human breast cancer spatial transcriptomics sample.
```{r}
fpath <- system.file("extdata", "humanBC.rda", package="stJoincount")
load(fpath)
head(humanBC)
```

Within the ‘extdata’ user can find a dataframe "humanBC.rda". This example data has come from a Seurat object of a human breast cancer sample. It contains the following information that is essential to this algorithm - barcode (index), cluster (they could either be categorical or numerical labels), central pixel location (imagerow and imagecol). This dataframe is simplified after combining metadata with spatial coordinates.

An example data preparation from Seurat:
```{r, eval=FALSE}
library(Seurat)
library(utils)
# Create a seurat object (this step can be skipped if user have customized label attached already, please start with "get metadata")
filename <- "../humanBC/"
spe <- Load10X_Spatial(filename)
sampleTransform <- SCTransform(spe, assay = "Spatial", verbose = FALSE)
clusterLabel <- read.csv(paste(filename, "/analysis/clustering/graphclust/clusters.csv", sep = ""), sep = ',')
if (colnames(clusterLabel)[1] == "barcode"){
  colnames(clusterLabel)[colnames(clusterLabel) == "barcode"] <- "Barcode"
  colnames(clusterLabel)[colnames(clusterLabel) == "cluster"] <- "Cluster"
}
rownames(clusterLabel) <- clusterLabel$Barcode
spe <- AddMetaData(object = sampleTransform, metadata = clusterLabel)

# Get metadata of Seurat object
metadata <- spe@meta.data
# Optional step: If user used their own labels, rename the column that containing cluster labels to "Cluster" 
names(metadata)[names(metadata) == 'label'] <- "Cluster"
coord <- GetTissueCoordinates(spe)
merged <-  merge(metadata, coord, by = 0)
humanBC <- merged[c("Barcode", "imagecol", "imagerow", "Cluster")]
rownames(humanBC) <- humanBC$"Barcode"
```

An example data preparation from SpatialExperiment object:
Standard analysis steps to get labels based on graph-based clustering can be found here [Chapter 13 Clustering | Orchestrating Spatially-Resolved Transcriptomics Analysis With Bioconductor](https://lmweber.org/OSTA-book/clustering.html)
Users can use their customized label as well.
```{r, eval=FALSE}
# Get metadata a SpaitialExperiment object
metadata <- colData(spe)
# If user used their own labels, rename the column that containing cluster labels to "Cluster" 
names(metadata)[names(metadata) == 'label'] <- "Cluster"
coord <- spatialCoords(spe)
merged <-  merge(metadata, coord, by = 0)

names(merged)[names(merged) == 'pxl_col_in_fullres'] <- "imagecol"
names(merged)[names(merged) == 'pxl_row_in_fullres'] <- "imagerow"
humanBC <- merged[c("barcode_id", "imagecol", "imagerow", "Cluster")]
rownames(humanBC) <- humanBC$barcode_id
```

## Raster processing

This tool first converts a labeled spatial tissue map into a raster object, in which each spatial feature is represented by a pixel coded by label assignment. This process includes automatic calculation of optimal raster resolution and extent for the sample.
```{r}
resolutionList <- resolutionCalc(humanBC)
resolutionList
```

```{r}
mosaicIntegration <- rasterizeEachCluster(humanBC)
```

## Join count analysis

A neighbors list is then created from the rasterized sample, in which adjacent and diagonal neighbors for each pixel are identified. After adding binary spatial weights to the neighbors list, a multi-categorical join count analysis is performed to tabulate "joins" between all possible combinations of label pairs. The function returns the observed join counts, the expected count under conditions of spatial randomness, and the variance calculated under non-free sampling.
```{r}
joincount.result <- joincountAnalysis(mosaicIntegration)
```

The z-score is then calculated as the difference between observed and expected counts, divided by the square root of the variance. A heatmap of z-scores represents the result from the join count analysis for all possible label pairs.

```{r}
matrix <- zscoreMatrix(humanBC, joincount.result)
```

## More examples

```{r}
browseVignettes("stJoincount")
```

### Contact

jiasong@coh.org
