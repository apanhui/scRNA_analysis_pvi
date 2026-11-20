library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
source("DEGene_lib.R")
load("obj.Rda")

## Find DE Gene
Idents(object) <- "seurat_clusters"
obj.markers <- FindAllMarkers(object = object, only.pos = TRUE,min.pct = 0.25, logfc.threshold = 0.25,return.thresh = 0.01, pseudocount.use = 0)
save( obj.markers, file = "markers.Rda" )

##output
ListMarker(object, obj.markers)
