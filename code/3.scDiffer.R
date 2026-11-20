library(Seurat); library(dplyr); library(ggplot2); library(cowplot)

source("scDiffer_lib.R", chdir = T)
file<-"scDiffer_config.yaml"
parameter <- yaml::yaml.load_file(file, handlers = list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x}))

load("obj.Rda")
group.data <- FindGroupIndex(object, parameter)
CalAvgByGrpInCls(object, group.data, group_names = unique(unlist(parameter$differ)), diff_type = parameter$diff_type)

marker <- DoFindGroupMarkers(object, group.data, parameter, diff_type = parameter$diff_type)
save(marker, file = "marker.Rda")

WriteDifferMarker(marker, "all", object = object)
StatDeGene(marker, group.by = "contrast", do.plot = FALSE)
PlotDeGeneVolcano(marker, parameter)
diff.marker <- lapply(marker, function(x) filter(x, sig != "nosig"))
WriteDifferMarker(diff.marker, "diff", object = object)

PlotContrastHeatmap(object, diff.marker, group.data, top.num = parameter$plots$top, cluster.name = "seurat_clusters")

