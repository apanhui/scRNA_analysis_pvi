library(Seurat); library(harmony); library(dplyr); library(ggplot2); library(patchwork)

#################
### rawdata ###
#################
#"barcodes.tsv"   "genes.tsv"      "matrix.mtx"

## Create Seurat object
data_name = c("CK","As")
object.list <- list()
for (i in seq(data_name)) {
	mat <- Read10X(data.dir = data_name[i], gene.column = 1)
	object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = "RNA")
}
## merge Seurat object
object <- merge(x = object.list[[1]], y = unlist(object.list[-1]),add.cell.ids = data_name)

##Filter
object <- subset(object,  subset = nFeature_RNA < 98000 &  nCount_RNA < 40000)
##Doublet
object[["DF"]] <- "Singlet"
filter.cells <- readLines("doublet.xls")
object@meta.data[filter.cells, "DF"] <- "Doublet"
object <- subset(object, subset = DF == "Singlet")

## integration
object <- NormalizeData(object, normalization.method = "LogNormalize",scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst",nfeatures = 2000)
object <- ScaleData( object, vars.to.regress = NULL,features = rownames(object))

object <- RunPCA(object,assay = "RNA",npcs = 50,features = VariableFeatures(object), verbose = FALSE)
object <- RunHarmony(object, group.by.vars = "orig.ident", project.dim = FALSE, assay.use = DefaultAssay(object))

object <- RunTSNE(object, dims = 1:50,dim.embed = 2L , perplexity = 30, reduction.name = "tsne",reduction="harmony")
object <- RunUMAP(object, dims = 1:50, umap.method = "uwot",n.neighbors = 30,n.components = 2L , reduction.name = "umap",reduction = "harmony")

## Find clusters
object <- FindNeighbors(object, reduction = "harmony", dims = 1:50,prune.SNN = 1/15)
object <- FindClusters(object, resolution = 0.5,algorithm = 1, temp.file.location = getwd(),method = "igraph")
object$orig.ident <- factor(object$orig.ident)
str(object@meta.data)
##color
n_clusters <- length(unique(object$seurat_clusters))
cluster_levels <- as.character(sort(unique(object$seurat_clusters)))
object@misc$color.cluster <- rainbow(n_clusters)
names(object@misc$color.cluster) <- cluster_levels

n_sample<-length(unique(object$orig.ident))
sample_levels <- as.character(sort(unique(object$orig.ident)))
object@misc$color.sample<-rainbow(n_sample)
names(object@misc$color.sample)<-sample_levels
## Draw t-SNE UMAP plot
p1 <- DimPlot(object, reduction = 'tsne', group.by = "orig.ident", cols=object@misc$color.sample, label = FALSE)
p2 <- DimPlot(object, reduction = 'tsne', group.by = "seurat_clusters", cols=object@misc$color.cluster, label = TRUE)
ggsave(p1+p2, file = "TSNE.pdf",width = 14, height = 6)
p1 <- DimPlot(object, reduction = 'umap', group.by = "orig.ident", cols=object@misc$color.sample, label = FALSE)
p2 <- DimPlot(object, reduction = 'umap', group.by = "seurat_clusters", cols=object@misc$color.cluster, label = TRUE)
ggsave(p1+p2, file = "UMAP.pdf",width = 14, height = 6)

## Save data object
DefaultAssay(object) <- "RNA"
save(object, file = "obj.Rda")

