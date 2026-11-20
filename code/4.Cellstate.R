args <- commandArgs(TRUE)
method <- args[1]
outdir<- args[2]
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
source("Cellstate_lib.R")
load("obj.Rda")
geneset_file<-"ref.yaml"
ref <- geneset_file %>% tools::file_path_as_absolute() %>% yaml::yaml.load_file()


dir.create(outdir, FALSE, TRUE)
outdir <- tools::file_path_as_absolute(outdir)
setwd(outdir)

Message("Filter invalid genesets...")
ref <- Filter(length, ref)
for (i in names(ref)) {
  ref[[i]] <- gsub("\\s*$", "", ref[[i]]) %>% 
    intersect(rownames(GetAssayData(object)))
}
ref <- Filter(length, ref)
print(str(ref))

Message("Run CellState. Method: ", method)
if (method == 'mean') {
  mat <- GetAssayData(object)
  score <- sapply(ref, function(g) Matrix::colMeans(expm1(mat[rownames(mat) %in% g, , drop = FALSE])))
} else if (method == "AddModuleScore") {
  object <- AddModuleScore(object, features = ref, name = "CellState_")
  score <- object[[]][, paste0("CellState_", seq_along(ref))] %>% 
    as.matrix()
  colnames(score) <- names(ref)
  object@meta.data[, paste0("CellState_", seq_along(ref))] <- NULL
} else if (method == "AUCell") {
  library(AUCell)
  library(BiocParallel)
  library(GSEABase)

  geneSets <- sapply(names(ref), function(i) GeneSet(ref[[i]], setName = i))

  BPPARAM <- MulticoreParam(workers = 8, force.GC = TRUE, RNGseed = 1234)
  exprMatrix <- GetAssayData(object, layer = "counts")
  cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE, BPPARAM = BPPARAM)
  cells_AUC <- AUCell_calcAUC(ref, cells_rankings, aucMaxRank = ceiling(top_gene * nrow(cells_rankings)), nCores = length(ref))
  score <- getAUC(cells_AUC) %>% t()
} else if (method == "GSVA") {
  library(GSEABase) 
  library(GSVA)
  exprMatrix <- GetAssayData(object, layer = "counts")
  score <- gsva(exprMatrix, ref, method = "ssgsea", kcdf = "Poisson")
  score <- as.matrix(Matrix::t(score))
} else if (method == "UCell") {
  library(UCell)
  library(BiocParallel)
  BPPARAM <- MulticoreParam(workers = 8, force.GC = TRUE, RNGseed = 1234)
  score <- ScoreSignatures_UCell(GetAssayData(object, layer = "counts"), features = ref, name = "", BPPARAM = BPPARAM, maxRank = max(lengths(ref)))
} else {
  stop("Unsupport method: ", method)
}

score.df <- as.data.frame(score) %>%
  mutate(CellState = names(ref)[apply(score, 1, which.max)])
saveRDS(score.df, file = "score_df.Rds")

object@meta.data[, colnames(score.df)] <- score.df
object$CellState <- factor(object$CellState, levels = names(ref))

if (!"Groups" %in% colnames(object@meta.data)) {
  object$Groups <- obj$group
}
if (!"Groups" %in% colnames(obj@meta.data)) {
  object$Groups <- object$orig.ident
  object@misc$color.group <- object@misc$color.sample
}
if ("group" %in% colnames(object@meta.data)) {
  object$Groups <- object$group
}

# Output
Message(">>>> Write Cells.CellState.list.xls")
data <- FetchData2(object, var = c("Clusters", "Samples", "Groups", names(ref), "CellState"))
print(head(data))

genesets <- setdiff(colnames(data), c("Cells","Clusters", "Samples", "Groups", "Cluster", "CellState"))
print(genesets)
for (i in genesets) {
  thres0 <- quantile(data[[i]], probs = 0.75)
  data[, paste(i, "Select")] <- ifelse(data[, i] > thres0, "Positive", "Negative")
}
data$CellState <- NULL
WriteTable(cbind(Cells = rownames(data), data), "Cells.CellSelect.list.xls")

