AddUnderscore <-function (data){
	if (!is.null(data)) {
		if (!exists("underscore", data) || !exists("dash", data)) {
			data$underscore <- rownames(data)
			data$dash <- gsub("_", "-", rownames(data))
		}
	}
	return(data)
}
ChangeOUTName <-function (features, fdata){
	features <- as.character(features)
	fdata <- AddUnderscore(fdata)
	if (!is.null(fdata) && all(features %in% fdata$dash)) {
		underscore_id <- fdata$underscore
		names(underscore_id) <- fdata$dash
		features <- underscore_id[features]
	}
	return(features)
}

CalAvgExp<-function(object, features = NULL, group.by = NULL, assay = NULL, slot = "data",
                      is.expm1 = ifelse(slot=="data", TRUE, FALSE), is.return = FALSE,
                      is.reverse = FALSE, is.bulk = FALSE, outfile = "AllGene.avg_exp.xls",
                      preftext = 'Cluster'
					  ){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, , drop = FALSE]
		if ( is.expm1 ) data <- expm1(data)
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		mean_exp <- do.call(cbind, by(colnames(object), Idents(object), function(x) {
							y <- if(is.reverse){ setdiff(colnames(object), x) } else { x }
							Matrix::rowMeans(data[, y, drop = F])
							}, simplify = FALSE))
		if ( is.bulk ) {
				mean_exp <- cbind(bulk = Matrix::rowMeans(data), mean_exp)
		}
		if ( is.return ) {
				return(mean_exp)
		}else{
				if ( ! is.null(preftext) ) {
						colnames(mean_exp) <- paste( preftext, colnames(mean_exp) )
				}
				Gene_ID   <- ChangeOUTName(rownames(mean_exp), object@misc$fdata)
				Gene_name <- FindFeaturesName(object, rownames(mean_exp), "name")
				mean_exp  <- cbind(Gene_ID = Gene_ID, Gene_name = Gene_name, mean_exp)
				write.table(mean_exp, file = outfile,quote = FALSE, sep = '\t', row.names = FALSE,col.names = TRUE)
		}
}
FindFeaturesName<- function(object, features, col = "merge_name", ...) {
		if ( ! exists("fdata", object@misc) )
				return(features)
		object@misc$fdata <- AddUnderscore(object@misc$fdata)
		new <- gsub("_", "-", features)
		if ( all(new %in% object@misc$fdata$dash) ) {
				rownames(object@misc$fdata) <- object@misc$fdata$dash
				features <- new
		}
		Name <- object@misc$fdata[features, col]
		Name[is.na(Name)] <- features[is.na(Name)]
		names(Name) <- features
		return(Name)
}


ListMarker<-function(object, object.markers, outfile = "DeGene.list.xls", is.return = FALSE, is.fast = FALSE, group.by = "seurat_clusters",assay = DefaultAssay(object), slot = "data", group.by.name = "Cluster", ...){
		Targets_name <- paste0( "Target_", group.by.name )
		Others_name  <- paste0( "Other_",  group.by.name )
		Targets_mean_name <- paste0( Targets_name, "_mean" )
		Others_mean_name  <- paste0( Others_name,  "_mean" )
		Targets_pct_name  <- paste0( Targets_name, "_pct" )
		Others_pct_name   <- paste0( Others_name,  "_pct" )
		Targets_mean <- CalAvgExp(object, unique(object.markers$gene), is.return = T, is.reverse = F, assay = assay, slot = slot, group.by = group.by, ...) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = Targets_mean_name) %>%
				mutate(cluster = factor(cluster))
		Others_mean  <- CalAvgExp(object, unique(object.markers$gene), is.return = T, is.reverse = T, assay = assay, slot = slot, group.by = group.by, ...) %>%
				reshape2::melt(varnames = c("gene", "cluster"), value.name = Others_mean_name) %>%
				mutate(cluster = factor(cluster))
		Name <- FindFeaturesName(object, unique(object.markers$gene), "name", is.fast = is.fast)
		marker_list <- object.markers %>% left_join(y = Targets_mean) %>% left_join(y = Others_mean) %>%
				mutate(Log2FC = log2(!! as.name(Targets_mean_name) / !! as.name(Others_mean_name)), name = Name[gene]) %>%
				select(!! Targets_name := cluster, "Gene ID" = gene, "Gene Name" = name,
						!! Targets_pct_name := pct.1, !! Others_pct_name := pct.2,
						!! Targets_mean_name, !! Others_mean_name, Log2FC,
						Pvalue = p_val, Qvalue = p_val_adj)
		if ( is.null(outfile) || is.return ){
				return(marker_list)
		}else{
				marker_list[["Gene ID"]] <- ChangeOUTName(marker_list[["Gene ID"]], object@misc$fdata)
				write.table(marker_list, file = outfile, quote = FALSE, sep = '\t', row.names = FALSE,col.names = TRUE)
		}
}

