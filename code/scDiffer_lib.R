dot_theme_default<-function(font_use = "Arial") {
	mytheme <- theme_bw() +
	theme(
		panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_rect(color = "#000000", size = 0.8, fill = NA),

		axis.text = element_text(color = "#000000", size = 11),
	    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
		axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
	    axis.title = element_text(color = "#000000", size = 14, face = "plain"),
		axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
	    axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
		axis.ticks = element_line(color = "#000000", size = 0.5),
	    axis.ticks.length = unit(0.1, 'cm'),

		legend.title = element_blank(),
	    legend.text = element_text(size = 11),
		plot.title = element_text(color = "#000000", size = 16, face = "plain", hjust = 0.5),
		plot.margin = unit(c(10,10,10,10), "mm")
	)

    if( font_use != "" & "extrafont" %in% installed.packages() ) {
        library(extrafont)
        library(extrafontdb)
        library(Rttf2pt1)
        # extrafont::font_import() to import fonts
        # extrafont::fonts() to check what fonts are imported
        if( font_use %in% fonts() ) {
            mytheme <- mytheme + theme(text = element_text(family = font_use))
        }
    }

    if( FALSE & font_use == "Arial" & "Cairo" %in% installed.packages() ) {
        library(Cairo)
        CairoFonts(
            regular = "Arial:style=Regular",
            bold = "Arial:style=Bold",
            italic = "Arial:style=Italic",
            bolditalic = "Arial:style=Bold Italic,BoldItalic"
        )
        source(ggsave_R)
    }

	if( FALSE & "showtext" %in% installed.packages() ) {
		library(showtext)
	    showtext_auto(enable = TRUE)
		font_add("Arial", regular = arial_ttf)
	    mytheme <- mytheme + theme(text = element_text(family = "Arial"))
	}

	mytheme

}

IfNull<-function(var, default = NULL){
		if ( is.null(var) ) {
				return(default)
		} else {
				return(var)
		}
}
AddUnderscore<-function(data){
		if ( ! is.null(data) ) {
				if ( ! exists("underscore", data) || ! exists("dash", data) ) {
						data$underscore <- rownames(data)
						data$dash <- gsub("_", "-", rownames(data))
				}
		}
		return(data)
}
FindFeaturesName<-function(object, features, col = "merge_name", fdata = NULL, ...) {
		if ( is.null(fdata) ) {
			fdata <- object@misc$fdata
		}
		if ( is.null(fdata) ) {
			return(features)
		}
		fdata <- AddUnderscore(fdata)
		new <- gsub("_", "-", features)
		if ( all(new %in% fdata$dash) ) {
				rownames(fdata) <- fdata$dash
				features <- new
		}
		Name <- fdata[features, col]
		Name[is.na(Name)] <- features[is.na(Name)]
		names(Name) <- features
		return(Name)
}
ChangeOUTName<-function(features, fdata) {
		features <- as.character(features)
		fdata <- AddUnderscore(fdata)
		if ( ! is.null(fdata) && all(features %in% fdata$dash) ) {
				underscore_id <- fdata$underscore
				names(underscore_id) <- fdata$dash
				features <- underscore_id[features]
		}
		return(features)
}
WriteTable<-function(x, file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, ...){
		write.table(x, file = file, quote = quote, sep = sep, row.names = row.names, col.names = col.names, ... )
}

FindFeaturesID<- function(object = NULL, features = NULL, unlist = TRUE, fdata = NULL, verbose = TRUE){
		if ( is.null(fdata) ) {
			fdata <- object@misc$fdata
		}
		if ( is.null(fdata) ) {
			return(NULL)
		}
		fdata <- AddUnderscore(fdata)
		if ( ! is.null(object) && all(rownames(object) %in% fdata$dash) ) {
			rownames(fdata) <- fdata$dash
		}
		fdata$id <- rownames(fdata)
		tmp <- sapply(fdata[, c('id', 'name', 'underscore'), drop = FALSE], toupper)
		tmp <- data.frame(tmp, row.names = rownames(fdata))
		features <- sapply(features, function(x) {
					x <- toupper(x)
					for ( y in tmp ) {
						i <- y == x
						if ( any(i) ) {
							return(rownames(tmp)[i])
						}
					}})
		if ( verbose ) {
			for ( i in names(features)[sapply(features, is.null)] ) {
				message("[WARNING] '", i, "' not found gene id.")
			}
		}
		if ( unlist ) features <- unlist(features)
		return(features)
}

CalAvgExp<-function(object, features = NULL, group.by = NULL, assay = NULL, slot = "data",
                      is.expm1 = ifelse(slot=="data", TRUE, FALSE), is.return = FALSE,
                      is.reverse = FALSE, is.bulk = FALSE,
                      row.type = 'Gene', outfile = paste0("All", row.type, ".avg_exp.xls"),
                      preftext = 'Cluster', output_name = TRUE
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
				if ( output_name ) {
						Gene_name <- FindFeaturesName(object, rownames(mean_exp), "name")
						mean_exp  <- cbind(Gene_ID = Gene_ID, Gene_name = Gene_name, mean_exp)
						colnames(mean_exp)[1:2] <- paste0(row.type, c("_ID", "_name"))
				} else {
						mean_exp  <- cbind(Gene_ID = Gene_ID, mean_exp)
						colnames(mean_exp)[1] <- paste0(row.type, c("_ID"))
				}
				WriteTable(mean_exp, file = outfile)
		}
}

CalPctExp<-function(object, features = NULL, group.by = NULL, assay = NULL, slot = "counts",
					  is.return = FALSE, is.reverse = FALSE, is.bulk = FALSE,
					  row.type = 'Gene', outfile = paste0("All", row.type, ".avg_pct.xls"),
					  preftext = 'Cluster', output_name = TRUE
					  ){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, , drop = FALSE]
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		mean_exp <- do.call(cbind, by(colnames(object), Idents(object), function(x) {
							y <- if(is.reverse){ setdiff(colnames(object), x) } else { x }
							Matrix::rowSums(data[, y, drop = F] > 0 ) / length(y)
							}))
		if ( is.bulk ) {
				mean_exp <- cbind(bulk = Matrix::rowSums(data > 0) / ncol(data), mean_exp)
		}
		if ( is.return ) {
				return(mean_exp)
		}else{
				if ( ! is.null(preftext) ) {
						colnames(mean_exp) <- paste( preftext, colnames(mean_exp) )
				}
				Gene_ID   <- ChangeOUTName(rownames(mean_exp), object@misc$fdata)
				if ( output_name ) {
						Gene_name <- FindFeaturesName(object, rownames(mean_exp), "name")
						mean_exp  <- cbind(Gene_ID = Gene_ID, Gene_name = Gene_name, mean_exp)
						colnames(mean_exp)[1:2] <- paste0(row.type, c("_ID", "_name"))
				} else {
						mean_exp  <- cbind(Gene_ID = Gene_ID, mean_exp)
						colnames(mean_exp)[1] <- paste0(row.type, c("_ID"))
				}
				WriteTable(mean_exp, file = outfile)
		}
}

FindGroupIndex <- function(object, parameter = list(), slot = "counts", assay = NULL,
				sample_col = IfNull(parameter$sample_col, "orig.ident"), cluster_col = IfNull(parameter$cluster_col, "seurat_clusters"),
				by_sample = parameter$group$by_sample,
				by_cluster = parameter$group$by_cluster,
				by_gene = parameter$group$by_gene
				){
		metadata <- object@meta.data %>% select(sample = !! sample_col)
		for ( i in unique(metadata$sample)) {
				metadata <- metadata %>% mutate(!! i := sample %in% i )
		}
		for ( i in names(by_sample) ){
				metadata <- metadata %>% mutate(!! i := sample %in% by_sample[[i]])
		}
		metadata <- metadata %>% select(- sample)

		if ( ! is.null(by_cluster) ) {
				metadata2 <- object@meta.data %>% select(cluster = !! cluster_col)
				for ( i in names(by_cluster) ){
						metadata2 <- metadata2 %>% mutate(!! i := cluster %in% by_cluster[[i]])
				}
				metadata2 <- metadata2 %>% select(- cluster)
				metadata <- cbind(metadata, metadata2)
		}

		for ( i in names(by_gene) ) {
				gene_id <- by_gene[[i]][[1]]
				gene_id <- FindFeaturesID(object, gene_id)
				gene_thres <- by_gene[[i]][[2]]
				data <- GetAssayData(object, slot = slot, assay = assay)
				metadata <- metadata %>% mutate(!! i := data[gene_id, ] >= gene_thres[[1]] & data[gene_id, ] <= gene_thres[[2]])
		}

		rownames(metadata) <- colnames(object)
		return(metadata)
}

FindGroupMarkers <- function (object, group.data, contrast = NULL,
		return.thresh = 0.01, use.qvalue = FALSE,
		assay = NULL, features = NULL, logfc.threshold = 0.25, base = 2,
		test.use = "MAST", slot = "data", min.pct = 0.1, min.pct.logic = 'OR', min.diff.pct = -Inf,
		min.exp = 0, min.exp.logic = 'OR',
		verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
		random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
		min.cells.group = 3, pseudocount.use = 1, is.bulk = TRUE, is.cluster = TRUE, ...) {
		
		if ((test.use == "roc") && (return.thresh == 0.01)) {
				return.thresh <- 0.7
		}

		cells.1 <- rownames(group.data)[group.data[[contrast[[1]]]]]
		cells.2 <- rownames(group.data)[group.data[[contrast[[2]]]]]
		if ( length(cells.1) == 0 || length(cells.2) == 0 ) {
				message( "===> cell num : \n",
						"        ", contrast[[1]], " : ", length(cells.1), "\n",
						"        ", contrast[[2]], " : ", length(cells.2), "\n",
						"        [SKIP]")
				return(NULL)
		}

		data.1 <- CalAvgExp(object[, cells.1], is.return = TRUE, is.bulk = is.bulk) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "mean.1") %>%
				mutate(cluster = as.factor(cluster))
		data.2 <- CalAvgExp(object[, cells.2], is.return = TRUE, is.bulk = is.bulk) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>% 
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "mean.2") %>%
				mutate(cluster = as.factor(cluster))
		pct.1 <- CalPctExp(object[, cells.1], is.return = TRUE, is.bulk = is.bulk) %>% round(digits = 3) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>% 
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "pct.1") %>%
				mutate(cluster = as.factor(cluster))
		pct.2 <- CalPctExp(object[, cells.2], is.return = TRUE, is.bulk = is.bulk) %>% round(digits = 3) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>% 
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "pct.2") %>%
				mutate(cluster = as.factor(cluster))
		c1 <- table(Idents(object)[cells.1])
		c2 <- table(Idents(object)[cells.2])
		if ( is.bulk ) {
				c1 <- c(bulk = length(cells.1), c1)
				c2 <- c(bulk = length(cells.2), c2)
		}
		data <- full_join(x = data.1, y = data.2, by = c("gene", "cluster")) %>%
				left_join(y = pct.1, by = c("gene", "cluster") ) %>%
				left_join(y = pct.2, by = c("gene", "cluster") ) %>%
				mutate(cells.1 = c1[as.character(cluster)], cells.2 = c2[as.character(cluster)]) %>%
				mutate(cluster = factor(cluster, levels = if ( is.bulk ) c("bulk", levels(object)) else levels(object))) %>%
				mutate(avg_logFC = log2((mean.2 + pseudocount.use)/(mean.1 + pseudocount.use)) / log2(base)) %>%
				mutate(p_val = 1, p_val_adj = 1, sig = "nosig")
				
		data[is.na(data)] <- 0
		rm(data.1, data.2, pct.1, pct.2, c1, c2)

		if ( is.bulk && is.cluster ) {
				idents.all <- c("bulk", levels(droplevels(object)))
		} else if ( is.bulk ) {
				idents.all <- "bulk"
		} else if ( is.cluster ) {
				idents.all <- levels(droplevels(object))
		}
		data <- data %>% dplyr::filter(cluster %in% idents.all)

		if ( is.null(features) ) {
				pop <- switch(min.pct.logic, AND = `&`, OR = `|`, XOR = xor, `|`)
				eop <- switch(min.exp.logic, AND = `&`, OR = `|`, XOR = xor, `|`)
				features.use <- data %>% dplyr::filter(
												mean.1 + mean.2 > 0 &
												pop(pct.1 > min.pct, pct.2 > min.pct) &
												eop(mean.1 > min.exp, mean.2 > min.exp) &
												abs(avg_logFC) > logfc.threshold / log(base))
				features.use <- split(features.use, features.use$cluster)
				features.use <- sapply(features.use, function(x) as.character(x$gene), simplify = FALSE)
				logfc.threshold <- 0
		}else{
				if ( inherits(features, 'list') ) {
					features.use <- features
				} else {
					features.use <- sapply(idents.all, function(x) features, simplify = FALSE)
				}
		}

		genes.de <- list()
		messages <- list()
		for (i in seq(idents.all)) {
				if (verbose) message("Calculating cluster ", idents.all[i])


				cells.idents <- if ( idents.all[[i]] == "bulk" ){ 
						WhichCells(object)
				} else {
						WhichCells(object, idents = idents.all[[i]])
				}
				cells.1.tmp <- intersect(cells.2, cells.idents) # switch order from T-vs-C to C-vs-T
				cells.2.tmp <- intersect(cells.1, cells.idents) # switch order from T-vs-C to C-vs-T
				if ( packageVersion("Seurat") >= as.numeric_version("4.0.0") ) {
						genes.de[[i]] <- tryCatch(expr = {
								FindMarkers(object = object, ident.1 = cells.1.tmp, ident.2 = cells.2.tmp,
											assay = assay,
											features = features.use[[idents.all[[i]]]], logfc.threshold = logfc.threshold, 
											test.use = test.use, slot = slot, min.pct = min.pct, 
											min.diff.pct = min.diff.pct, verbose = verbose, 
											only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
											random.seed = random.seed, latent.vars = latent.vars, 
											min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
											pseudocount.use = pseudocount.use, 
											base = exp(1), ## keep same as Seurat < 4.0 
											...)
								}, error = function(cond) { return(cond$message) })
				} else {
						genes.de[[i]] <- tryCatch(expr = {
								FindMarkers(object = object, ident.1 = cells.1.tmp, ident.2 = cells.2.tmp,
											assay = assay,
											features = features.use[[idents.all[[i]]]], logfc.threshold = logfc.threshold, 
											test.use = test.use, slot = slot, min.pct = min.pct, 
											min.diff.pct = min.diff.pct, verbose = verbose, 
											only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
											random.seed = random.seed, latent.vars = latent.vars, 
											min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
											pseudocount.use = pseudocount.use, 
											...)
								}, error = function(cond) { return(cond$message) })
				}
				if (class(x = genes.de[[i]]) == "character") {
						messages[[i]] <- genes.de[[i]]
						genes.de[[i]] <- NULL
				}
		}

		gde.all <- data.frame()
		for (i in seq(idents.all)) {
				if (is.null(x = unlist(x = genes.de[i]))) {
						next
				}
				gde <- genes.de[[i]]
				if (nrow(x = gde) > 0) {
						if (test.use == "roc") {
								gde <- subset(x = gde, subset = (myAUC > return.thresh | myAUC < (1 - return.thresh)))
						} else {
								gde <- gde[order(gde$p_val_adj, gde$p_val, -gde[, 2]), ]
								if ( use.qvalue ) {
										gde <- subset(x = gde, subset = p_val_adj < return.thresh)
								} else {
										gde <- subset(x = gde, subset = p_val < return.thresh)
								}
						}
				
#						if (nrow(x = gde) > 0) {
								genes.de[[i]]$sig <- "nosig"
								genes.de[[i]][rownames(gde)[gde$avg_logFC > 0], "sig"] <- "up"
								genes.de[[i]][rownames(gde)[gde$avg_logFC < 0], "sig"] <- "down"	
								genes.de[[i]]$cluster <- idents.all[i]
								genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
								gde.all <- rbind(gde.all, genes.de[[i]])
#						}
				}
		}
		if ( nrow(x = gde.all) == 0) {
				warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
				data <- data %>%
						mutate(sig = factor(sig, levels = c("up", "nosig", "down")), 
						       cluster = factor(cluster, levels = levels(data$cluster)), 
							   contrast = factor(paste0(contrast, collapse = "-vs-")) ) %>%	
						select(cluster, gene, everything())
				return(data)
		}
		gde.all <- gde.all %>%
				select(cluster, gene, p_val, p_val_adj, sig) %>%
				full_join(y = data, by = c("gene", "cluster")) %>%
				mutate( p_val     = if_else(is.na(p_val.x),     p_val.y,     p_val.x),
						p_val_adj = if_else(is.na(p_val_adj.x), p_val_adj.y, p_val_adj.x),
						sig       = if_else(is.na(sig.x),       sig.y,       sig.x)) %>%						
				mutate( sig     = factor(sig, levels = c("up", "nosig", "down")),
						cluster = factor(cluster, levels = levels(data$cluster)),
						contrast  = factor(paste0(contrast, collapse = "-vs-")) ) %>%
				select(-p_val.x, -p_val.y, -p_val_adj.x, -p_val_adj.y, -sig.x, -sig.y)

		if ((only.pos) && nrow(x = gde.all) > 0) {
				return(subset(x = gde.all, subset = gde.all[, "avg_logFC"] > 0))
		}
		if (nrow(x = gde.all) == 0) {
				warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
		}
		if (length(x = messages) > 0) {
				warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
				for (i in 1:length(x = messages)) {
						if (!is.null(x = messages[[i]])) {
								warning("When testing cluster ", idents.all[i], " :\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
						}
				}
		}
		return(gde.all)
}

DoFindGroupMarkers <- function(object, group.data, parameter = list(),
				differ_groups   = parameter$differ,
				return.thresh   = eval(parse(text = IfNull(parameter$FindMarkers$pvalue, 0.01))),
				use.qvalue      = IfNull(parameter$FindMarkers$use.q, FALSE),
				logfc.threshold = eval(parse(text = IfNull(parameter$FindMarkers$log2fc, 0.3606738))) / log2(exp(1)),
				test.use        = IfNull(parameter$FindMarkers$method, "MAST"),
				min.pct         = eval(parse(text = IfNull(parameter$FindMarkers$min_pct, 0.1))),
				min.pct.logic   = IfNull(parameter$FindMarkers$min.pct.logic, 'OR'),
				min.exp         = eval(parse(text = IfNull(parameter$FindMarkers$min_exp, 0))),
				min.exp.logic   = IfNull(parameter$FindMarkers$min.exp.logic, 'OR'),
				pseudocount.use = eval(parse(text = IfNull(parameter$FindMarkers$pseudocount.use, 1))),
				diff_type       = IfNull(parameter$diff_type, "all"),
				...){
		if ( diff_type == "all" )  {
				is.bulk <- TRUE
				is.cluster <- TRUE
		} else if ( diff_type == "no_bulk" ) {
				is.bulk <- FALSE
				is.cluster <- TRUE
		} else if ( diff_type == "only_bulk" ) {
				is.bulk <- TRUE
				is.cluster <- FALSE
		} else {
				stop()
		}
		if ( future::nbrOfWorkers() > 1 ) {
				marker <- future.apply::future_lapply(differ_groups, function(i)
								FindGroupMarkers(object, group.data, contrast = i,
										return.thresh = return.thresh,
										use.qvalue = use.qvalue,
										logfc.threshold = logfc.threshold,
										test.use = test.use,
										min.pct = min.pct, 
										min.pct.logic = min.pct.logic,
										min.exp = min.exp,
										min.exp.logic = min.exp.logic,
										pseudocount.use = pseudocount.use,
										is.bulk = is.bulk, is.cluster = is.cluster,
										...
										)
						)
				names(marker) <- sapply(differ_groups, paste, collapse = "-vs-")
				index <- sapply(marker, is.null)
				marker <- marker[!index]
		}else{
				marker <- list()
				for ( i in differ_groups ) {
						name <- paste(i, collapse = "-vs-")
						marker[[name]] <- FindGroupMarkers(object, group.data, contrast = i,
										return.thresh = return.thresh,
										use.qvalue = use.qvalue,
										logfc.threshold = logfc.threshold,
										test.use = test.use,
										min.pct = min.pct, 
										min.pct.logic = min.pct.logic,
										min.exp = min.exp,
										min.exp.logic = min.exp.logic,
										pseudocount.use = pseudocount.use,
										is.bulk = is.bulk, is.cluster = is.cluster,
										...
										)
				}
		}
		return(marker)
}


WriteDifferMarker <- function(marker, add_name = NULL, object = NULL) {
		for ( contrast_name in names(marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]
				if ( ! is.null(object) ) {
						marker[[contrast_name]]$gene <- ChangeOUTName(marker[[contrast_name]]$gene, object@misc$fdata)
						marker[[contrast_name]]$name <- FindFeaturesName(object, marker[[contrast_name]]$gene, "name")
				} else {
						marker[[contrast_name]]$name <- marker[[contrast_name]]$gene
				}
				data <- marker[[contrast_name]] %>%
						select(Cluster = cluster, GeneID = gene, GeneName = name,
								!! contrast[[1]] := mean.1, !! contrast[[2]] := mean.2,
								cells.1, cells.2, pct.1, pct.2,
								log2FC = avg_logFC, p_value = p_val, p_val_adjust = p_val_adj, significance = sig)
				WriteTable(data, file = paste(c("DifferMarker", contrast_name, add_name, "xls"), collapse = "." ))
		}
		all.list <- do.call(rbind, lapply(names(marker), function(i) { if ( nrow(marker[[i]]) == 0 ) return(NULL); cbind(contrast = i, marker[[i]][,c("cluster", "gene")]) } ) )
		WriteTable(all.list, file = "all.differ.list")
}

StatDeGene <- function(marker, group.by = c("contrast", "cluster"), do.plot = TRUE){
		group.by <- match.arg(group.by)
		if ( is.null(marker) || length(marker) == 0 ) return(NULL)
		stat <- do.call(rbind, lapply(marker, function(x) x %>%
								group_by(cluster, contrast, cells.1, cells.2, sig) %>%
								summarise(count = n()) %>%
								as.data.frame() ) )
		stat2 <- reshape2::dcast(stat, ... ~ sig, value.var = "count", fill = 0, drop = T) %>%
				mutate(up = if(exists("up")) up else 0, down = if(exists("down")) down else 0) %>% 
				select(-nosig) %>% mutate(total = up + down) %>% dplyr::filter(cells.1 > 0 & cells.2 > 0)
		for ( i in levels(stat2[[group.by]]) ) {
				name <- i
				if ( group.by == "contrast" ) {
#						name <- paste0("Contrast.", name)
				} else if ( group.by == "cluster" ) {
#						name <- paste0("Cluster", name)
				}
				contrast <- strsplit(name, "-vs-")[[1]]
				dt <- stat2 %>% dplyr::filter(.data[[group.by]] == i) %>% select(- !! group.by) %>%
					rename(!! contrast[[1]] := cells.1, !! contrast[[2]] := cells.2)
				WriteTable(dt, file = paste0("Stat.", name, ".xls"))
		}

		if ( do.plot ) {
				stat <- reshape2::melt(stat2, id.var = c("cluster", "contrast"), measure.vars = c("up", "down"), variable.name = "sig", value.name = "count" )
				plist <- list()
				for ( i in levels(stat[[group.by]]) ) {
						name <- i
						if ( group.by == "contrast" ) {
#								name <- paste0("Contrast.", name)
								group.by.inv <- "cluster"
						} else if ( group.by == "cluster" ) {
#								name <- paste0("Cluster", name)
								group.by.inv <- "contrast"
						}
						dt <- stat %>% dplyr::filter(.data[[group.by]] == i) %>%
								droplevels() %>% 
								dplyr::filter(count != 0)
						if ( nrow(dt) > 0 ) {
						plist[[i]] <- ggplot(dt, aes_string(x = group.by.inv, y = "count", fill = "sig")) + 
								geom_bar(stat = "identity", position = "dodge") +
								geom_text(aes(label = count), hjust = 0.5, vjust = -0.5, position = position_dodge(0.9)) + 
								scale_fill_discrete(NULL, drop = FALSE) +
								scale_x_discrete(NULL, drop = FALSE) + 
								ylab("Number of Genes") + theme_light() +
								theme(axis.text.x = element_text(angle=45, hjust=1))
						plist[[i]] <- plist[[i]] + bar_theme_default()
#						x.text <- unique(p$data[[as.character(p$mapping$x)[2]]])
#						w <- grid::convertWidth(unit(1, "strwidth", x.text[which.max(nchar(x.text))]), "in", T) * length(x.text)
						w <- max(7, length(unique(stat[[group.by.inv]])) * 0.8)
						ggsave(plist[[i]], file = paste0("Stat.", name, ".pdf"), width = w, height = 7)
						}
				}
		}
}


CalAvgByGrpInCls <- function(object, group.data, group_names,
				features = NULL, group.by = NULL, assay = NULL, slot = "data",
				is.expm1 = ifelse(slot=="data", TRUE, FALSE), diff_type = "all",
				outpref = "AllGene.avg_exp"
				){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, ]
		if ( is.expm1 ) data <- expm1(data)
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		all_mean_exp <- NULL
		clusters <- if ( diff_type == "all" )  {
				c("bulk", levels(Idents(object)))
			} else if ( diff_type == "no_bulk" ) {
				levels(Idents(object))
			} else if ( diff_type == "only_bulk" ) {
				"bulk"
			} else {
				stop()
			}


		Gene_ID   <- ChangeOUTName(rownames(data), object@misc$fdata)
		Gene_name <- FindFeaturesName(object, rownames(data), "name")
		for ( i in clusters ) {
				cell.cluster <- if ( i == "bulk" ) !is.na(Idents(object)) else Idents(object) == i
				mean_exp <- do.call(cbind, 
								lapply(group_names, function(j){
										cell.group <- group.data[[j]]
										if ( sum(cell.cluster & cell.group) )
											Matrix::rowMeans(data[, cell.cluster & cell.group, drop = F])
										else
											rep(0, nrow(data))
								})
							)
				colnames(mean_exp) <- group_names
				if ( ! is.null(outpref) ) {						
						out_mean_exp <- cbind(Gene_ID = Gene_ID, Gene_name = Gene_name, mean_exp)
						WriteTable(out_mean_exp, file = paste(outpref, gsub("[/ ,()]", "_", i), "xls", sep = "."))
				}
				mean_exp <- data.frame(Cluster = i, mean_exp)
				all_mean_exp <- rbind(all_mean_exp, mean_exp)
		}
		WriteTable(clusters, file = "enrich.list", col.names = FALSE)

		invisible(all_mean_exp)
}

PlotDeGeneVolcano<-function(marker, parameter = list(),
				log2fc = eval(parse(text = IfNull(parameter$FindMarkers$log2fc, 1))),
				pvalue = IfNull(parameter$FindMarkers$pvalue, 0.05),
				use.q  = IfNull(parameter$FindMarkers$use.q, FALSE),
				show_genes = NULL, gene_col = 'gene'
				){
		plist <- list()
		flist <- list(order = names(marker))
		for ( i in names(marker) ) {
				marker[[i]] <- marker[[i]] %>% filter(cells.1 > 0 & cells.2 > 0)
				marker[[i]]$cluster <- droplevels(marker[[i]]$cluster)
				flist[[i]] <- levels(marker[[i]]$cluster)
				for ( j in levels(marker[[i]]$cluster) ) {
						name <- paste0(i, "_", j)
						marker.tmp <- subset(marker[[i]], cluster == j)
						if ( ! is.null(show_genes[[i]][[j]]) ) {
								genes <- show_genes[[i]][[j]]
								marker.tmp$show_gene <- NA
								marker.tmp$show_gene[marker.tmp[[gene_col]] %in% genes]<- if(!is.null(names(genes))) names(genes) else genes
						}
						yc  <- ifelse(use.q, "p_val_adj", "p_val")
						y   <- paste0("-log10(", yc, ")")
						ylb <- ifelse(use.q, expression(-log[10](P.adjust)), expression(-log[10](Pvalue)) )
						xc  <- "avg_logFC"
						marker.tmp[[yc]][marker.tmp[[yc]] < 10^-50] <- 10^-50
						marker.tmp[[yc]][marker.tmp[[yc]] > 1] <- 1
						marker.tmp[[xc]][marker.tmp[[xc]] > 10] <- 10
						marker.tmp[[xc]][marker.tmp[[xc]] < -10] <- -10
						plist[[name]] <- ggplot(marker.tmp, aes_string(x = xc, y = y)) +
								geom_point(aes(color = sig)) +
								scale_color_manual(values = c("up" = "#F8766D", "nosig"= "grey60", "down" = "#619CFF"), name = NULL) +
								geom_hline(yintercept = -log10(pvalue),   linetype = "dashed", color = "grey40") +
								geom_vline(xintercept = c(-1,1) * log2fc, linetype = "dashed", color = "grey40") +
								labs(x = expression(log[2](FC)), y = ylb) +
								xlim(-10,10) + ylim(0,50) +
								theme_light() +
								theme(legend.key = element_rect(color = "grey70"))
						if ( ! is.null(show_genes[[i]][[j]]) ) {
								plist[[name]] <- plist[[name]] + ggrepel::geom_text_repel(aes(label = show_gene))
						}
						plist[[name]] <- plist[[name]] + dot_theme_default()
						name2 <- gsub("/", "_", name)
						ggsave(plist[[name]], file = paste0("Volcano.", name2, ".pdf"), width =7, height = 6)
				}
		}
		yaml::write_yaml(flist, file = "Volcano.list", indent.mapping.sequence = TRUE)
		invisible(plist)
}
GetAnnotation<-function(group.data, contrast, meta.data = NULL, add.group = "seurat_clusters", is.sort.by.contrast = TRUE, keep.group = FALSE) {
		if ( is.null(group.data) ) return(NULL)
		cell.1 <- group.data[[contrast[1]]]
		cell.2 <- group.data[[contrast[2]]]
		data <- data.frame(cell = rownames(group.data), row.names = rownames(group.data))
		if ( sum(cell.1&cell.2) == 0 ){
				data$contrast <- contrast[1]
				data$contrast[cell.2] <- contrast[2]
				data$contrast <- factor(data$contrast, levels = c(contrast[1], contrast[2]))
				data$group <- data$contrast
		} else {
				data[[contrast[1]]] <- sapply(cell.1, function(x) ifelse(x, contrast[1], NA))
				data[[contrast[2]]] <- sapply(cell.2, function(x) ifelse(x, contrast[2], NA))
				data$group <- "both"
				data$group[is.na(data[[contrast[1]]])] <- contrast[2]
				data$group[is.na(data[[contrast[2]]])] <- contrast[1]
				data$group <- factor(data$group, levels = c(contrast[1], "both", contrast[2]))
		}
		data <- data[cell.1|cell.2, , drop = F]
		if ( ! is.null(meta.data) && ! is.null(add.group) ) {
				name <- switch(add.group, "seurat_clusters" = "cluster", "orig.ident" = "samples", add.group)
				data[[name]] <- meta.data[rownames(data), add.group]
				if ( is.sort.by.contrast ) {
						data <- data %>% arrange(group, !! as.name(name))
				} else {
						data <- data %>% arrange(!! as.name(name), group)
				}
		} else {
				data <- data %>% arrange(group)
		}
		if (  packageVersion("dplyr") >= as.numeric_version("1.0.3") ) { ## damn it
				data <- data %>% select(-cell)
		} else {
				data <- tibble::column_to_rownames(data, var = "cell")
		}
		if ( ! keep.group ) {
				data$group <- NULL
		}
		return(data)
}
MinMax<-function (data, min, max)
{
	data2 <- data
	data2[data2 > max] <- max
	data2[data2 < min] <- min
	return(data2)
}


PlotContrastHeatmap<-function(object, object.marker, group.data = NULL,
				top.num = 5, cluster.name = "seurat_clusters", do.return = FALSE,
				is.sort.by.contrast = TRUE,
				assay = NULL, slot = "scale.data",
				scale.min = -2.5, scale.max = 2.5,
				color = c("purple", "black", "yellow"), ...
				) {
		for ( contrast_name in names(object.marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]

				annotation <- GetAnnotation(group.data, contrast, meta.data = object@meta.data, add.group = cluster.name, is.sort.by.contrast = is.sort.by.contrast)
				cells <- rownames(annotation)
				features <- object.marker[[contrast_name]] %>% group_by(cluster) %>% top_n(top.num, abs(avg_logFC))
				features <- unique(features$gene)
				if ( length(features) == 0 || length(cells) == 0 ) {
					message("PlotContrastHeatmap skip [", contrast_name, "] : no differ gene")
					next
				}

				if ( slot == 'scale.data' ) {
						if ( ! is.null(features) ) {
								if ( is.null(assay) ) assay <- DefaultAssay(object)
								DefaultAssay(object) <- assay
								if ( ! all(features %in% rownames(GetAssayData(object[[assay]], slot))) ) {
										message('Below features arenot in scale.data, recal scale.data:')
										message(paste0(features[! features %in%rownames(GetAssayData(object[[assay]], slot))], collapse = ", "))
										object <- ScaleData(object, features = features)
								}
						}
				}

				mat <- GetAssayData(object[[assay]], slot)[features, cells, drop = FALSE]
				mat <- MinMax(mat, scale.min, scale.max)
				rownames(mat) <- FindFeaturesName(object, rownames(mat), "name")
				if ( do.return ){
				} else {
						annotation_colors <- list(switch(cluster.name,
								"Groups" = object@misc[["color.group"]],
								"orig.ident" = object@misc[["color.sample"]],
								"seurat_clusters" = object@misc[["color.cluster"]]))
						names(annotation_colors) <- switch(cluster.name, "seurat_clusters" = "cluster", "orig.ident" = "samples", cluster.name)
						if ( is.null(annotation_colors[[1]]) ) annotation_colors[[1]] <- SetColor(unique(annotation[[names(annotation_colors)]]))

						color <- colorRampPalette(color)(100)
						MAX <- max(abs(mat))
						breaks <- seq(-1 * MAX, MAX, length.out = length(color))
						cluster_rows <- nrow(mat) > 1
						pheatmap::pheatmap(mat, cluster_cols = F, show_colnames = F, cluster_rows = cluster_rows,
								annotation_col = annotation, annotation_colors = annotation_colors,
								scale = "none", color = color, breaks = breaks,
								cellheight = 10, #cellwidth = 10,
								filename = paste0("Heatmap.", contrast_name, ".pdf"), ...)
				}

		}
}
