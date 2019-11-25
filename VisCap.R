# VisCap: Visualize normalized capture coverage
# Generates heatmap of exon coverage from a directory of sample_interval_summary files
# By Trevor Pugh, March 2012 - April 2013

###########
#Libraries#
###########

library("cluster")
library("gplots")
library("zoo")

##########
#Defaults#
##########

svn.revision                <- "$Id: VisCap.R 1372 2013-04-25 13:53:52Z tp908 $"
source('VisCap.cfg')
#lane_dir                    <- "\\\\rfanfs.research.partners.org\\NGS"
#out_dir                     <- "\\\\rfanfs.research.partners.org\\gigpad_clinical\\Facilities\\Laboratory of Molecular Medicine_4171303\\VisCap_Outputs"

#lane_dir                    <- "\\\\rfanfs.research.partners.org\\NGS"
#out_dir                     <- "\\\\rfanfs.research.partners.org\\gigpad_clinical\\Facilities\\Laboratory of Molecular Medicine_4171303\\VisCap_Outputs"
# interval_list_dir           <- "\\\\Sfa6\\lmm$\\DEVELOPMENT\\ACTIVE_DEVELOPMENT\\NEXT_GEN_COPY_NUMBER\\VisCap\\interval_lists"   
#explorer_file               <- "C:\\Windows\\explorer.exe"
#cov_file_pattern            <- ".target.cov.sample_interval_summary$"
#cov_field                   <- "_total_cvg"
#interval_file_pattern       <- ".interval_list$"
#ylimits                     <- c(-2, 2)
#iqr_multiplier              <- 3
#threshold.min_exons         <- 1
#iterative.calling.limit     <- 0        #Set to 0 to iterate until all failed samples are removed
#infer.batch.for.sub.out_dir <- TRUE     #Set to FALSE to prompt users for output directory
#clobber.output.directory    <- FALSE    #Set to FALSE to stop run when output directory already exists

#Setting a path to data skips prompts. Set to FALSE for deployment.
#dev_dir    <- "\\\\Sfa6\\lmm$\\DEVELOPMENT\\ACTIVE_DEVELOPMENT\\NEXT_GEN_COPY_NUMBER\\VisCap\\test_data"
dev_dir     <- FALSE

###########
#Functions#
###########

#Modified winDialog function to run in non-interactive mode
winDialog.nonint <- function (type = c("ok", "okcancel", "yesno", "yesnocancel"), message) 
    {
        #if (!interactive())
        #    stop("winDialog() cannot be used non-interactively")
        type <- match.arg(type)
        res <- .Internal(winDialog(type, message))
        if (res == 10L) 
            return(invisible(NULL))
        c("NO", "CANCEL", "YES", "OK")[res + 2L]
    }

make_matrix_from_cov_files <- function(lane_dir, cov_file_pattern, cov_field) {
    #Join all cov files into a single matrix
    filenames <- list.files(lane_dir, full.names=TRUE, pattern=cov_file_pattern, recursive=TRUE)
    print(filenames)
    for(file in filenames) {
        tab      <- read.table(file, header=TRUE, sep ="\t", row.names=1)
        col_name <- colnames(tab)[grep(cov_field, colnames(tab))]

        #Make new matrix labelled with appropriate identifier
        #new_head <- gsub(cov_field, "", col_name)                  #sample identifier only
        new_head  <- basename(gsub(cov_file_pattern, "", file))      #derive name from original filename
        new       <- matrix(tab[,col_name], dimnames=list(rownames(tab), new_head) )
        
        #If this is the first file, make a new matrix. Otherwise, join it to the existing matrix using genome coordinates
        if(file == filenames[1]) {
            mat.cov <- new
        } else {
            mrg               <- merge(new, mat.cov, by = "row.names", all = TRUE)
            mat.cov           <- as.matrix(mrg[-1])
            rownames(mat.cov) <- mrg[,1]
        }
    }

    #Remove any rows with NA entries
    mat.cov <- mat.cov[complete.cases(mat.cov),]
    #TODO: Provide warning if n% of probes are discarded (i.e. multiple panels are present)
    
    # Extract fraction of total coverage assigned to each exon in each sample
    # Note two possible methods, depending on preference for X-chromosome handling.
    # mat.cov.totals <- colSums(mat.cov) #use X-chromosome for normalization
    mat.cov.totals <- colSums(mat.cov[grep("X",rownames(mat.cov), invert=TRUE),]) # do not use X-chromosome for normalization
    mat.frac_cov <- sweep(mat.cov, 2, mat.cov.totals, "/")

    #Order columns alphanumerically
    mat.frac_cov <- mat.frac_cov[,order(colnames(mat.frac_cov))]
    
    return(mat.frac_cov)
}

load_interval_names <- function(interval_list_dir, interval_file_pattern) {
    #Read interval list files and make lookup_table for genome coords and interval names
    lookup = c()
    filenames = list.files(interval_list_dir, full.names=TRUE, pattern=interval_file_pattern, recursive=FALSE)
    for(file in filenames) {
        tab = read.table(file, header=FALSE, comment.char = "@", sep="\t", stringsAsFactors=FALSE)
	if(ncol(tab) < 4)
	{stop("The interval file should have at least 4 columns, chromosome, start, end, and interval name, all tab separated.")}
        colnames(tab) <- c("chr", "start", "end", "interval_name")
        rownames(tab) <- paste(tab$chr, ":", tab$start, "-", tab$end, sep = "")
        tab$interval_file <- file
        lookup <- rbind(lookup, tab)
    }
    return(lookup)
}

annotate_interval_names <- function(coords, interval_lookup) {
    #Add name column to mat containing interval name from lookup table, if provided
    if(is.null(interval_lookup)) {
    labels <- coords
    } else {
        matches <- match(coords, rownames(interval_lookup))
        labels <- interval_lookup[matches, "interval_name"]
    }
    return(labels)
}

divide_by_batch_median <- function(mat) {
    #Fractional coverage normalization by median across each exon
    rmeds <- apply(mat, 1, median)         ## calculate row medians
    mat_norm <- sweep(mat, 1, rmeds, "/")  ## divide each entry by row median
    mat_norm <- mat_norm[complete.cases(mat_norm),] ##remove NA entries
	return(mat_norm)
}

heatmap_by_chrom <- function(mat, analysis_name, ylimits, out_dir) {
    #Limit matrix for plotting
    mat[mat < min(ylimits)] <- min(ylimits)
    mat[mat > max(ylimits)] <- max(ylimits)

	# Plot heatmap for each chromosome, compare each exon's value/median ratio across samples, save to files
	pdf(file=paste(out_dir, "/", analysis_name, ".pdf", sep=""))

    #Per-chromosome plots
	if(length(interval_lookup$chr[grep("chr",interval_lookup$chr)])!=0){chroms <- c("all",paste("chr",1:22,sep=""),"chrX","chrY", "chrMT", "chrM")}else{chroms = c("all", 1:22,"X","Y", "MT", "M")}
	for(chr in chroms){
        if(chr == "all") {
            main_title <- "All Chromosomes"
            matches <- rownames(mat)
        } else {
            main_title <- paste("Chromosome",chr)
            matches <- grep(paste("^",chr,":",sep=""), rownames(mat))
        }
		if(length(matches) > 1) {
			par("cex.main" = 0.5)
            #Set heatmap color scale
            steps = 100
            color_scale = bluered(steps)
            color_breaks = seq(ylimits[1], ylimits[2], by=(ylimits[2] - ylimits[1])/steps)
            #color_breaks = seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out=steps+1)

            #exons by samples
            #heatmap.2(mat[matches,], Rowv=NA, Colv=NA, scale="none", cexCol=0.5, cexRow=0.3, col=color_scale, breaks=color_breaks, dendrogram="none", main=title, symbreaks=TRUE, trace="none")
            #samples vs exons
            xlabels <- gsub("__", "\n", colnames(mat))
            heatmap.2(
                        t(mat[matches,]),

                        Rowv       = NA,
                        Colv       = NA,
                        scale      = "none",
                        cexCol     = 0.5,
                        cexRow     = 0.5,
                        col        = color_scale,
                        breaks     = color_breaks,
                        dendrogram = "none",
                        symbreaks  = TRUE,
                        trace      = "row",
                        tracecol   = "black",
                        main       = main_title,
                        labRow     = xlabels,
                        labCol     = NA
            )
		}
	}
    dev.off()
}

score_boxplot <- function(bplot, llim= log2(1/2), ulim=log2(3/2)) {
    fail_lower_thresh <- bplot$stats[1,] < llim
    fail_upper_thresh <- bplot$stats[5,] > ulim
    fail <- as.logical(fail_lower_thresh + fail_upper_thresh)
    qc_string <- gsub(TRUE, "FAIL", gsub(FALSE, "PASS", fail))
    return(qc_string)
}

boxplot_cnv_matrix <- function(nmat, bplot_name, out_dir, ylimits, iqr_multiplier) {
    plot_nmat <- nmat
    plot_nmat[nmat < min(ylimits)] <- min(ylimits)
    plot_nmat[nmat > max(ylimits)] <- max(ylimits)
    pdf(file=paste(out_dir, "/", bplot_name, ".pdf", sep=""))
    par(las=2, mar=c(12,4,4,2))
    xlabels = sub("__", "\n", colnames(plot_nmat))
    bplot <- boxplot(plot_nmat, range=iqr_multiplier, ylim=ylimits, srt=90, pch=16, cex.axis=0.6, names=xlabels, ylab="Log2 ratio sample/batch median")
    bplot$names <- sub("\n", "__", bplot$names)
    dev.off()
    
    #Write boxplot information to separate file
    bplot$qc <- score_boxplot(bplot)
    bplot.out <- cbind(round(t(bplot$stats), 3), bplot$qc)
    rownames(bplot.out) <- bplot$names
    colnames(bplot.out) <- c("del_threshold", "Q1", "median", "Q3", "amp_threshold", "qc")
    write.table(bplot.out, paste(out_dir, "/", bplot_name, ".xls", sep=""), quote = FALSE, sep = "\t", col.names = NA)
    return(bplot)
}

filter_cnvs <- function(segs, threshold.min_exons, threshold.cnv_log2_cutoffs) {
    #Apply consecutive exon filter
    fsegs <- segs[(segs[,"CNV"] != 0) & (segs[,"Interval_count"] >= threshold.min_exons), , drop = FALSE]

    #Apply zero in normal range filter
    fsegs <- fsegs[(as.numeric(fsegs[,"Loss_threshold"]) < 0) & (as.numeric(fsegs[,"Gain_threshold"]) > 0),, drop = FALSE]

    #Apply hard CNV threshold
    fsegs <- fsegs[(as.numeric(fsegs[,"Median_log2ratio"]) < min(threshold.cnv_log2_cutoffs)) | (as.numeric(fsegs[,"Median_log2ratio"]) > max(threshold.cnv_log2_cutoffs)),, drop = FALSE]

    return(fsegs)
}

call_cnvs <- function(nmat, ylimits, interval_lookup, threshold.min_exons, iqr_multiplier, threshold.cnv_log2_cutoffs, out_dir) {
    #Plot boxplots to visualize ranges used to detect copy number variation
    bplot <- boxplot_cnv_matrix(nmat, "QC_cnv_boxplot", out_dir, ylimits, iqr_multiplier)
    
    #Set outlier thresholds by distribution test
    #Use boxplot$stats so whiskers in pdf accurately represent thresholds used
    lbound <- round(bplot$stats[1,], 3)
    ubound <- round(bplot$stats[5,], 3)
    batch_size <- dim(bplot$stats)[2]

    #Use hard-threshold, if it is greater than boxplot thresholds
    lbound[lbound > min(threshold.cnv_log2_cutoffs)] <- min(threshold.cnv_log2_cutoffs)
    ubound[ubound < max(threshold.cnv_log2_cutoffs)] <- max(threshold.cnv_log2_cutoffs)

    #Make matrix of thresholds
    lbound_mat <-  matrix(rep(lbound, dim(nmat)[1]), ncol=length(lbound), byrow=TRUE, dimnames=list(rownames(nmat),colnames(nmat)))
    ubound_mat <-  matrix(rep(ubound, dim(nmat)[1]), ncol=length(lbound), byrow=TRUE, dimnames=list(rownames(nmat),colnames(nmat)))

    #Flag outliers
    nmat_loutliers <- (nmat < lbound_mat) + 0 #Adding zero converts TRUE/FALSE to 1/0
    nmat_uoutliers <- (nmat > ubound_mat) + 0 #Adding zero converts TRUE/FALSE to 1/0
    #TODO: Estimate number of copies gained or lost (i.e. support -2 and +n)
    
    #Make tracking matrix of all zero values then subtract copies lost and add copies gained
    nmat_cnvs <- matrix(data = 0, nrow = dim(nmat)[1], ncol = dim(nmat)[2], dimnames = list(rownames(nmat),colnames(nmat)))
    nmat_cnvs <- nmat_cnvs - nmat_loutliers + nmat_uoutliers
    gene_exon <- annotate_interval_names(rownames(nmat_cnvs), interval_lookup)
    write.table(cbind.data.frame(gene_exon, nmat_cnvs), paste(out_dir, "/", "cnv_boxplot_outliers", ".xls", sep=""), quote = FALSE, sep = "\t", col.names = NA)
    
    #Merge consecutive calls and write out to file
    all_fsegs <- c()
    for(id in colnames(nmat_cnvs)) {
        segs <- c()
        segs.header <- c("Sample", "CNV", "Genome_start_interval", "Genome_end_interval", "Start_interval","End_interval", "Interval_count", "Min_log2ratio", "Median_log2ratio", "Max_log2ratio", "Loss_threshold", "Gain_threshold", "Batch_size")
        if(length(interval_lookup$chr[grep("chr",interval_lookup$chr)])!=0){chroms <- c(paste("chr",1:22,sep=""),"chrX","chrY", "chrMT", "chrM")}else{chroms = c(1:22,"X","Y", "MT", "M")}
        for(chr in chroms){
            matches <- grep(paste("^",chr,":",sep=""), rownames(nmat))
            if(length(matches) > 1) {
                #Segmentation: Detect runs of consecutive copy number calls
                rl <- rle(nmat_cnvs[matches, id])
                values <- rl$values
                lengths <- rl$lengths
                starts <- c(rownames(nmat_cnvs[matches,])[1], names(rl$lengths[1:length(rl$lengths)-1]))
                ends <- names(rl$values)

                #Calculate rounded log2 ratios of intervals involved
                log2s <- lapply(1:length(starts), function(x) round(nmat[(match(starts[x], names(nmat[,id])):match(ends[x], names(nmat[,id]))), id], 3))
                log2s_min <- unlist(lapply(log2s, min))
                log2s_med <- unlist(lapply(log2s, median))
                log2s_max <- unlist(lapply(log2s, max))

                #Report genome coordinate range
                coordinates_part1 <- data.frame(strsplit(starts, "-"), stringsAsFactors = FALSE)[1,]
                coordinates_part2 <- data.frame(strsplit(ends, "-"), stringsAsFactors = FALSE)[2,]                
                coordinates <- paste(coordinates_part1, coordinates_part2, sep="-")
                
                #Lookup interval names
                start_names <- annotate_interval_names(starts,  interval_lookup)
                end_names <- annotate_interval_names(ends,  interval_lookup)
                
                #Add threshold columns, ensure consistent header
                segs <- rbind(segs, cbind(rep(id, length(values)), values, starts, ends, start_names, end_names, lengths, log2s_min, log2s_med, log2s_max, lbound_mat[1,id], ubound_mat[1,id], batch_size))
                colnames(segs) <- segs.header
            }
        }

        #Handle case where there are no cnvs on any chromosomes in any samples
        if(is.null(segs)) {
            segs <- matrix(data=rep(NA,length(segs.header)), ncol=length(segs.header), dimnames=list("", segs.header))[-1,,drop=FALSE]
        }
        
        #Filter segments
        fsegs <- filter_cnvs(segs, threshold.min_exons, threshold.cnv_log2_cutoffs)
        
        #Convert CNV type values to text
        fsegs[,"CNV"] <- gsub("^1$", "Gain", fsegs[,"CNV"])
        fsegs[,"CNV"] <- gsub("^-1$", "Loss", fsegs[,"CNV"])

        #Convert infinite values to large, numerical value
        large_value <- "10"
        num_cols <- c("Min_log2ratio", "Median_log2ratio", "Max_log2ratio", "Loss_threshold", "Gain_threshold")
        fsegs[,num_cols] <- gsub("Inf", large_value, fsegs[,num_cols])

        #Number CNVs for later labeling on visual output
        if(dim(fsegs)[1] > 0) {
            CNV_id <- seq(1, dim(fsegs)[1])
        } else {
            CNV_id <- c()
        }
        fsegs <- cbind(fsegs[,1, drop=FALSE],
                       CNV_id,
                       fsegs[,2:dim(fsegs)[2], drop=FALSE])

        #Write output for each sample
        write.table(fsegs, paste(out_dir, "/", id, ".cnvs.xls", sep=""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        #Add fsegs to master fsegs tracking matrix
        all_fsegs <- rbind(all_fsegs, fsegs)
    }
    return(list(bplot, all_fsegs))
}

space_label_names <- function(label_names) {
    spacer <- "     "
    if(length(label_names) > 1) {
        alt_positions <- seq(from=2, to=length(label_names), by=2)
        label_names[alt_positions] <- paste(label_names[alt_positions], spacer)
    }
    return(label_names)
}

exon_plot_per_sample <- function(nmat, ylimits, interval_lookup, cnv_bplot_and_calls, out_dir) {
    #Extract thresholds and cnvs, determine intervals from genome coordinates
    nmat_bplot <- cnv_bplot_and_calls[[1]]
    nmat_fsegs <- cnv_bplot_and_calls[[2]]

    #Limit matrix for plotting
    nmat[nmat < ylimits[1]] <- ylimits[1]
    nmat[nmat > ylimits[2]] <- ylimits[2]

    for(name in colnames(nmat)) {
        #Get exon names from interval lists
        interval_names <- annotate_interval_names(rownames(nmat), interval_lookup)

        #Plot all intervals, then plot by chromosome
        pdf(file=paste(out_dir, "/", name, ".plot.pdf", sep=""))                 
		if(length(interval_lookup$chr[grep("chr",interval_lookup$chr)])!=0){chroms <- c(".*",paste("chr",1:22,sep=""),"chrX","chrY")}else{chroms = c(".*", 1:22,"X","Y")}
		for(chr in chroms) {
            #Grep rownames with genome coordinates for desired chromosome, support chr1: and 1: formats
			matches <- grep(paste("(^chr", chr, "|^", chr, "):", sep=""), rownames(nmat))
			if(length(matches) > 1) {
                nmat_chr <- nmat[matches, name, drop=FALSE]
                #Different title and axes labeling for All vs individual chromosomes
                if(chr == ".*") {
                    title <- "All chromosomes"
					chrom_names <- sapply(rownames(nmat_chr), function(x) strsplit(x, ":")[[1]][1])
					labels_rle <- rle(chrom_names)
                    #TODO: Find better method for spacing chromosome labels
                    labels_rle$values <- space_label_names(labels_rle$values)
                    #labels_rle$values <- rep("", length(labels_rle$values-1))
                    labels_size <- 0.6
                    draw_grid_lines <- FALSE
                    label_called_CNVs <- FALSE
                } else {
                    title <- paste("Chromosome", chr)
					#Convert genome coordinates to unique gene names
					exon_names <- annotate_interval_names(rownames(nmat_chr), interval_lookup)
					gene_names <- sapply(exon_names, function(x) strsplit(x, "_")[[1]][1])
					#Replace NA gene_names with nearest non-NA value
					gene_names <- na.locf(gene_names)
					labels_rle <- rle(gene_names)
                    labels_size <- 1
                    draw_grid_lines <- TRUE
                    label_called_CNVs <- TRUE
                }

                #Plot formatting and command
                shift <- -0.5 #shift for plotting gene guide lines
                par(pch=16)
				plot(nmat_chr, xlim=c(1 + shift,length(nmat_chr)), ylim=ylimits, main=title, ylab = "Log2 ratio sample/batch median", xaxt="n", xlab="")
                mtext(name, line=3, cex=0.5, adj=0)
                
                #Plot gene or chomosome grid lines
				section_ends <- (unlist(lapply(seq(1:length(labels_rle$lengths)), function(x) sum(labels_rle$lengths[1:x]))))
                section_starts <- c(1, section_ends[-length(section_ends)] + 1)
                grid_marks <- section_starts + shift
                grid_labels <- c(labels_rle$values)
                if(draw_grid_lines == TRUE) {
                    abline(v=grid_marks, col="grey")
                }
				axis(1, las=3, at=grid_marks, labels=grid_labels, cex.axis=labels_size)

                #Plot guidelines and thresholds
                name_bplot <- nmat_bplot$stats[,nmat_bplot$names == name]
                abline(h=0, col="black")                          #zero line
                abline(h=c(log2(1/2), log2(3/2)), col="grey")     #expected log2 ratios for single copy loss and gain
                abline(h=name_bplot[c(1,5)], lty=2, col="black")  #boxplot thresholds
                abline(h=threshold.cnv_log2_cutoffs, col="black") #hard cnv log2 ratio thresholds
                
                #Mark data points outside of ylimits and scaled to fit plot
                off_scale <- as.numeric(which((nmat_chr <= ylimits[1]) | (nmat_chr >= ylimits[2])))
                points(off_scale, nmat_chr[off_scale], pch=22)
                                            
                #Color CNV data points
                types     <- matrix(byrow=TRUE, ncol=4,
                                    dimnames=list(c(),
                                        c("type", "col", "ypos.line", "ypos.label")),
                                    data=c(
                                          "Loss", "blue", -2,         -2.1,
                                          "Gain", "red",   2,          2.1))

                segs <- nmat_fsegs[(nmat_fsegs[,"Sample"] == name),,drop=FALSE]
                segs <- segs[(rownames(segs) %in% rownames(nmat_chr)),,drop=FALSE] #Filter to chromosome being plotted
                start_indexes <- which(rownames(nmat_chr) %in% segs[,"Genome_start_interval"])
                end_indexes   <- which(rownames(nmat_chr) %in% segs[,"Genome_end_interval"])
                if(length(start_indexes) > 0) {
                    segs <- cbind(segs, types[match(segs[,"CNV"], types[,"type"]),,drop=FALSE])
                    for(i in seq(1:length(start_indexes))) {
                        indexes <- seq(start_indexes[i], end_indexes[i])
                        #Plot colored points
                        points(indexes, nmat_chr[indexes], col=segs[i, "col"])

                        #Label called CNVs with identifier
                        if(label_called_CNVs == TRUE) {
                        segments(start_indexes[i] - 0.25,
                                 as.numeric(segs[i, "ypos.line"]),
                                 end_indexes[i] + 0.25,
                                 as.numeric(segs[i, "ypos.line"]),
                                 col="orange", lwd=5)
                            text(mean(c(start_indexes[i], end_indexes[i])),
                                 as.numeric(segs[i, "ypos.label"]),
                                 segs[i, "CNV_id"],
                                 col="orange")
                        }
                    }
                }
            }
        }
	dev.off()
    }
}

rescale_chrX <- function(nmat_badX, ylimits, iqr_multiplier, out_dir) {
    #Isolate chromosome X, remove outlier probes    
    nmatX <- nmat_badX[grep("^X",rownames(nmat_badX)),]
    if(dim(nmatX)[1] == 0) {
        return(nmat_badX) # return original matrix if no probes are found on X chromosome
    }

    #Assign cases to clusters and scale medians to zero
    #    DEPRECATED METHOD: kmeans is mislead by outliers. multiple solutions possible.    
    #    clusters <- kmeans(t(nmatX),2, nstart=1)
    #    cluster_med <- apply(clusters$centers, 1, median)
    #    scaling_factors <- cluster_med[clusters$cluster]

    #CURRENT METHOD: "Partitioning Around Medoids" followed by calculation of cluster medians
    clusters <- pam(t(nmatX), 2, cluster.only=TRUE)
    cluster_med <- c()
    for(i in 1:max(clusters)) {
        cluster_med[i] <- median(nmatX[,clusters == i])
    }
    scaling_factors <- cluster_med[clusters]
    nmatX_scaled <- sweep(nmatX, 2, scaling_factors, "-")  ## subtract cluster median

    #Plot pre- and post-scaled values, write to file
    bplotX <- boxplot_cnv_matrix(nmatX, "QC_chrX_pre-scale", out_dir, ylimits, iqr_multiplier)
    bplotX_scaled <- boxplot_cnv_matrix(nmatX_scaled, "QC_chrX_post-scale", out_dir, ylimits, iqr_multiplier)

    #Infer sexes from two clusters of fractional coverage
    farthest_from_zero <- abs(cluster_med) == max(abs(cluster_med))
    if(cluster_med[farthest_from_zero] > 0) {
        cluster_female <- farthest_from_zero
    } else {
        cluster_female <- !farthest_from_zero
    }
    sexes <- sapply(cluster_female[clusters], function(x) if(x){"Female"} else {"Male"})
    names(sexes) <- names(clusters)
    write.table(sexes, paste(out_dir, "sexes.xls", sep="/") , quote = FALSE, sep = "\t", col.names = FALSE)

    #Overwrite unscaled X-chromosome values with new, scaled values
    nmat <- nmat_badX
    nmat[grep("^X",rownames(nmat)),] <- nmatX_scaled
    
    return(nmat)
}

remove_failed_samples <- function(mat, cnv_bplot_and_calls) {
    passes <- score_boxplot(cnv_bplot_and_calls)
    mat.trimmed <- mat[,passes]
    return(mat.trimmed)
}

###########
#Arguments#
###########

#Argument collection and parsing
arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) == 1) {  
    if(dev_dir != FALSE) {
        #Skips prompts if dev_dir is set
        lane_dir <- dev_dir
        out_dir <- dev_dir
    } else {
        #Collect input and output information from user

        #Input directory
        lane_dir <- choose.dir(caption = "Select a lane directory (e.g. L001):", default = lane_dir)
        if(is.na(lane_dir)) {
            try( winDialog.nonint(type="ok", "Run canceled. No input lane directory provided."), silent=TRUE)
            q(save="no")
        }
        
        #Output diretory: Attempt to derive batch information from file name, prompt user if unsuccessful
        file1 <- list.files(lane_dir, full.names=TRUE, pattern=cov_file_pattern, recursive=TRUE)[1]
        batch.regex <- "__(B*[0-9]+)"
        batch.match <- regexec(batch.regex, file1)[[1]]
        batch <- substring(file1, batch.match[2], batch.match[2] + attr(batch.match, "match.length")[2] - 1)
        if((infer.batch.for.sub.out_dir == FALSE) || is.na(batch)) {
            out_dir <- choose.dir(caption = "Select an output directory:", default = out_dir)
            batch   <- basename(out_dir)
        } else {
            out_dir <- paste(out_dir, batch, sep="/")
        }
        if(is.na(out_dir)) {
            try( winDialog.nonint(type="ok", "Run canceled. No output directory provided."), silent=TRUE)
            q(save="no")
        }

        #If output directory already exists, prompt user to overwrite
        if((clobber.output.directory == FALSE) & (file.exists(out_dir))) {
            overwrite <- try( winDialog.nonint(type="yesno", "Output directory already exists. Overwrite?"), silent=TRUE)
            if(overwrite == "NO") {
                shell(paste(explorer_file, out_dir, sep=" "), wait=FALSE)
                q(save="no")
            }
        }
		
		interval_list_dir <-arguments[1]
	}
} else if(length(arguments) == 3) {
    #Uses provided command line arguments
    lane_dir <- arguments[1]
    out_dir  <- arguments[2]
    viscap.cfg <- arguments[3]
    batch    <- basename(out_dir)
} else {
    #Usage statement
    try( winDialog.nonint(type="ok", "Usage: VisCap.R lane_directory output_directory interval_lists_directory"), silent=TRUE)
    q(save="no")
}

######
#Main#
######

# Read coverage tables
mat.all <- make_matrix_from_cov_files(lane_dir, cov_file_pattern, cov_field)
mat.all[which(mat.all==0)]=0.00001

# Read interval name files
interval_lookup <- load_interval_names(interval_list_dir, interval_file_pattern)

# Sort matrix by genome coordinates found in rownames
chroms <- c(paste("chr",1:22,sep=""),"chrX","chrY", "chrMT", "chrM")}else{chroms <- c(1:22,"X","Y", "MT", "M")}
chroms <- factor(chroms, levels=chroms, labels=chroms, ordered=TRUE)
coords <- matrix(unlist(strsplit(rownames(mat.all), ":|-")), ncol=3, byrow=TRUE, dimnames=list(rownames(mat.all)))
coords <- coords[order(match(coords[,1], chroms), as.numeric(coords[,2]), as.numeric(coords[,3])),]
mat <- mat.all[rownames(coords),]

#Iteratively run VisCap algorithm, removing bad samples after each run
if(iterative.calling.limit == 0) {
    iterative.calling.limit <- dim(mat)[2]
}

for(iteration in 1:iterative.calling.limit) {
    if(iterative.calling.limit == 1) {
        out_dir.iteration <- out_dir
    } else {
        out_dir.iteration <- paste(out_dir, "/", batch, "_run", iteration, sep="")
    }
    dir.create(out_dir.iteration, showWarnings = FALSE, recursive=TRUE)
    
    # Normalize exon coverage by exon
    nmat_badX <- log2(divide_by_batch_median(mat))
    nmat <- rescale_chrX(nmat_badX, ylimits, iqr_multiplier, out_dir.iteration)

    # Call cnvs then plot heatmaps by chromosome and per-sample exon coverage
    heatmap_by_chrom(nmat, "QC_cnv_heatmap", ylimits, out_dir.iteration)
    cnv_bplot_and_calls <- call_cnvs(nmat, ylimits, interval_lookup, threshold.min_exons, iqr_multiplier, threshold.cnv_log2_cutoffs, out_dir.iteration)
    exon_plot_per_sample(nmat, ylimits, interval_lookup, cnv_bplot_and_calls, out_dir.iteration)

    # Write out matrix of log2 ratios to file
    outfile <- paste(out_dir.iteration, "/", "log2_ratio_table", ".xls", sep="")
    gene_exon <- annotate_interval_names(rownames(nmat), interval_lookup)
    nmat.with.gene_exon <- cbind(gene_exon, nmat)
    write.table(nmat.with.gene_exon, outfile, , quote = FALSE, sep = "\t", col.names = NA)

    #Write out VisCap run information
    run_info_table <- matrix(ncol = 2, byrow=TRUE, data = c(
            "Date",                                         date(),
            "VisCap command",                               paste(commandArgs(), collapse=" "),
            "Subversion revision information",              svn.revision,
            "Batch directory",                              lane_dir,      
            "Coverage file pattern",                        cov_file_pattern,
            "Field within coverage file",                   cov_field,              
            "Output directory",                             out_dir.iteration,                
            "Interval name lookup files",                   interval_list_dir,
            "Interval name lookup file pattern",            interval_file_pattern,
            "Exons used for CNV detection",                 dim(nmat)[1],
            "Samples used for CNV detection",               dim(nmat)[2],
            "Samples not used for CNV detection",           paste(c(colnames(mat.all)[!(colnames(mat.all) %in% colnames(nmat))], ""), sep=","),
            "Plot y-axis limits",                           paste(ylimits, collapse=","),
            "Minimum consecutive exons to call CNV",        threshold.min_exons,
            "IQR multiplier used for boxplots",             iqr_multiplier,
            "Static log2 ratio thresholds to call CNVs",    paste(threshold.cnv_log2_cutoffs, collapse=","),
            "Iteration",                                    iteration
            ))
    run_info_outfile <- paste(out_dir.iteration, "/", "VisCap_run_info", ".xls", sep="")
    write.table(run_info_table, run_info_outfile, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    save.image(file=paste(out_dir.iteration, "/", "session", ".Rdata", sep=""))

    #Remove failed samples from matrix for next run
    passes <- cnv_bplot_and_calls[[1]]$names[cnv_bplot_and_calls[[1]]$qc == "PASS"]
    if(length(passes) == dim(mat)[2]) {
        break
    } else {
        #Restrict mat only to samples that pass boxplot qc
        mat <- mat[,passes]
    }
}

# Open output directory and quit R
if((dev_dir == FALSE) && (length(arguments) == 0)) {
    shell(paste(explorer_file, out_dir, sep=" "), wait=FALSE)
}
quit(save="no")
