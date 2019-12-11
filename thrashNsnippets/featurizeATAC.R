################################################################
### Functions for gene-level featurization of ATAC-seq peaks ###
################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Organism.dplyr)
  library(pbapply)
  library(future)
  library(reshape2)
  library(cicero)
})

### Gene-level featurization ###

#' Create Gene Activity Matrix (Seurat style)
#' Creates matrix of gene activity by counting peaks within
#' gene bodies and promoters
#' 
#' @param peak.matrix Matrix of peaks (as read by filtered_peak_bc_matrix output from cellranger)
#' @param annotation.file gtf file from ensembl with gene annotations (transformed to conform with UCSC 
#' notation with chr)
#' @param seq.levels chromosome names 
#' @param include.body logical indicating whether gene bodies should be included
#' @param upstream no. of base pairs to add upstream
#' @param downstream no. of base pairs to add downstream
#' @param verbose
#' 
#' @return dcgMatrix object of gene activities
#' 
makeGeneActivity_Seurat <- function( peak.matrix,
                                        annotation.file,
                                        seq.levels = c(1:22, "X", "Y"),
                                        include.body = TRUE,
                                        upstream = 2000,
                                        downstream = 0,
                                        verbose = TRUE){
  
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  
  # if any peaks start at 0, change to 1
  # otherwise GenomicRanges::distanceToNearest will not work 
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
  
  # get annotation file, select genes
  gtf <- rtracklayer::import(con = annotation.file)
  
  # change seqlevelsStyle if not the same
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  gtf.genes <- gtf[gtf$type == 'gene']
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse') ## I basically just swiched position of this line
  
  # Extend definition up/downstream
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  
  # Some GTF rows will not have gene_name attribute
  # Replace it by gene_id attribute
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
  
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c('peak', 'gene.name')]
  colnames(x = annotations) <- c('feature', 'new_feature')
  
  # collapse into expression matrix
  peak.matrix <- as(object = peak.matrix, Class = 'matrix')
  all.features <- unique(x = annotations$new_feature)
  
  if (nbrOfWorkers() > 1) {
    mysapply <- future_sapply
  } else {
    mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  }
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){
    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      return(Matrix::colSums(x = submat))
    } else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = 'dgCMatrix'))
}

matrixToCDS <- function(peak.matrix){
  ## Check for right separator in peak annotation
  if (length(strsplit(rownames(peak.matrix)[1], '-')[[1]])){
    rownames(peak.matrix) <- str_replace(rownames(peak.matrix), ":", "_") %>%
      str_replace("-", "_")
  }
  long.peaks <- melt(as.matrix(peak.matrix)) 
  atac_cds <- cicero::make_atac_cds(long.peaks)
  return(atac_cds)
}


# # load in your data using rtracklayer
# # annotation_file <- "~/annotations/Homo_sapiens.GRCh38.86.gtf"
# gene_anno <- rtracklayer::readGFF(annotation_file)
# 
# # rename some columns to match plotting requirements
# gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
# gene_anno$gene <- gene_anno$gene_id
# gene_anno$transcript <- gene_anno$transcript_id
# gene_anno$symbol <- gene_anno$gene_name
# 
# #### Add a column for the pData table indicating the gene if a peak is a promoter ####
# pos <- subset(gene_anno, strand == "+")
# pos <- pos[order(pos$start),] 
# # remove all but the first exons per transcript
# pos <- pos[!duplicated(pos$transcript),] 
# # make a 1 base pair marker of the TSS
# pos$end <- pos$start + 1 
# 
# neg <- subset(gene_anno, strand == "-")
# neg <- neg[order(neg$start, decreasing = TRUE),] 
# # remove all but the first exons per transcript
# neg <- neg[!duplicated(neg$transcript),] 
# neg$start <- neg$end - 1
# 
# gene_annotation_sub <- rbind(pos, neg)
# 
# # Make a subset of the TSS annotation columns containing just the coordinates 
# # and the gene name
# gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]
# 
# # Rename the gene symbol column to "gene"
# names(gene_annotation_sub)[4] <- "gene"
# 
# ### Annotate CDS by promoters ###
# atac_cds <- annotate_cds_by_site(atac_cds, gene_annotation_sub)
# 
# ### Run cicero co-accessibility scores ###
# hg38.genome <- read.table("~/annotations/hg38.genome")
# sample_genome <- subset(hg38.genome, V1 == "chr1")
# sample_genome$V2[1] <- 10000000
# # Estimate distance parameter
# dist_param <- estimate_distance_parameter(atac_cds, genomic_coords = sample_genome)
# 
# 
# #### Generate gene activity scores ####
# # generate unnormalized gene activity matrix
# unnorm_ga <- build_gene_activity_matrix(atac_cds, conns)
# 
# # remove any rows/columns with all zeroes
# unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
#                        !Matrix::colSums(unnorm_ga) == 0]
# 
# # make a list of num_genes_expressed
# num_genes <- pData(atac_cds)$num_genes_expressed
# names(num_genes) <- row.names(pData(atac_cds))
# 
# # normalize
# cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
# 
# # if you had two datasets to normalize, you would pass both:
# # num_genes should then include all cells from both sets
# unnorm_ga2 <- unnorm_ga
# cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), 
#                                                     num_genes)
# 
# filtered_peak_h5_file = "~/my_data/filtered_peak_bc_matrix.h5"
# annotation.file <- "~/annotations/Homo_sapiens.GRCh38.86.gtf"
# 

makeGeneActivity <- function(filtered_peak_h5_file, annotation_file, type="Seurat", save=FALSE, sample_name=NULL){
  peaks <- Read10X_h5(filtered_peak_h5_file)
  chromosomes = c(1:22, "X", "Y")
  chromosomes = paste("chr", chromosomes, sep="")
  if (type=="Seurat") {
    myMakeGeneActivity <- makeGeneActivity_Seurat
  }
  countMat <- myMakeGeneActivity(peak.matrix = peaks, annotation.file = annotation.file, seq.levels = chromosomes, upstream = 2000)
  # if (save) {
  #   if (is.null(sample_name)) {
  #     warning()
  #   }
  #   outfile_name <- str_remove(filtered_peak_h5_file, 'filtered_peak_bc.matrix.h5') %>% str_c(sample_name, ".geneActivity.tsv")
  #   write.table(countMat, file=)
  # }
  return(countMat)
  }

### Utility functions

# Resize GenomicRanges upstream and or downstream
# from https://support.bioconductor.org/p/78652/
# (Seurat function but not exported)
Extend <- function(x, upstream = 0, downstream = 0) {
  if (any(GenomicRanges::strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
  new_start <- GenomicRanges::start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
  new_end <- GenomicRanges::end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- GenomicRanges::trim(x = x)
  return(x)
}