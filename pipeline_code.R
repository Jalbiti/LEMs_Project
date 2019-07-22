library(SingleCellExperiment)
library(nem)
library(graph)
library(Rgraphviz)
library(BiocSingular)
library(scater)
library(scran)

base_dir <- "~/Desktop/Gasperini Data"
setwd(base_dir)
input_dir <- file.path(base_dir)

# Data preparation------------------------------------------------

get_data <- function(){
  screen <- readRDS(file.path(input_dir, "GSE120861_pilot_highmoi_screen.cds.rds"))
  colmetadata <- pData(screen)
  rowmetadata <- fData(screen)
  countsdata <- exprs(screen)
  
  map <- as.data.frame(screen@phenoData@data[,15:1575])
  map <- as.data.frame(lapply(map, as.numeric)) #true input if gRNA is knocked out (true means perturbed)
  rownames(map) <- rownames(screen@phenoData@data)
  
  return(list(colmetadata = colmetadata, rowmetadata = rowmetadata,
              countsdata = countsdata, map = map))
}

## create sce object and sample
create_sce <- function(colmetadata, rowmetadata, countsdata) {
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=as.matrix(countsdata)), colData=colmetadata, rowData=rowmetadata)
  set.seed(1234)
  subsample <- sample(ncol(sce), 5000)
  sce <- sce[,subsample]
}

# Data processing-------------------------------------------------

#' Do quality control on a SingleCellExperiment
#'
#' @param ps A SingleCellExperiment object
#' @return SingleCellExperiment object after quality control trimming
#' @export
qualitycontrol <- function(ps, plots=TRUE) {
  # get rid of rows that have never been observed
  keep_feature <- rowSums(DESeq2::counts(ps)) > 0
  summary(keep_feature) # all features are observed
  # but some are NA's
  ps <- ps[! is.na(keep_feature), ]
  
  # tag mitochondrial genes
  mt_genes <- grepl("MT", SummarizedExperiment::rowData(ps)$chromosome_name)
  ps <- scater::calculateQCMetrics(ps, feature_controls = list(MT = mt_genes), percent_top = 100)
  colnames(SummarizedExperiment::colData(ps))
  
  # experiment had no spike ins, so don't need to mark those
  
  if(plots) {
    p1 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = log10_total_counts_endogenous, y = 1, fill = 1)) +
      ggridges::geom_density_ridges(scale = 4) + 
      ggplot2::ggtitle(expression(log[10]*" library size")) + 
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::scale_y_discrete(expand = c(0.01, 0)) + # will generally have to set the `expand` option
      ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1)) 
    
    # total detected genes
    p2 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = total_features_by_counts, y = 1,
                                                                                         fill = 1)) +
      ggridges::geom_density_ridges(scale = 4) + 
      ggplot2::ggtitle("total features by counts detected") +
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::scale_y_discrete(expand = c(0.01, 0)) +
      ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1))
    
    # percentage reads in mitochondrial genes
    p3 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = pct_counts_MT, y = 1,
                                                                                         fill = 1)) +
      ggridges::geom_density_ridges(scale = 4) +
      ggplot2::ggtitle("% reads in mitochondrial genes") +
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::scale_y_discrete(expand = c(0.01, 0)) +
      ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1))
    
    # percentage of counts taken by the top 100 expressed features
    p4 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = pct_counts_in_top_100_features, y = 1,
                                                                                         fill = 1)) +
      ggridges::geom_density_ridges(scale = 4) + 
      ggplot2::ggtitle("% counts from top 100 expressed genes") +
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::scale_y_discrete(expand = c(0.01, 0)) +
      ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1))
    
    print(scater::multiplot(plotlist = list(p1, p2, p3, p4), cols = 2))
  }
  
  # QC data by removing bad quality cells by removing outliers
  # the data looks like it might have already been qc'ed?
  bad_qual_mt <- scater::isOutlier(ps$pct_counts_MT, nmads=3, type="both")
  bad_qual_tf <- scater::isOutlier(ps$total_features_by_counts, nmads=3, type="both")
  bad_qual_tc <- scater::isOutlier(ps$log10_total_counts_endogenous, nmads=3, type="both")
  bad_qual_pc <- scater::isOutlier(ps$pct_counts_in_top_100_features, nmads=3, type="both")
  ps_GQ <- ps[,!(bad_qual_mt | bad_qual_tf | bad_qual_tc | bad_qual_pc)]
  # do we lose any features because cells have been removed?
  ps_GQ <- ps_GQ[rowMeans(DESeq2::counts(ps_GQ))>0,] # no op
  
  return(ps_GQ)
  
}

sce_processing <- function(sce) {
  #Before clustering we perform normalization, feature selection and dimensionality reduction
  # do quality control and normalize data
  sce <- qualitycontrol(sce)
  sce <- normalizeSCE(sce)
  
  ## feature selection
  fit <- scran::trendVar(sce, use.spikes = FALSE)
  dec <- scran::decomposeVar(sce, fit)
  hvg <- rownames(dec[dec$bio > 0, ])
  
  ## dimensionality reduction
  sce <- runPCA(sce, ncomponents = 30,
                feature_set = hvg)
  
  sce <- runUMAP(sce, use_dimred = 'PCA', n_dimred = 20)
  
  return(sce)
}

## Clustering
library(igraph)

## replace `use.dimred` with 'PCA' if no integration was performed
## this will be automatically added via `runUMAP` above

sce_clustering <- function(sce) {
  set.seed(1234)
  g <- buildSNNGraph(sce,
                     # PCA = 'MNN',
                     k = 30,  # higher = bigger clusters
                     BNPARAM = BiocNeighbors::AnnoyParam(),
                     BPPARAM = BiocParallel::MulticoreParam())
  clusters <- as.factor(igraph::cluster_louvain(g)$membership)
  sce$clusters <- clusters
  
  return(sce)
}

# Subsetting------------------------------------------------------

## subset enhancers (to identified ones for now, then hamming)
identified_enh <- function(map) {
  eg_pairings <- read.csv(file.path(input_dir, "145_enhancer_gene.csv"), header = TRUE)
  eg_pairings$enhancer <- paste0(eg_pairings$chr.candidate_enhancer, ".",
                                 eg_pairings$start.candidate_enhancer)
  
  reduced_map <- data.frame(row.names = rownames(map))
  #reduce map to candidate enhancer
  map <- map[, which((grepl("^chr", colnames(map))))]
  enhancers <- unique(paste0(eg_pairings$chr.candidate_enhancer, ".",
                             eg_pairings$start.candidate_enhancer)) #as we observed that candidate enhancers only matched the starting column not the stopping
  
  for (i in 1:ncol(map)) {
    column <- colnames(map)[i]
    for (j in enhancers) {
      if(grepl(j, column)) { #to limit columns to enhancer-gene identified pairs
        reduced_map <- cbind(reduced_map, map[,i])
        colnames(reduced_map)[ncol(reduced_map)] <- j #rename columns to identified enhancer
      }
    }
  }
  
  #check that we have at least one perturbation per each enhancer
  stopifnot(colSums(reduced_map) >= 1)
  
  #keep only cells where at least one perturbation has occurred
  reduced_map <- reduced_map[-which(rowSums(reduced_map) < 1),]
  
  return(reduced_map)
}


sce_subsetting <- function(sce, coldata = colmetadata,
                           countsData = countsdata,
                           rowmetaData = rowmetadata,
                           subset_cells = TRUE,
                           subset_genes = TRUE, cluster_spec = FALSE,
                           pathway, boot_cells = FALSE, boot_genes = FALSE) {
  if (subset_cells) {
    if (sum(as.numeric(sce$clusters)) == 0) {
      return("Need to cluster before subsetting.")
    } else {
      n_clusters <- length(unique(sce$clusters))
      if (!(cluster_spec)) {
        #select a random cluster_spec
        if (!(boot_cells)) {
          set.seed(1234)
        }
        cluster_spec <- sample(c(1:n_clusters), 1)
        sce <- sce[,which(sce$clusters == cluster_spec)]
        colmetadata <- coldata[which(sce$clusters == cluster_spec), ]
      } else {
        #get proportions of clusters
        prop <- c()
        for (i in 1:n_clusters) {
          denom <- length(which(sce$clusters == i))
          prop[i] <- denom / length(sce$clusters)
        }
        #subsample to 5000 cells according to clusters
        total <- 5000
        indices <- c()
        for (i in 1:n_clusters) {
          set.seed(1234)
          subsample <- round(total * prop[i])
          indices <- append(indices,
                            sample(which(sce$clusters == i), subsample))
        }
        indices <- unique(indices)
        sce <- sce[, indices]
      }
    }
  }
  if (subset_genes) {
    genes <- DelayedMatrixStats::rowVars(logcounts(sce))
    names(genes) <- rownames(sce)
    genes <- sort(genes, decreasing = TRUE)
    
    #Specify the top 1000 most highly variable genes (hvg) by name #for future use
    metadata(sce)$hvg_genes <- names(genes)[1:1000]
    # hvg_names <- rowData(sce)$gene_short_name[which(rowData(sce)$id %in% names(genes)[1:1000])]
    
    #Select specific pathway from pathway enrichment
    screen <- readRDS(file.path(input_dir, "GSE120861_pilot_highmoi_screen.cds.rds"))
    all_genes <- fData(screen)$gene_short_name
    selected_genes <- unique(pathway[which(pathway[,1] %in% all_genes), 1])
    sce <- sce[which(rowData(sce)$gene_short_name %in% selected_genes),]
    
    #Look for paired gene-enhancers
    eg_pairings <- read.csv(file.path(input_dir, "145_enhancer_gene.csv"), header = TRUE)
    eg_pairings$enhancer <- paste0(eg_pairings$chr.candidate_enhancer, ".",
                                   eg_pairings$start.candidate_enhancer)
    paired_enh <- eg_pairings$enhancer[which(eg_pairings$target_gene_short %in% selected_genes)]
    paired_genes <- data.frame(enhancer = paired_enh)
    paired_genes$gene <- eg_pairings$target_gene_short[match(paired_enh, eg_pairings$enhancer)]
    
    #SINGLE PERTURB
    map <- as.data.frame(screen@phenoData@data[,15:1575])
    map <- as.data.frame(lapply(map, as.numeric)) #true input if gRNA is knocked out (true means perturbed)
    rownames(map) <- rownames(screen@phenoData@data)
    single_pb <- FALSE
    if (single_pb) {
      colmetadata <- coldata[which(rownames(coldata) %in% rownames(single_perturb)),]
      countsdata <- countsdata[which(rownames(countsdata) %in% genes),
                               which(colnames(countsdata) %in% rownames(single_perturb))]
    }
    
    single_pb_E <- FALSE
    if (single_pb_E) {
      #Subset cells to single perturbed
      single_perturb <- map[which(rowSums(map) == 1),]
      sp_cell_enh <- rownames(map)[which(rowSums(map) == 1)]
      single_perturb <- single_perturb[, which((grepl("^chr", colnames(single_perturb))))]
      
      non_perturb <- colnames(single_perturb[, which(colSums(single_perturb) == 0)])
      single_perturb <- single_perturb[, which(colSums(single_perturb) != 0)] #keep enhancers that have been perturbed at least once
      
      E <- single_perturb
      E <- identified_enh(E)
      
      E <- E[, which(colnames(E) %in% paired_enh)]
    } else {
      non_single_perturb <- map
      non_single_perturb <- non_single_perturb[, which((grepl("^chr", colnames(non_single_perturb))))]
      non_single_perturb <- non_single_perturb[, which(colSums(non_single_perturb) != 0)] #keep enhancers that have been perturbed at least once
      
      E <- non_single_perturb
      E <- identified_enh(E)
      
      E <- E[, which(colnames(E) %in% paired_enh)]
    }
    
    Y <- t(countsData)
    colnames(Y) <- rowmetaData$gene_short_name[which(rowmetaData$id == colnames(Y))]
    Y <- Y[which(rownames(Y) %in% rownames(E)), which(colnames(Y) %in% selected_genes)]
    if (boot_genes) {
      indices <- sample(c(1:ncol(Y)), round(ncol(Y) * 0.5))
      Y <- Y[,indices]
    }
    Y <- Y[which(rowSums(Y) != 0), which(colSums(Y) != 0)]
    E <- E[which(rownames(E) %in% rownames(Y)),] #double check
    
    return(list(sce = sce, Y = Y, E = E))
  }
}

# LEM-------------------------------------------------------------

final_lem <- function(Y, E) {
  sapply(list.files("./R_lem/", full.names = TRUE), source)
  set.seed(1234)
  lem_res <- lem(as.matrix(Y), as.matrix(E), inference = "greedy",
                 parameter.estimation = "linear.reg", verbose = TRUE)
  
  return(lem_res)
}

# print.lem(lem_res)
# plot.lem(lem_res)
# plot(lem_res$all.scores)


# Bootstrapping---------------------------------------------------

## we will check the network's stability by bootstrapping

bootstrapping_cells <- function() {
  data_sce <- get_data()
  print("Got Data")
  countsdata <- data_sce$countsdata
  colmetadata <- data_sce$colmetadata
  rowmetadata <- data_sce$rowmetadata
  
  sce <- create_sce(countsdata = countsdata,
                    colmetadata = colmetadata,
                    rowmetadata = rowmetadata)
  print("Created SCE")
  
  sce <- sce_processing(sce)
  print("Processed SCE")
  sce <- sce_clustering(sce)
  print("Clustered SCE")
  pathway_GO0050900 <- read.csv("./GO_0050900.txt", header = FALSE)
  data_sce <- sce_subsetting(sce, coldata = colmetadata,
                             pathway = pathway_GO0050900,
                             boot_cells = TRUE)
  print("Subsetted SCE")
  # plotUMAP(sce, colour_by = "clusters")
  
  sce <- data_sce[[1]]
  Y <- data_sce[[2]]
  E <- data_sce[[3]]
  
  lem_res <- final_lem(Y, E)
  return(lem_res)
}

boot_cells_lem <- list()
for (i in 1:5) {
  boot_cells_lem[[i]] <- bootstrapping_cells()
  print(i/5, "boots of cells")
}

bootstrapping_genes <- function() {
  data_sce <- get_data()
  
  countsdata <- data_sce$countsdata
  colmetadata <- data_sce$colmetadata
  rowmetadata <- data_sce$rowmetadata
  
  sce <- create_sce(countsdata = countsdata,
                    colmetadata = colmetadata,
                    rowmetadata = rowmetadata)
  
  sce <- sce_processing(sce)
  sce <- sce_clustering(sce)
  pathway_GO0050900 <- read.csv("./GO_0050900.txt", header = FALSE)
  data_sce <- sce_subsetting(sce, pathway = pathway_GO0050900,
                             boot_genes = TRUE)
  
  # plotUMAP(sce, colour_by = "clusters")
  
  sce <- data_sce[[1]]
  Y <- data_sce[[2]]
  E <- data_sce[[3]]
  
  lem_res <- final_lem(Y, E)
  return(lem_res)
}

boot_genes_lem <- list()
for (i in 1:5) {
  boot_genes_lem[[i]] <- bootstrapping_genes()
  print(i/5, "boots of genes")
}

