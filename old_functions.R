#!/usr/bin/Rscript

# old version of this function
caf_plot_tumour_juxtatumour <- function(annotated_df_for_plotting, gene_of_interest, label_yaxis = gene_of_interest, info){
  library(ggplot2)
  
  if (is.null(label_yaxis)){
    label_yaxis <- gene_of_interest
  }
  annotated_df_for_plotting$gene_of_interest <- annotated_df_for_plotting[, gene_of_interest]
  ggplot(annotated_df_for_plotting, aes(Subpopulation, gene_of_interest, colour = Tumor_JuxtaTumor)) +
    geom_point(size = 1,  # reduce point size to minimize overplotting 
               position = position_jitter(
                 width = 0.15,  # amount of jitter in horizontal direction
                 height = 0     # amount of jitter in vertical direction (0 = none)
               )
    ) +
    scale_y_continuous(gene_of_interest) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())  +
    ylab(label_yaxis)  #, axis.title.x = element_blank())
}

caf_plot_function_study_old <- function(input_table, gene_of_interest){
  library(ggplot2)
  input_table$gene_of_interest <- input_table[, gene_of_interest]
  ggplot(input_table, aes(Subpopulation, gene_of_interest, colour = Study)) +
    geom_point(size = 1,  # reduce point size to minimize overplotting 
               position = position_jitter(
                 width = 0.15,  # amount of jitter in horizontal direction
                 height = 0     # amount of jitter in vertical direction (0 = none)
               )
    ) +
    scale_y_continuous(gene_of_interest) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())
}

caf_plot_function_tumour_juxtatumour_old <- function(input_table, gene_of_interest){
  library(ggplot2)
  input_table$gene_of_interest <- input_table[, gene_of_interest]
  ggplot(input_table, aes(Subpopulation, gene_of_interest, colour = Tumor_JuxtaTumor)) +
    geom_point(size = 1,  # reduce point size to minimize overplotting 
               position = position_jitter(
                 width = 0.15,  # amount of jitter in horizontal direction
                 height = 0     # amount of jitter in vertical direction (0 = none)
               )
    ) +
    scale_y_continuous(gene_of_interest) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank())
}

# Function to take in dds object, ensembl - hgnc info and genes interest and output df for plotting
# if this function doesn't work, try running the above chunk again (getBM function)
create_annotated_df_for_plotting <- function(dds, info, genes_interest){
  library(stringr)
  # put check for type of gene symbol in dds object, assuming ensembl gene version atm
  # assume info have column called hgnc_symbol
  info <- info[!duplicated(info$hgnc_symbol),]
  rownames(info) <- info$hgnc_symbol
  # assuming genes interest are in info, put in check for this to deal with missing gene symbol
  genes_interest_info <- info[genes_interest,]
  # again, here assuming ENSG1234.1 etc
  ensembl_genes_dds <- str_split_fixed(rownames(assay(dds)), pattern = "\\.", n = 2)[,1]
  ids <- which(ensembl_genes_dds %in% genes_interest_info$ensembl_gene_id)
  dds_genes_interest <- dds[ids,]
  dds_genes_interest_exprs <- data.frame(assay(dds_genes_interest))
  metadata <- as.data.frame(colData(dds))
  colnames(dds_genes_interest_exprs) <- rownames(metadata)
  rownames(dds_genes_interest_exprs) <- genes_interest
  dds_genes_interest_exprs_t <- as.data.frame(t(dds_genes_interest_exprs))
  dds_genes_interest_exprs_t$Subpopulation <- metadata$Subpopulation
  dds_genes_interest_exprs_t$Study <- metadata$Study
  dds_genes_interest_exprs_t$Tumor_JuxtaTumor <- metadata$Tumor_JuxtaTumor
  return(dds_genes_interest_exprs_t)
}

create_annotated_df_for_plotting_2 <- function(dds, genes_interest){
  library(DESeq2)
  # put check for type of gene symbol in dds object, assuming ensembl gene version atm
  # assume info have column called hgnc_symbol
  info <- info[!duplicated(info$hgnc_symbol),]
  rownames(info) <- info$hgnc_symbol
  # assuming genes interest are in info, put in check for this to deal with missing gene symbol
  genes_interest_info <- info[genes_interest,]
  outputs <- list()
  for (i in 1:length(genes_interest)){
    gene_hgnc <- genes_interest[i]
    gene_ensg <- info$ensembl_gene_id[which(info$hgnc_symbol == genes_interest[i])]
    # find ENSG1234.1 in dds object given ENGS1234
    idx <- str_detect(rownames(assay(dds)), paste0(gene_ensg, "\\.."))
    gene_ensg_version <- rownames(assay(dds))[idx] 
    geneCounts <- plotCounts(dds, gene = gene_ensg_version, intgroup = c("Study", "Subpopulation", "Tumor_JuxtaTumor"), returnData = TRUE)
    outputs <- append(outputs, geneCounts)
  }
}

caf_plot_study <- function(annotated_df_for_plotting, gene_of_interest, label_yaxis = gene_of_interest, info){
  library(ggplot2)
  info <- info[!duplicated(info$hgnc_symbol),]
  rownames(info) <- info$hgnc_symbol
  # assuming genes interest are in info, put in check for this to deal with missing gene symbol
  genes_interest_info <- info[genes_interest,]
  if (is.null(label_yaxis)){
    label_yaxis <- gene_of_interest
  }
  annotated_df_for_plotting$gene_of_interest <- annotated_df_for_plotting[, gene_of_interest]
  ggplot(annotated_df_for_plotting, aes(Subpopulation, gene_of_interest, colour = Study)) +
    geom_point(size = 1,  # reduce point size to minimize overplotting 
               position = position_jitter(
                 width = 0.15,  # amount of jitter in horizontal direction
                 height = 0     # amount of jitter in vertical direction (0 = none)
               )
    ) +
    scale_y_continuous(gene_of_interest) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank()) +
    ylab(label_yaxis)#
}