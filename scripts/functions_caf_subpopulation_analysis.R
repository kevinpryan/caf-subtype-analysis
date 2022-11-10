#!/usr/bin/Rscript

# a set of functions called from the main Rmd documents in this repo.

get_upregulated <- function(df){
  
  key <- intersect(rownames(df)[which(df$log2FoldChange>=2)], rownames(df)[which(df$padj<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

anti_join_rownames <- function(df1, df2){
  anti_join((df1 %>% mutate(Symbol = rownames(df1))),
            (df2 %>% mutate(Symbol = rownames(df2))),
            by = 'Symbol')
}

filter_dfs_antijoin_rownames <- function(df_list){
  output <- list()
  for (i in 1:length(dfs_to_filter)){
    df_interest <- dfs_to_filter[[i]]
    if (nrow(df_interest) == 0) stop("one of your dataframes has no rows")
    not_i <- seq(1,length(dfs_to_filter))[seq(1,length(dfs_to_filter)) != i]
    for (j in not_i){
      df_filtered <- anti_join_rownames(df_interest, dfs_to_filter[[j]])
      if (nrow(df_filtered) == 0){
        break
      }
    }
    df_filtered <- subset(df_filtered, select = -c(Symbol)) 
    output <- c(output, list(df_filtered))
    print(dim(output)[[i]])
  }
  return(output)
}

annotate_de_genes <- function(df, filter_by){
  # if your df has hgnc_symbol as rownames, filter by that, if it is the ENSG1234.12, use "ensembl_gene_id_version", if it is the regular engs, filter by "ensembl_gene_id"
  filter_by_string <- as.character(filter_by)
  df$gene_symbol <- rownames(df)
  colnames(df)[6] <- filter_by_string
  #print(df)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",  #host="https://www.ensembl.org")
                  host="uswest.ensembl.org")
  info <- getBM(attributes=c("hgnc_symbol",
                             "ensembl_gene_id_version",
                             "chromosome_name",
                             "start_position",
                             "end_position",
                             "strand",
                             "entrezgene_description",
                             "entrezgene_id"),
                filters = c(filter_by_string),
                values = df[,6],
                mart = mart,
                useCache=FALSE)
  
  tmp <- merge(df, info, by=filter_by_string)
  tmp$strand <- gsub("-1", "-", tmp$strand)
  tmp$strand <- gsub("1", "+", tmp$strand)
  #tmp$hgnc_symbol <- make.names(tmp$hgnc_symbol, unique = T)
  tmp <- tmp[!grepl("CHR", tmp$chromosome_name),]
  
  output_col <- c("Gene", "Ensembl ID", "Chromosome", "Start", "Stop", "Strand", "Description", "Log2FC", "P-value", "Adj P-value", "Entrez ID")
  tmp <- subset(tmp, select=c(hgnc_symbol, ensembl_gene_id_version, chromosome_name, start_position, end_position, strand, entrezgene_description, log2FoldChange, pvalue, padj, entrezgene_id))
  colnames(tmp) <- output_col
  
  if(min(tmp$Log2FC) > 0){
    tmp <- tmp[order(-tmp$Log2FC),]
  }else{
    tmp <- tmp[order(tmp$Log2FC),]
  }
  
  return(tmp)
  
}

create_df_plotcounts <- function(dds, gene_interest, info, symbol_type_dds = "hgnc_symbol"){
  library(DESeq2)
  library(biomaRt)
  # put check for type of gene symbol in dds object, assuming ensembl gene version atm
  # assume info have column called hgnc_symbol
  info = info[!duplicated(info$hgnc_symbol),]
  rownames(info) = info$hgnc_symbol
  # assuming genes interest are in info, put in check for this to deal with missing gene symbol
  if (!(gene_interest %in% info$hgnc_symbol)){
    print(paste(gene_interest, "not found in table of hgnc symbols", sep = " "))
    geneCounts <- NA
  } else {
    gene_interest_info <- info[gene_interest,]
    #print(gene_interest_info)
    if (symbol_type_dds == "ensembl_gene_id_version"){
      gene_ensg <- info$ensembl_gene_id[which(info$hgnc_symbol == gene_interest)]
      # find ENSG1234.1 in dds object given ENGS1234
      idx <- str_detect(rownames(assay(dds)), paste0(gene_ensg, "\\.."))
      gene_correct_version <- rownames(assay(dds))[idx]
    } else if (symbol_type_dds == "hgnc_symbol"){
      #idx <- str_detect(rownames(assay(dds)), gene_interest)
      #gene_correct_version <- rownames(assay(dds))[idx]
      gene_correct_version <- gene_interest
    } else {
      stop("symbol_type_dds must be one of hgnc_symbol or ensembl_gene_id_version")
    }
    geneCounts <- plotCounts(dds, gene = gene_correct_version, intgroup = c("Study", "Subpopulation", "Tumor_JuxtaTumor"), returnData = TRUE)
  }
}


caf_plot_tumour_juxtatumour <- function(df_for_plotting, gene = NULL){
  library(ggplot2)
  gene = gene
  if (is.list(df_for_plotting)){
    cols <- colnames(df_for_plotting)
    df_for_plotting <- as.data.frame(df_for_plotting)
    colnames(df_for_plotting) <- cols
  }
  ggplot(df_for_plotting, aes(x = Subpopulation, y = count,  colour = Tumor_JuxtaTumor)) +
    geom_point(size = 1,  # reduce point size to minimize overplotting 
               position = position_jitter(
                 width = 0.15,  # amount of jitter in horizontal direction
                 height = 0     # amount of jitter in vertical direction (0 = none)
               )
    ) +
    scale_y_log10() +
    labs(title = gene) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ylab("Normalised counts")  #, axis.title.x = element_blank())
}

caf_plot_study <- function(df_for_plotting, gene = NULL){
  library(ggplot2)
  gene = gene
  if (is.list(df_for_plotting)){
    cols <- colnames(df_for_plotting)
    df_for_plotting <- as.data.frame(df_for_plotting)
    colnames(df_for_plotting) <- cols
  }
  ggplot(df_for_plotting, aes(x = Subpopulation, y = count,  colour = Study)) +
    geom_point(size = 1,  # reduce point size to minimize overplotting 
               position = position_jitter(
                 width = 0.15,  # amount of jitter in horizontal direction
                 height = 0     # amount of jitter in vertical direction (0 = none)
               )
    ) +
    scale_y_log10() +
    labs(title = gene) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    ylab("Normalised counts")  #, axis.title.x = element_blank())
}

