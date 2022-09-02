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
