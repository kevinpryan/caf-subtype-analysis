filter_out_low_expressed <- function(dds){
  library(DESeq2)
  print(paste("no of genes before filtering...", nrow(dds)))
  # returns a vector of whether the total count of each gene is >= 10 (True or false)
  keep <- rowSums(counts(dds)) >= 10
  # only keep rows (genes) for which keep is TRUE
  dds <- dds[keep,]
  # at least X samples with a count of 10 or more, where X is 5% of samples
  X <- round(0.05*ncol(dds))
  keep <- rowSums(counts(dds) >= 10) >= X
  dds <- dds[keep,]
  print(paste("no of genes after filtering...", nrow(dds)))
  return(dds)
}

# define function to remove batch effect
limma_batch_effect <- function(binary_vector, vsd, model.matrix, coldata, extra_info = ""){
      library(limma)
      library(tibble)
      if (binary_vector[1] == T){
        batch <- coldata$Study
        outfile_study <- "remove_study_"
      } else {
        batch <- NULL
        outfile_study <- ""
      }
      if (binary_vector[2] == T){
        batch2 <- coldata$Patient
        outfile_patient <- "remove_patient_"
      } else {
        batch2 <- NULL
        outfile_patient <- ""
      }
      vsd.batch.corrected <- limma::removeBatchEffect(x = vsd, 
                                                      batch = batch, 
                                                      batch2 = batch2, 
                                                      design = model.matrix
                                                     )
      outfile <- paste("~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction/", Sys.Date(), "_limma_remove_batch_effect_", outfile_study, outfile_patient, extra_info, ".txt", sep = "")
      non.log.vsd.batch.corrected <- 2^vsd.batch.corrected
      non.log.vsd.batch.corrected <- tibble::rownames_to_column(as.data.frame(non.log.vsd.batch.corrected), var = "Gene")
      colnames(non.log.vsd.batch.corrected) <- gsub("\\.", "-", colnames(non.log.vsd.batch.corrected))
      stopifnot(all(coldata$names == colnames(non.log.vsd.batch.corrected[2:ncol(non.log.vsd.batch.corrected)])))
      non.log.vsd.batch.corrected.for.cibersort <- non.log.vsd.batch.corrected
      non.log.vsd.batch.corrected.for.cibersort <- rbind(c("Gene", coldata$Subpopulation), non.log.vsd.batch.corrected.for.cibersort)
      
      s1_cols <- which(non.log.vsd.batch.corrected.for.cibersort[1,] == "S1")
      s3_cols <- which(non.log.vsd.batch.corrected.for.cibersort[1,] == "S3")
      s4_cols <- which(non.log.vsd.batch.corrected.for.cibersort[1,] == "S4")
      
      s1_cols_length <- length(s1_cols)
      s3_cols_length <- length(s3_cols)
      s4_cols_length <- length(s4_cols)
      
      # Put signature data in the right order ie all the S1s together etc
      non.log.vsd.batch.corrected.for.cibersort.order <- cbind.data.frame(non.log.vsd.batch.corrected.for.cibersort$Gene,
                                                                          non.log.vsd.batch.corrected.for.cibersort[,s1_cols],
                                                                          non.log.vsd.batch.corrected.for.cibersort[,s3_cols],
                                                                          non.log.vsd.batch.corrected.for.cibersort[,s4_cols]
      )
      colnames(non.log.vsd.batch.corrected.for.cibersort.order) <- NULL
      rownames(non.log.vsd.batch.corrected.for.cibersort.order) <- NULL
      print("outfile...")
      print(outfile)
      return(list(outfile, non.log.vsd.batch.corrected.for.cibersort.order))
}

# define function to remove batch effect
limma_batch_effect_include_inhouse <- function(binary_vector, vsd, model.matrix, coldata, extra_info = ""){
  library(limma)
  library(tibble)
  if (binary_vector[1] == T){
    batch <- coldata$Study
    outfile_study <- "remove_study_"
  } else {
    batch <- NULL
    outfile_study <- ""
  }
  if (binary_vector[2] == T){
    batch2 <- coldata$Patient
    outfile_patient <- "remove_patient_"
  } else {
    batch2 <- NULL
    outfile_patient <- ""
  }
  vsd.batch.corrected <- limma::removeBatchEffect(x = vsd, 
                                                  batch = batch, 
                                                  batch2 = batch2, 
                                                  design = model.matrix
  )
  outfile <- paste("~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction/", Sys.Date(), "_limma_remove_batch_effect_", outfile_study, outfile_patient, extra_info, ".txt", sep = "")
  outfile_inhouse <- paste("~/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction/", Sys.Date(), "_limma_remove_batch_effect_", outfile_study, outfile_patient, extra_info, "inhouse_data", ".txt", sep = "")
  non.log.vsd.batch.corrected <- 2^vsd.batch.corrected
  non.log.vsd.batch.corrected <- tibble::rownames_to_column(as.data.frame(non.log.vsd.batch.corrected), var = "Gene")
  colnames(non.log.vsd.batch.corrected) <- gsub("\\.", "-", colnames(non.log.vsd.batch.corrected))
  stopifnot(all(coldata$names == colnames(non.log.vsd.batch.corrected[2:ncol(non.log.vsd.batch.corrected)])))
  non.log.vsd.batch.corrected.for.cibersort <- non.log.vsd.batch.corrected
  non.log.vsd.batch.corrected.for.cibersort <- rbind(c("Gene", coldata$Subpopulation), non.log.vsd.batch.corrected.for.cibersort)
  
  s1_cols <- which(non.log.vsd.batch.corrected.for.cibersort[1,] == "S1")
  s3_cols <- which(non.log.vsd.batch.corrected.for.cibersort[1,] == "S3")
  s4_cols <- which(non.log.vsd.batch.corrected.for.cibersort[1,] == "S4")
  inhouse_cols <- which(non.log.vsd.batch.corrected.for.cibersort[1,] == "Unknown") 
  # s1_cols_length <- length(s1_cols)
  # s3_cols_length <- length(s3_cols)
  # s4_cols_length <- length(s4_cols)
  # inhouse_cols_length <- length(inhouse_cols)
  # Put signature data in the right order ie all the S1s together etc
  non.log.vsd.batch.corrected.for.cibersort.order <- cbind.data.frame(non.log.vsd.batch.corrected.for.cibersort$Gene,
                                                                      non.log.vsd.batch.corrected.for.cibersort[,s1_cols],
                                                                      non.log.vsd.batch.corrected.for.cibersort[,s3_cols],
                                                                      non.log.vsd.batch.corrected.for.cibersort[,s4_cols]
  )
  non.log.vsd.batch.corrected.for.cibersort.order.inhouse <- cbind.data.frame(non.log.vsd.batch.corrected.for.cibersort$Gene,
                                                                              non.log.vsd.batch.corrected.for.cibersort[,inhouse_cols])
  non.log.vsd.batch.corrected.for.cibersort.order.inhouse[1,2:ncol(non.log.vsd.batch.corrected.for.cibersort.order.inhouse)] <- coldata$names[which(coldata$Study == "InHouse")]
  colnames(non.log.vsd.batch.corrected.for.cibersort.order) <- NULL
  rownames(non.log.vsd.batch.corrected.for.cibersort.order) <- NULL
  colnames(non.log.vsd.batch.corrected.for.cibersort.order.inhouse) <- NULL
  rownames(non.log.vsd.batch.corrected.for.cibersort.order.inhouse) <- NULL
  print("outfile...")
  print(outfile)
  return(list(outfile, non.log.vsd.batch.corrected.for.cibersort.order, outfile_inhouse, non.log.vsd.batch.corrected.for.cibersort.order.inhouse))
}
# define function to convert inhouse dds to CIBERSORT format
inhouse_vsd_to_cibersort_format <- function(dds){
  library(DESeq2)
  library(tibble)
  dds <- filter_out_low_expressed(dds)
  vsd <- vst(dds, blind = TRUE)
  vsd.inhouse.mat <- as.data.frame(assay(vsd))
  vsd.inhouse.mat.non.log <- 2^(vsd.inhouse.mat)
  non.log.vsd <- tibble::rownames_to_column(as.data.frame(vsd.inhouse.mat.non.log), var = "Gene")
  non.log.vsd.inhouse.for.cibersort <- rbind(colnames(non.log.vsd), non.log.vsd)
  colnames(non.log.vsd.inhouse.for.cibersort) <- NULL
  return(non.log.vsd.inhouse.for.cibersort)
}
