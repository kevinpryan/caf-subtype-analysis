metadata <- read.table("~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/metadata/metadata_full.txt", row.names = 1)
hiseq <- which(metadata$Study == "EGAD00001004810" | metadata$Study == "EGAD00001003808" | metadata$Study == "EGAD00001005744")
novaseq <- which(metadata$Study == "EGAD00001006144")
miseq <- which(metadata$Study == "InHouse")
metadata$Sequencer <- rep("a", nrow(metadata))
metadata$Sequencer[hiseq] <- "hiseq"
metadata$Sequencer[miseq] <- "miseq"
metadata$Sequencer[novaseq] <- "novaseq"
write.table(metadata, file = "~/Documents/PhD/subtypes/caf-subtype-analysis/nf-subpop/outdir/metadata/metadata_with_sequencer.txt", quote = F, sep = "\t")
