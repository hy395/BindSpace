install.packages("rpackage/BindSpace/", repos=NULL, type="source")
library(BindSpace)

model.path <- system.file("extdata","model.rds" ,package="BindSpace")
probe.path <- system.file("extdata", "test_Alx1_SELEX.fasta", package="BindSpace")
chip.path <- system.file("extdata", "test_MAFK_ChIP.fasta", package="BindSpace")
atac.path <- system.file("extdata", "test_K562_ATAC.fasta", package="BindSpace")

model <- readRDS(model.path)
probes <- readDNAStringSet(probe.path)
chip <- readDNAStringSet(chip.path)
atac <- readDNAStringSet(atac.path)

#############################
# prediction on 20bp probes #
#############################
res1 <- eval_probe(probes, model) # prediction scores for each probe for each TF

# simple performance check
sort(table(colnames(res1)[apply(res1, 1, which.max)])) # treated as a multiclass task, 1058 sequences are predicted to be Alx1

##########################
# prediction on ChIP-seq #
##########################
res2 <- eval_peak(chip, model, "MAFK")

# 1. score table shows predicted for each PEAK for each TF
a <- res2$score_table
sort(table(colnames(a)[apply(a,1,which.max)])) # 1621/2000 peaks are predicted to be MAFK peak

# 2. TF_score shows predicted score for each 20bp WINDOW for each PEAK for the TF "MAFK"
b <- res2$TF_score
plotHeatmap(b)  # plot binding signals

##########################
# prediction on ATAC-seq #
##########################
res3 <- eval_peak(atac, model)
a <- res3$score_table # predicted score for each peak for each TF
b <- res3$binary_table # binary preditions


