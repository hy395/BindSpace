#' Compute kmer frequency.
#'
#' A function to predict specificity of 20bp sequences.
#' @keywords BindSpace
#' @param seqs DNAStringSet object. query sequences.
#' @param k int value. kmer length.
#' @param rc logical value. whether reverseComplement sequence should also be counted.
#' @return a kmer frequency table.
#' @import Biostrings
#' @export 
makeFeats <- function(seqs, k, rc=F) {
  a <- oligonucleotideFrequency(seqs, k)
  feats <- a
  if (rc) {
    seqs_rc <- reverseComplement(seqs)
    b <- oligonucleotideFrequency(seqs_rc, k)
    feats <- a+b
  }
  return(feats)
}

#' plot binding score Heatmap.
#' 
#' A function to plot binding score.
#' @keywords BindSpace
#' @param m matrix of per window per peak binding score for a given TF
#' @return a ComplexHeatmap object
#' @import ComplexHeatmap
#' @import circlize
#' @export
plotHeatmap <- function (m) {
  ht1 <- Heatmap(m, cluster_columns=F, cluster_rows=F, show_row_names=F, show_heatmap_legend=F, show_column_names=F, 
                 col=colorRamp2(c(0, max(m)/2, max(m)), c("blue","white", "red")))
  return(ht1)
}

#' Compute F1 score.
#' 
#' A function to compute F1 score.
#' @keywords BindSpace
#' @param pred numeric vector of predcted labels (0 and 1).
#' @param obs numeric vector of true labels (0 and 1).
#' @return numeric vector of precision, recall, and F1 score.
#' @export
F1 <- function(pred, obs) {
  if (any(!pred%in%c(0,1)) | any(!obs%in%c(0,1)) ) {
    stop("coding must be 0/1")
  }
  prec <- sum(pred==1&obs==1)/sum(pred==1)
  reca <- sum(pred==1&obs==1)/sum(obs==1)
  F1 <- 2*prec*reca/(prec+reca)
  return(list(prec=prec,reca=reca,F1=F1))
}


#' Compute Area Under Prediction Recall Curve.
#' 
#' A function to compute AUPRC.
#' @keywords BindSpace
#' @param pred numeric vector of predcted scores.
#' @param obs numeric vector of true labels (0 and 1).
#' @return numeric value of Area Under Precision Recall Curve.
#' @import PRROC
#' @export
aupr <- function(pred, obs) {
  # scores.class0 is postive class
  # score.class1 is negative class
  pr <- pr.curve(scores.class0 = pred[obs==1], scores.class1 = pred[obs==0], curve = T)
  return(pr$auc.integral)
}

#' Get current system time.
#' 
#' A function to get current system time. This is the same function as from SeqGL.
#' @keywords BindSpace
#' @return None.
#' @export
get.time <- function () {
  return(as.numeric(format(Sys.time(), "%s")))
}