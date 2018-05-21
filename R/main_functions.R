#' Evaluate on probe sequence.
#'
#' A function to predict specificity of 20bp sequences.
#' @keywords BindSpace
#' @param input DNAStringSet object of query sequences.
#' @param model list object of embedding model. Contains four elements: kmer_embed, label_embed, kmer2kmerwc and TF 0.05 empirical thresholds
#' @return a matrix of predicted specificity for all sequences and all labels.

#' @import Biostrings
#' @import Matrix
#' @export
eval_probe <- function(input, model) {
  if (any(nchar(input)!=nchar(input[1]))) stop("All input sequences need to be of the same length!") # make sure all inputs have same length
  if (length(grep("N", input)) != 0) stop("input must not contain N!") # make sure no N in inputs
  
  kmer_m <- model$kmer_embed
  label_m <- model$label_embed
  pairwise.kmers <- model$kmer2kmerwc
  
  res <- matrix(NA, length(input), nrow(label_m), dimnames=list(names(input), rownames(label_m)))
  
  N <- length(input) # input sequences size
  batch_size <- 5000
  peak_size <- nchar(input[1])
  batches <- ceiling(N/batch_size)
  p <- 0.5 # normalization factor, same as training
  print(paste0("total ",N, " sequences."))
  time.start <- get.time()
  
  for (i in 1:batches) {
    # generate batch index
    start <- (i-1)*batch_size + 1
    end <- min(i*batch_size, N)
    # generate feature matrix and labels
    tmp <- Matrix(makeFeats(input[start:end], 8), sparse=T)
    feats_tmp <- tmp %*% pairwise.kmers
    feats_tmp[feats_tmp!=0] <- 1 # only consider unique kmers
    tmp_count <- rowSums(feats_tmp) # number of unique kmers
    tmp_embed <- as.matrix(feats_tmp %*% kmer_m / (tmp_count^p)) # embedding for probe, normalized by # kmers and factor p, this is same as starspace
    tmp_res <- tmp_embed %*% t(label_m)# dot product similarity
    res[start:end, ] <- tmp_res
    
    time.end <- get.time()
    print(paste0("Process ",min(i*batch_size, N), " sequences takes ", round( (time.end-time.start)/60, 3), " min."))
  }
  return(res)
}

#' Evaluate on ChIP-seq / ATAC-seq peaks.
#'
#' A function to predict binding affinity to peaks. A score will be computed for every 20bp window of the peak.
#' @keywords BindSpace
#' @param input DNAStringSet object of query sequences (all sequences must be of the same length and don't have N).
#' @param model list object of embedding model. Contains four elements: kmer_embed, label_embed, kmer2kmerwc and TF 0.05 empirical thresholds
#' @param TF a character value. If a particular TF is specsified, a full score matrix will be returned for this TF. (A matrix of dimension #peaks * #windows).
#' @return a list of two elements. 
#' First element is a matrix of size #peaks * #labels. For each peak, a single score is assigned for each label, which is the max score across all 20bp windows.
#' Second element is a binary thresholded matrix of first.
#' Third element is not Null only when TF is specified. It is a matrix of size #peaks * #windows, that is the score matrix for specified TF.

#' @import data.table
#' @import Biostrings
#' @import Matrix
#' @import parallel
#' @export
eval_peak <- function(input, model, TF=NULL) {
  if (any(nchar(input)!=nchar(input[1]))) stop("All input sequences need to be of the same length!") # make sure all inputs have same length
  if (length(grep("N", input)) != 0) stop("input must not contain N!") # make sure no N in inputs
  
  kmer_m <- model$kmer_embed
  label_m <- model$label_embed
  thres <- model$TF_threshold
  pairwise.kmers <- model$kmer2kmerwc
  
  N <- length(input) # input sequences size
  batch_size <- 5000
  bin_size <- 20
  peak_size <- nchar(input[1])
  batches <- ceiling(N/batch_size)
  n_windows <- peak_size - bin_size + 1
  p <- 0.5 # normalization factor, same as training
  print(paste0("total ",N, " sequences."))
  time.start <- get.time()
  res_combine <- list()
  tf_score_combine <- list()
  for (i in 1:batches) {
    # generate batch index
    start <- (i-1)*batch_size + 1
    end <- min(i*batch_size, N)
    tmp_seqs <- input[start:end]
    
    ###############################
    # bin the sequence to windows #
    ###############################
    # for all peaks, for each windows, generate:
    # 1. a matrix showing sum of embedding per peak
    # 2. a vector showing kmer count per peak
    features <- lapply(1 : n_windows, function(x) {
      return(sparseMatrix(length(tmp_seqs), ncol(pairwise.kmers), x = 0))
    })
    for (j in 1:(peak_size-8+1)){
      # take the 8mer at that position
      sub <- as.character(substring(tmp_seqs, j, j + 8 - 1))
      tmp <- pairwise.kmers[sub, , drop=FALSE]; dimnames(tmp) <- list(NULL, NULL)
      features_with_tmp <- max(j-bin_size+8,1):min(j,peak_size-bin_size+1) # decide which 20mer windows include this 8mer
      # add the 8mer to those 20mer windows
      for (k in features_with_tmp) {
        features[[k]] <- features[[k]] + tmp
      }
      if (j%%20==0) print(paste0("scanning 8mer #",j, " out of ", peak_size-8+1," total 8mers..."))
    }
    
    print("computing embedding of sequences...")
    res <- mclapply(features, function(x) {
      x[x!=0] <- 1 # only consider unique kmers in each 20bp window (this is consistent with training)
      tmp_count <- rowSums(x) # # unique kmers per windows
      tmp_embed <- as.matrix(x %*% kmer_m / (tmp_count^p)) # embedding for each window, normalized by # kmers and factor p, this is same as starspace
      return(tmp_embed %*% t(label_m)) # dot product similarity
    }, mc.cores=8)
    
    reduce_res <- matrix(0, length(tmp_seqs), nrow(label_m), dimnames=list(NULL, rownames(label_m)))
    for (j in 1:ncol(reduce_res)) {
      tmp <- do.call(cbind, lapply(res, function(x) return(x[,j])))
      reduce_res[,j] <- apply(tmp, 1, max)
    }
    res_combine[[i]] <- data.table(reduce_res)
    
    if (!is.null(TF)) {
      reduce_TF <- matrix(0, length(tmp_seqs), length(res))
      for (j in 1:ncol(reduce_TF)) {
        reduce_TF[,j] <- res[[j]][,paste0("#TF_",TF)]
      }
      tf_score_combine[[i]] <- data.table(reduce_TF)
    }
    
    time.end <- get.time()
    print(paste0("Process ",min(i*batch_size, N), " sequences takes ", round( (time.end-time.start)/60, 3), " min."))
  }
  
  # 1. output the predicted score matrix
  full_table <- data.matrix(data.frame(rbindlist(res_combine)))
  full_table <- full_table[, grep("TF_",colnames(full_table))]
  colnames(full_table) <- gsub("X\\.TF_","",colnames(full_table))
  
  # 2. output the binary prediction matrix
  names(thres) <- gsub("-",".",names(thres))
  thres <- thres[colnames(full_table)]
  thres_m <- matrix(rep(thres, times=nrow(full_table)), nrow(full_table), ncol(full_table), byrow = T)
  binary_table <- 1 * (full_table>thres_m)
  
  if (!is.null(TF)) {
    TF_m <- data.matrix(data.frame(rbindlist(tf_score_combine)))
  } else {
    TF_m <- NULL
  }
  return(list(score_table=full_table, binary_table=binary_table, TF_score=TF_m))
}

