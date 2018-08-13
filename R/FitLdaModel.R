################################################################################
# This file turns the Gibbs sampler cobbled from various sources into an R function
# Once I am sure that it is (a) calculating properly and (b) efficient, I will
# roll it over to C++
################################################################################



Dtm2Lexicon <- function(dtm, ...) {
  # # get doc lengths
  # doc_lengths <- Matrix::rowSums(dtm)
  # 
  # # do in parallel in batches of about 3000 if we have more than 3000 docs
  # if(nrow(dtm) > 3000){
  #   
  #   batches <- seq(1, nrow(dtm), by = 3000)
  #   
  #   dtm_list <- lapply(batches, function(x) dtm[ x:min(x + 2999, nrow(dtm)) , ])
  #   
  #   out <-TmParallelApply(X = dtm_list, FUN = function(x){
  #     Dtm2LexiconC(dtm = x, doc_lengths = doc_lengths)
  #   }, ...)
  #   
  # }else{
  #   out <- Dtm2LexiconC(dtm = dtm, doc_lengths = doc_lengths)
  # }
  # 
  # out <- unlist(out)
  # 
  # names(out) <- rownames(dtm)
  # 
  # out
  
  # above is real code, but I have C++ stuff to sort out
  docs <- apply(X, 1, function(x){
    out <- c()
    for (j in seq_along(x)) {
      out <- c(out, rep(j, x[j]))
    }
    out
  })
  
  docs
}

#' @param dtm A document term matrix of class dgCMatrix
#' @param k Integer number of topics
#' @param alpha Vector of length \code{k} for asymmetric or a number for symmetric.
#'        This is the prior for topics over documents
#' @param beta Vector of length \code{ncol(dtm)} for asymmetric or a number for symmetric.
#'        This is the prior for words over topics.
#' @param iterations Integer number of iterations for the Gibbs sampler to run. A
#'        future version may include automatic stopping criteria.
#' @param seed If not null (the default) then the random seed you wish to set
#' @param ... Other arguments to be passed to textmineR::TmParallelApply
FitLdaModel <- function(dtm, k, iterations = NULL, alpha = 0.1, beta = 0.05, 
                        seed = NULL, ...){
  
  ### Check inputs are of correct dimensionality ----
  
  # dtm of the correct format?
  if (! "dgCMatrix" %in% class(dtm)) {
    message("dtm is not of class dgCMatrix, attempting to convert...")
    
    dtm <- try(methods::as(dtm, "dgCMatrix", strict = TRUE))
    
    if (! "dgCMatrix" %in% class(dtm))
      stop("conversion failed. Please pass an object of class dgCMatrix for dtm")
  }
  
  # is k formatted correctly?
  if (! is.numeric(k))
    stop("k must be an integer")
  
  k <- floor(k) # in case somebody is cheeky and passes a decimal
  
  # iterations?
  if (is.null(iterations))
    stop("You must specify number of iterations")
  
  # alpha and beta?
  if (! is.numeric(alpha) | sum(is.na(alpha)) > 0)
    stop("alpha must be a numeric scalar or vector with no missing values")
  
  if (length(alpha) == 1) {
    alpha <- rep(alpha, k)
  } else if (length(alpha) != k){
    stop("alpha must be a scalar or vector of length k")
  }
  
  if (! is.numeric(beta) | sum(is.na(beta)) > 0)
    stop("beta must be a numeric scalar or vector with no missing values")
  
  if (length(beta) == 1) {
    beta <- rep(beta, ncol(dtm))
  } else if (length(beta) != ncol(dtm)){
    stop("beta must be a scalar or vector of length ncol(dtm)")
  }
  
  if (! is.null(seed)){
    if (! is.numeric(seed)){
      stop("if seed is not NULL, then it must be numeric")
    } else {
      set.seed(seed)
    }
  }
  
  
  ### Format inputs ----
  
  docs <- Dtm2Lexicon(dtm)
  
  Nd <- nrow(dtm)
  
  Nk <- 2
  
  Nv <- ncol(dtm)
  
  ### Declare data structures ----
  z_dn <- lapply(docs, function(x) x * 0) # count of topic/term assignments by document, ta
  
  theta_counts <- matrix(0, nrow = Nd, ncol = Nk) # count of topics over documents, n_m_z
  
  phi_counts <- matrix(0, nrow = Nk, ncol = Nv) # count of terms over topics, wt
  
  # n_d <- numeric(Nd) # count of term totals, I don;t know that I need this...
  
  n_z <- numeric(Nk) # count of topic totals
  
  ### Assign initial values ----
  
  for (d in seq_along(docs)) {
    for (v in seq_along(docs[[d]])) {
      
      z <- sample(seq_len(Nk), 1) # sample a topic
      
      z_dn[[d]][v] <- z # word-topic assignment for the document
      
      phi_counts[z, docs[[d]][v]] <- phi_counts[z, docs[[d]][v]] + 1 # count that word and topic
      
      # This is different from Andrew's, check to see if they calculate the same
      theta_counts[d,z] <- theta_counts[d,z] + 1 # count that topic in the document
      
    }
  }
  
  
  ### Burn in iterations ----
  # this is going to take some figuring out...
  
  ### Gibbs iterations ----
  for (i in seq_len(iterations)) {
    for (d in seq_along(docs)) {
      for (v in seq_along(docs)[[d]]) {
        
        # Get current values
        w <- docs[[d]][v] # index of the current word
        
        z <- z_dn[[d]][v] # topic assignment of the current word
        
        # discount current values from sampling
        theta_counts[d,z] <- theta_counts[d,z] - 1
        
        phi_counts[z,w] <- phi_counts[z,w] - 1
        
        # sample a new topic
        d_a <- sum(theta_counts[d,] + alpha)
        
        d_b <- rowSums(phi_counts + beta[w])
        
        p_z <- (phi_counts[,w] + beta[w]) / (d_b) * (theta_counts[d,] + alpha) / (d_a)
        
        z <- sample(seq_len(Nk), 1, prob = p_z)
        
        # update topic assignments
        theta_counts[d,z] <- theta_counts[d,z] + 1 # count that topic in the document
        
        z_dn[[d]][n] <- z
        
        phi_counts[z,docs[[d]][n]] <- phi_counts[z,docs[[d]][n]] + 1 # count that word and topic
        
      }
    }
  }
  
  ### Format posteriors correctly ----
  # calculate output parameters
  phi <- t(t(phi_counts) + beta)
  
  phi <- phi_counts / rowSums(phi_counts, na.rm = TRUE)
  
  phi[ is.na(phi) ] <- 0
  
  theta <- t(t(theta_counts) + alpha)
  
  theta <- theta_counts / rowSums(theta_counts, na.rm = TRUE)
  
  theta[ is.na(theta) ] <- 0
  
  ### return the result ----
  
  result <- list(phi = phi, theta = theta, dtm = dtm) # add other things here
  
  class(result) <- c("LDA", "TopicModel")
  
  result
}