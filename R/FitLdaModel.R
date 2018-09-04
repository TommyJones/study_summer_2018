################################################################################
# This file turns the Gibbs sampler cobbled from various sources into an R function
# Once I am sure that it is (a) calculating properly and (b) efficient, I will
# roll it over to C++
################################################################################


#' Turn a document term matrix into a list for LDA Gibbs sampling
#' @description Represents a document term matrix as a list. 
#' @param dtm A document term matrix (or term co-occurrence matrix) of class 
#' \code{dgCMatrix}. 
#' @param ... Other arguments to be passed to \code{\link[textmineR]{TmParallelApply}}.
#' @return Returns a list. Each element of the list represents a row of the input
#' matrix. Each list element contains a numeric vector with as many entries as
#' tokens in the original document. The entries are the column index for that token, minus 1. 
#' @examples
#' \dontrun{
#' # Load pre-formatted data for use
#' data(nih_sample_dtm)
#' 
#' result <- Dtm2Lexicon(dtm = nih_sample_dtm, 
#'                       cpus = 2)
#' }
#' @export 
Dtm2Lexicon <- function(dtm, ...) {

  # do in parallel in batches of about 3000 if we have more than 3000 docs
  if(nrow(dtm) > 3000){

    batches <- seq(1, nrow(dtm), by = 3000)

    dtm_list <- lapply(batches, function(x) dtm[ x:min(x + 2999, nrow(dtm)) , ])

    out <-textmineR::TmParallelApply(X = dtm_list, FUN = function(y){
      dtm_to_lexicon_c(x = y)
    }, ...)
    
    out <- do.call(c, out)

  }else{
    out <- dtm_to_lexicon_c(x = dtm)
  }

  names(out) <- rownames(dtm)

  out
  
}

#' Fit a Latent Dirichlet Allocation topic model
#' @description Fit a Latent Dirichlet Allocation topic model using collapsed Gibbs sampling. 
#' @param dtm A document term matrix or term co-occurrence matrix of class dgCMatrix
#' @param k Integer number of topics
#' @param alpha Vector of length \code{k} for asymmetric or a number for symmetric.
#'        This is the prior for topics over documents
#' @param beta Vector of length \code{ncol(dtm)} for asymmetric or a number for symmetric.
#'        This is the prior for words over topics.
#' @param iterations Integer number of iterations for the Gibbs sampler to run. A
#'        future version may include automatic stopping criteria.
#' @param seed If not null (the default) then the random seed you wish to set. This
#' will return the same outputs for the same inputs. Useful for diagnostics.
#' @param ... Other arguments to be passed to textmineR::TmParallelApply
#' @return Returns an S3 object of class c("LDA", "TopicModel"). DESCRIBE MORE
#' @details EXPLAIN IMPLEMENTATION DETAILS
#' @examples GIVE EXAMPLES
#' @export
FitLdaModel <- function(dtm, k, iterations = NULL, burnin = -1, alpha = 0.1, beta = 0.05, 
                        optimize_alpha = FALSE, calc_likelihood = FALSE, 
                        calc_coherence = TRUE, calc_r2 = FALSE, seed = NULL, ...){
  
  ### Check inputs are of correct dimensionality ----
  
  # iterations and burnin acceptable?
  if (burnin >= iterations) {
    stop("burnin must be less than iterations")
  }
  
  # dtm of the correct format?
  if (! "dgCMatrix" %in% class(dtm)) {
    message("dtm is not of class dgCMatrix, attempting to convert...")
    
    dtm <- try(methods::as(dtm, "dgCMatrix", strict = TRUE)) # requires Matrix in namespace
    
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
    alpha <- numeric(k) + alpha
  } else if (length(alpha) != k){
    stop("alpha must be a scalar or vector of length k")
  }
  
  if (! is.numeric(beta) | sum(is.na(beta)) > 0)
    stop("beta must be a numeric scalar or vector with no missing values")
  
  if (length(beta) == 1) {
    beta <- numeric(ncol(dtm)) + beta
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
  
  if (! is.logical(calc_coherence))
    stop("calc_coherence must be logical")
  
  if (! is.logical(calc_r2))
    stop("calc_r2 must be logical")
  
  
  ### Format inputs ----
  
  docs <- Dtm2Lexicon(dtm)
  
  Nd <- nrow(dtm)
  
  Nk <- k
  
  Nv <- ncol(dtm)

  ### run C++ gibbs sampler ----
  
  result <- fit_lda_c(docs = docs, Nk = Nk, Nd = Nd, Nv = Nv, 
                      alph = alpha, beta = beta,
                      iterations = iterations, burnin = burnin,
                      optimize_alpha = optimize_alpha,
                      calc_likelihood = calc_likelihood)
  
  
  ### Format posteriors correctly ----
  
  phi <- result$phi
  
  theta <- result$theta
  
  phi <- t(t(phi) + beta)
  
  phi <- phi / rowSums(phi, na.rm = TRUE)
  
  phi[ is.na(phi) ] <- 0
  
  theta <- t(t(theta) + alpha)
  
  theta <- theta / rowSums(theta, na.rm = TRUE)
  
  theta[ is.na(theta) ] <- 0
  
  colnames(phi) <- colnames(dtm)
  
  rownames(phi) <- paste0("t_", seq_len(Nk))
  
  colnames(theta) <- rownames(phi)
  
  rownames(theta) <- rownames(dtm)
  
  ### collect the result ----
  gamma <- textmineR::CalcPhiPrime(phi = phi, theta = theta, 
                                   p_docs = Matrix::rowSums(dtm))
  
  result <- list(phi = phi, theta = theta, gamma = gamma,
                 dtm = dtm, alpha = result$alpha, beta = result$beta,
                 log_likelihood = data.frame(result$log_likelihood)) # add other things here
  
  names(result$log_likelihood) <- c("iteration", "log_likelihood")
  
  class(result) <- c("LDA", "TopicModel")
  
  ### calculate additional things ----
  if (calc_coherence) {
    result$coherence <- textmineR::CalcProbCoherence(result$phi, dtm, M = 5)
  }
  
  if (calc_r2) {
    result$r2 <- textmineR::CalcTopicModelR2(dtm, result$phi, result$theta, ...)
  }
  
  if (! calc_likelihood) {
    result$log_likelihood <- NULL
  }
  
  ### return result ----
  result
}

### Predict method for LDA objects
predict.LDA <- function(object, newdata, method = c("gibbs", "dot"), 
                        iterations = NULL, burnin = -1, seed = NULL, ...) {
  
  ### Check inputs ----
  if (method[1] == "gibbs") {
    
    if (is.null(iterations)) {
      stop("when using method 'gibs' iterations must be specified.")
    }
    
    if (burnin >= iterations) {
      stop("burnin must be less than iterations")
    }
    
  }
  
  
  if (sum(c("LDA", "TopicModel") %in% class(object)) < 2) {
    stop("object must be a topic model object of class c('LDA', 'TopicModel')")
  }
  
  if (sum(c("dgCMatrix", "character") %in% class(newdata)) < 1) {
    stop("newdata must be a matrix of class dgCMatrix or a character vector")
  }
  
  if (sum(c("gibbs", "dot") %in% method) == 0) {
    stop("method must be one of 'gibbs' or 'dot'")
  }
  
  if (! is.null(seed)) {
    if (! is.numeric(seed)){
      stop("seed must be NULL or numeric")
    }
    set.seed(seed)
  }
  
  ### If newdata is a character vector, convert to dgCMatrix ----
  # TODO: requires re-write of outputs of CreateDtm/CreateTcm
  
  if ("dgCMatrix" %in% class(newdata)) {
    dtm_newdata <- newdata
  } else {
    # code to convert goes here
    stop("newdata not of class character is not yet supported")
  }
  
  ### Align vocabulary ----
  vocab1 <- colnames(object$dtm)
  # vocab2 <- setdiff(colnames(dtm_newdata), colnames(object$dtm)) # for now unused
  
  ### Get predictions ----
  
  if (method[1] == "dot") { # dot product method
    
    result <- dtm_newdata[ ,vocab1]
    result <- (result / Matrix::rowSums(result)) %*% t(object$gamma[ ,vocab1])
    result <- as.matrix(result)
    
  } else { # gibbs method
    # format inputs
    docs <- Dtm2Lexicon(dtm_newdata[,vocab1])
    
    Nd <- nrow(dtm_newdata)
    
    Nk <- nrow(object$phi)
    
    sum_alpha <- sum(object$alpha)
    
    # pass inputs to C function
    theta <- predict_lda_c(docs = docs, Nk = Nk, Nd = Nd, 
                           alpha = object$alpha, phi = object$phi,
                           iterations = iterations, burnin = burnin)
    
    theta <- theta$theta
    
    # format outputs
    theta <- t(t(theta) + object$alpha)
    
    theta <- theta / rowSums(theta, na.rm = TRUE)
    
    theta[ is.na(theta) ] <- 0
    
    result <- theta
  }
  
  # return result
  
  rownames(result) <- rownames(dtm_newdata)
  colnames(result) <- rownames(object$phi)
  
  result
  
}
