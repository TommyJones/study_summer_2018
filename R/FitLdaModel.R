################################################################################
# This file turns the Gibbs sampler I found in Agustinus into an R function
# Once I am sure that it is (a) calculating properly and (b) efficient, I will
# roll it over to C++
################################################################################

#' @param dtm A document term matrix of class dgCMatrix
#' @param k Integer number of topics
#' @param alpha Vector of length \code{k} for asymmetric or a number for symmetric.
#'        This is the prior for topics over documents
#' @param beta Vector of length \code{ncol(dtm)} for asymmetric or a number for symmetric.
#'        This is the prior for words over topics.
#' @param iterations Integer number of iterations for the Gibbs sampler to run. A
#'        future version may include automatic stopping criteria.
#' @param seed If not null (the default) then the random seed you wish to set
FitLdaModel <- function(dtm, k, alpha, beta, iterations, seed) {
  
  # Make sure inputs are formatted correctly ----
  
  # 
  
}
