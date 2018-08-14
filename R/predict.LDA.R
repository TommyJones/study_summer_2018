

predict.LDA <- function(object, newdata, method = c("gibbs", "dot"), 
                        iterations = NULL, seed = NULL, ...) {
  
  ### Check inputs ----
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
    stop("newdata not of class dgCMatrix is not yet supported")
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
    
    Nv <- ncol(object$phi)
    
    sum_alpha <- sum(object$alpha)
    
    # declare data structures
    z_dn <- lapply(docs, function(x) numeric(length(x))) # count of topic/term assignments by document, z_m_n
    
    theta_counts <- matrix(0, nrow = Nd, ncol = Nk) # count of topics over documents, n_m_z
    
    n_d <- numeric(Nd) # count of term totals
    
    # random inital values
    for (d in seq_along(docs)) {
      for (n in seq_along(docs[[d]])) {
        # sample a topic
        z <- sample(seq_len(Nk),1)
        
        # update counts, keeping phi objects static
        theta_counts[d,z] <- theta_counts[d,z] + 1 # count that topic in the document
        
        n_d[d] <- n_d[d] + 1 # count that topic in that document overall
        
        z_dn[[d]][n] <- z # count that topic for that word in the document
        
      }
    }
    
    # gibbs iterations
    for (i in seq_len(iterations)) {
      for (d in seq_along(docs)) {
        for(n in seq_along(docs[[d]])){
          # discount for the n-th word with topic z
          z <- z_dn[[d]][n]
          
          theta_counts[d,z] <- theta_counts[d,z] - 1
          
          n_d[d] <- n_d[d] - 1
          
          # sample topic index
          d_a <- n_d[d] + sum_alpha # denominator
          
          p_z <- object$phi[ ,docs[[d]][n]] * (theta_counts[d,] + object$alpha) / (d_a)
          
          z <- sample(seq_len(Nk), 1, prob = p_z)
          
          # update counts
          theta_counts[d,z] <- theta_counts[d,z] + 1 # count that topic in the document
          
          n_d[d] <- n_d[d] + 1 # count that topic in that document overall
          
          z_dn[[d]][n] <- z
          
        }
      }
    }
    
    # format result
    theta <- t(t(theta_counts) + object$alpha)
    
    theta <- theta_counts / rowSums(theta_counts, na.rm = TRUE)
    
    theta[ is.na(theta) ] <- 0
    
    result <- theta
  }

  # return result
  result
  
}