#### Declare some data and input parameters ----
# X <- rbind(c(0, 0, 1, 2, 2),
#            c(0, 0, 1, 1, 1),
#            c(0, 1, 2, 2, 2),
#            c(4, 4, 4, 4, 4),
#            c(3, 3, 4, 4, 4),
#            c(3, 4, 4, 4, 4))

# X <- textmineR::nih_sample_dtm

X <- rbind(c(0, 0, 1, 2, 2),
           c(0, 0, 1, 1, 1),
           c(0, 1, 2, 2, 2),
           c(4, 4, 4, 0, 0),
           c(3, 3, 4, 0, 0),
           c(3, 4, 4, 1, 0))


docs <- apply(X, 1, function(x){
  out <- c()
  for (j in seq_along(x)) {
    out <- c(out, rep(j, x[j]))
  }
  out
})

Nd <- nrow(X)

Nk <- 2

Nv <- ncol(X)

iterations <- 1000

alpha <- numeric(Nk) + 1

beta <- numeric(Nv) + 0.001

################################################################################
# Prototype collapsed gibbs sampler
################################################################################

# initial inputs
# docs - a list like what is returned from LDA's lexicalize function
#   each element is a vector of word *indices* e.g. "3" refers to the 3rd column in the DTM
# k - number of topics (converted to Nk)
# iterations - maximum number of gibbs iterations
# burnin - number of burnin interations
# seed - random seed

# Declare data structures
z_dn <- lapply(docs, function(x) x * 0) # count of topic/term assignments by document, z_m_n

theta_counts <- matrix(0, nrow = Nd, ncol = Nk) # count of topics over documents, n_m_z

phi_counts <- matrix(0, nrow = Nk, ncol = Nv) # count of terms over topics, n_z_t

n_d <- numeric(Nd) # count of term totals, I don;t know that I need this...

n_z <- numeric(Nk) # count of topic totals

# Initialize values
# Note: why initialize at completely random values? Why not use priors?
for (d in seq_along(docs)) { # for every document
  for (n in seq_along(docs[[d]])) { # for every word in that document
    
    z <- sample(seq_len(Nk), 1) # sample a topic at random
    
    theta_counts[d,z] <- theta_counts[d,z] + 1 # count that topic in the document
    
    n_d[d] <- n_d[d] + 1 # count that topic in that document overall
    
    z_dn[[d]][n] <- z # count that topic for that word in the document
    
    phi_counts[z,docs[[d]][n]] <- phi_counts[z,docs[[d]][n]] + 1 # count that word and topic
    
    n_z[z] <- n_z[z] + 1 # count the that topic overall
  }
}

# Burn-in iterations

# Gibbs iterations

# converge <- FALSE

it <- 1

while (it < iterations) { #  | ! converge
  
  # for each document
  for (d in seq_along(docs)) {
    # for each word in that document
    for (n in seq_along(docs[[d]])) {
      # discount for the n-th word with topic z
      z <- z_dn[[d]][n]
      
      theta_counts[d,z] <- theta_counts[d,z] - 1
      
      phi_counts[z,docs[[d]][n]] <- phi_counts[z,docs[[d]][n]] - 1
      
      n_z[z] <- n_z[z] - 1
      
      n_d[d] <- n_d[d] - 1
      
      # sample topic index
      d_a <- sum(theta_counts[d,] + alpha) # n_d[d] + Nk * alpha # 
      
      d_b <- rowSums(phi_counts + beta[docs[[d]][n]]) # n_z[z] + Nv *  beta[docs[[d]][n]]# 
      
      p_z <- (phi_counts[,docs[[d]][n]] + beta[docs[[d]][n]]) / (d_b) * (theta_counts[d,] + alpha) / (d_a)
      
      z <- sample(seq_len(Nk), 1, prob = p_z)
      
      # update counts
      theta_counts[d,z] <- theta_counts[d,z] + 1 # count that topic in the document
      
      # n_d[d] <- n_d[d] + 1 # count that topic in that document overall
      
      z_dn[[d]][n] <- z
      
      phi_counts[z,docs[[d]][n]] <- phi_counts[z,docs[[d]][n]] + 1 # count that word and topic
      
      n_z[z] <- n_z[z] + 1 # count the word in that topic overall
      
    }
  }
  
  # check convergence criteria and increase iteration count
  # if (a < b)
  #   converge <- TRUE
  
  it <- it + 1
}

# calculate output parameters
phi <- t(t(phi_counts) + beta)

phi <- phi_counts / rowSums(phi_counts, na.rm = TRUE)

phi[ is.na(phi) ] <- 0

theta <- t(t(theta_counts) + alpha)

theta <- theta_counts / rowSums(theta_counts, na.rm = TRUE)

theta[ is.na(theta) ] <- 0

# if algorithm did not converge, issue a warning
# if (! converge)
#   warning("Algorithm did not converge. Results may not be reliable.")
