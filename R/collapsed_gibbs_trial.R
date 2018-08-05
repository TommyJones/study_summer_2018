################################################################################
# Prototype collapsed gibbs sampler
################################################################################

# docs - list of word indices, akin to Chang's documents
# K - number of topics
# V - number of unique words (vocabulary size)
# iterations - number of Gibbs iterations

### Input formatting ----

# dtm to docs object

K <- k

Nd <- nrow(dtm)

V <- ncol(dtm)

### Initialize variables ----

# start with zero values
z_dn <- vector(mode = "list", length = Nd) # topics of words of documents
# note most efficient would be to just copy docs and overwrite, 
# but leaves possiblity of undetected error

theta_c <- matrix(0, ncol = K, nrow = Nd) # matrix of doc-topic counts

phi_c <- matrix(0, ncol = V, nrow = K) # matrix of topic-word counts

n_z <- numeric(K) # word count of each topic

# randomly assign initial values
for (d in seq_along(docs)) { # for each document
  z_n <- numeric(K)
  
}

