// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <R.h>
#include <Rcpp.h>
#define ARMA_64BIT_WORD
using namespace Rcpp;

// This turns a dgCMatrix DTM into a list of indices for looping over in
// an LDA gibbs sampler
// [[Rcpp::export]]
List dtm_to_lexicon_c(arma::sp_mat x) {
  
  // initialize some variables
  
  int d,v; // iteration indices
  
  int Nd = x.n_rows; // number of docs
  
  int Nw = x.n_cols; // number of words
  
  IntegerVector doc_sums(Nd); // gives length of each output vector
  
  List out(Nd); // list to hold output vectors
  
  // loop over each document
  for (d = 0; d < Nd; d++) {
    
    // count tokens in each document for length of output vector
    for (v = 0; v < Nw; v++) {
      doc_sums[d] += x(d,v);
    }
    
    IntegerVector tmp(doc_sums[d]); // temporary vector to be inserted into output
    
    int k = 0; // index of tmp, advances when we have non-zero entries 
    
    for (v = 0; v < Nw; v++) {
      
      if (x(d,v) > 0) { // if non-zero, add elements to tmp
        
        int idx = k + x(d,v); // where to stop the loop below
        
        while (k < idx) {
          
          tmp[k] = v;
          
          k += 1;
        }
        
      }
      
    }
    
    out[d] = tmp;
  }
  
  return out;
  
}



// This is a collapsed gibbs sampler for LDA.
// Pre-processing and post-processing is assumed in R
// [[Rcpp::export]]
List fit_lda_c(List &docs, int &Nk, int &Nd, int &Nv, NumericVector &alpha, 
               NumericVector &beta, int &iterations, int &burnin) {
  
  // print a status so we can see where we are
  // std::cout << "declaring variables\n";
  
  // Declare some initial variables
  double sum_alpha = sum(alpha); // rcpp sugar here, I guess

  NumericVector k_beta = beta * Nk; // rcpp sugar here, I guess

  // Declare data structures
  int i, d, n, k; // indices for loops
  
  NumericVector p_z(Nk);
  
  IntegerVector z1; // placeholder for topic sampled

  int z = 0; // placeholder for topic sampled

  int v = 0; // placeholder for word index

  IntegerVector topic_sample = seq_len(Nk) -1;

  List z_dn(Nd) ; // count of topic/term assignments by document

  IntegerMatrix theta_counts(Nd, Nk); // count of topics over documents

  IntegerMatrix phi_counts(Nk, Nv); // count of terms over topics

  IntegerVector n_d(Nd); // total count of term document totals

  IntegerVector n_z(Nk); // total count of topic totals

  IntegerMatrix theta_sums(Nd, Nk); // initialize matrix for averaging over iterations

  IntegerMatrix phi_sums(Nk, Nv); // initialize matrix for averaging over iterations
  
  // Assign initial values at random
  // std::cout << "assigning initial values \n";
  
  for(d = 0; d < Nd; d++){
    IntegerVector doc = docs[d];
    
    IntegerVector z_dn_row(doc.length());
    
    for(n = 0; n < doc.length(); n++){
      
      // sample a topic at random
      z1 = RcppArmadillo::sample(topic_sample, 1, false, NumericVector::create());
      
      z = z1[0];
      
      theta_counts(d,z) += 1;
      
      v = doc[n];
      
      phi_counts(z,v) += 1;
      
      n_d[d] = n_d[d] + 1; // count the number of tokens in the document
      
      n_z[z] += 1; // count the that topic overall
      
      z_dn_row[n] = z; // # count that topic for that word in the document
      
    }
    
    z_dn[d] = z_dn_row; // update topic-doc-word tracking
    
  }
  
  // Gibbs iterations
  // std::cout << "beginning Gibbs \n";
  for (i = 0; i < iterations; i++) { // for each iteration
    
    for (d = 0; d < Nd; d++) { // for each document
      // std::cout << "document " << d << "\n";
      
      R_CheckUserInterrupt();
      
      IntegerVector doc = docs[d]; // placeholder for a document
      
      IntegerVector z_dn_row = z_dn[d]; // placeholder for doc-word-topic assigment
      
      for (n = 0; n < n_d[d]; n++) { // for each word in that document
        
        // discount for the n-th word with topic z
        z = z_dn_row[n];
        
        theta_counts(d,z) -= 1; 
        
        phi_counts(z,doc[n]) -= 1;
        
        n_z[z] -= 1;
        
        // n_d[d] -= 1;
        
        // sample topic index
        for (k = 0; k < Nk; k++) {
          
          p_z[k] = (phi_counts(k,doc[n]) + beta[doc[n]]) / (n_z[k] + k_beta[doc[n]]) *
            (theta_counts(d,k) + alpha[k]) / (n_d[d] + sum_alpha);
          
          // p_z[k] = (phi_counts(k,doc[n]) + beta[doc[n]]) * (theta_counts(d,k) + alpha[k]);
          
        }
        
        
        // update counts
        z1 = RcppArmadillo::sample(topic_sample, 1, false, p_z);
        
        z = z1[0];
        
        theta_counts(d,z) += 1; // update document topic count
        
        phi_counts(z,doc[n]) += 1; // update topic word count
        
        // n_d[d] = n_d[d] + 1; // count that topic in that document overall
        
        n_z[z] += 1; // count the that topic overall
        
        z_dn_row[n] = z; // # count that topic for that word in the document
        
      }
    }
    
    // if using burnin, update sums
    if (burnin > -1 && i >= burnin) {
      
      for (k = 0; k < Nk; k++) {
        
        for (d = 0; d < Nd; d++) {
          theta_sums(d,k) += theta_counts(d,k);
        }
        
        for (v = 0; v < Nv; v++) {
          phi_sums(k,v) += phi_counts(k,v);
        }
      }
      
    }
    
  }
  
  // return the result
  // std::cout << "prepare result\n";
  
  if (burnin > -1) {
    int i_diff = iterations - burnin;
    
    NumericMatrix theta(Nd,Nk);
    
    NumericMatrix phi(Nk,Nv);
    
    // average over chain after burnin 
    for (k = 0; k < Nk; k++) {
      
      for (d = 0; d < Nd; d++) {
        theta(d,k) = (theta_sums(d,k) / i_diff);
      }
      
      for (v = 0; v < Nv; v++) {
        phi(k,v) = (phi_sums(k,v) / i_diff);
      }
    }
    
    return List::create(Named("theta") = theta,
                        Named("phi") = phi);
    
  } else {
    return List::create(Named("theta") = theta_counts,
                        Named("phi") = phi_counts);
  }
  
  
}

/*** R

source("~/Documents/gmu/study_summer_2018/R/FitLdaModel.R")

dtm <- textmineR::nih_sample_dtm

k <- 15

alpha <- rgamma(k, 5, 20)

beta <- Matrix::colSums(dtm) / sum(Matrix::colSums(dtm)) * 256

m <- FitLdaModel(dtm = dtm, 
                 k = k, 
                 iterations = 1000,
                 burnin = 900,
                 alpha = alpha, 
                 beta = beta,
                 calc_coherence = TRUE,
                 calc_r2 = TRUE,
                 seed = 1234)
  

microbenchmark::microbenchmark(FitLdaModel(dtm = dtm, 
                                           k = k, 
                                           iterations = 100,
                                           burnin = -1,
                                           alpha = alpha, 
                                           beta = beta,
                                           calc_coherence = FALSE,
                                           calc_r2 = FALSE,
                                           seed = 1234))
*/

