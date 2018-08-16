// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
using namespace Rcpp;


// This is a collapsed gibbs sampler for LDA.
// Pre-processing and post-processing is assumed in R

// [[Rcpp::export]]
List lda_c(List docs, int Nk, int Nd, int Nv, NumericVector alpha, 
           NumericVector beta, int iterations, int burnin) {

  // Declare some initial variables
  double sum_alpha = sum(alpha); // rcpp sugar here, I guess
  
  NumericVector k_beta = beta * Nk; // rcpp sugar here, I guess
  
  // Declare data structures
  NumericVector p_z(Nk);
  
  IntegerVector z1 = 0; // placeholder for topic sampled
  
  int z2 = 0; // converted type for topic sampled
  
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
  for(int d = 0; d < docs.length(); d++){
    IntegerVector doc = docs[d];
    
    IntegerVector z_dn_row(doc.length());
    
    for(int n = 0; n < doc.length(); n++){
      
      // sample a topic at random
      z1 = RcppArmadillo::sample(topic_sample, 1, false, NumericVector::create());
      
      z2 = z1[1];
      
      theta_counts(d,z2) = theta_counts(d,z2) + 1;
      
      v = doc[n];
      
      phi_counts(z2,v) = phi_counts(z2,v) + 1;
      
      n_d[d] = n_d[d] + 1; // count that topic in that document overall
      
      n_z[z2] = n_z[z2] + 1; // count the that topic overall
      
      z_dn_row[n] = z2; // # count that topic for that word in the document
      
    }
    
    z_dn[d] = z_dn_row; // update topic-doc-word tracking
    
  }
  
  
  // Gibbs iterations
  for (int i = 0; i < iterations; i++) { // for each iteration
    
    for (int d = 0; d < Nd; d++) { // for each document
      
      IntegerVector doc = docs[d]; // placeholder for a document
      
      IntegerVector z_dn_row = z_dn[d]; // placeholder for doc-word-topic assigment
      
      for (int n = 0; n < doc.length(); n++) { // for each word in that document
        
        // discount for the n-th word with topic z
        z2 = z_dn_row[n];
        
        theta_counts(d,z2) = theta_counts(d,z2) - 1; 
        
        phi_counts(z2,doc[n]) = phi_counts(z2,doc[n]) - 1;
        
        n_z[z2] = n_z[z2] - 1;
        
        n_d[d] = n_d[d] - 1;
        
        // sample topic index
        for (int k = 0; k < Nk; k++) {
          
          p_z[k] = (phi_counts(k,doc[n]) + beta[doc[n]]) / (n_z[k] + k_beta[doc[n]]) *
            (theta_counts(d,k) + alpha[k]) / (n_d[d] + sum_alpha);
          
        }
        
        
        // update counts
        z1 = RcppArmadillo::sample(topic_sample, 1, false, p_z);
        
        z2 = z1[1];
        
        theta_counts(d,z2) = theta_counts(d,z2) + 1; // update document topic count
        
        phi_counts(z2,doc[n]) = phi_counts(z2,doc[n]) + 1; // update topic word count
        
        n_d[d] = n_d[d] + 1; // count that topic in that document overall
        
        n_z[z2] = n_z[z2] + 1; // count the that topic overall
        
        z_dn_row[n] = z2; // # count that topic for that word in the document
        
      }
    }
    
    // if using burnin, update sums
    if (burnin > -1 && i >= burnin) {
      
      for (int k = 0; k < Nk; k++) {
        
        for (int d = 0; d < Nd; d++) {
          theta_sums(d,k) = theta_sums(d,k) + theta_counts(d,k);
        }
        
        for (int v = 0; v < Nv; v++) {
          phi_sums(k,v) = phi_sums(k,v) + phi_counts(k,v);
        }
      }
    }
    
  }
  
  // return the result

  if (burnin > -1) {
    int i_diff = iterations - burnin;
    
    NumericMatrix theta(Nd,Nk);
    
    NumericMatrix phi(Nk,Nv);
    
    // average over chain after burnin 
    for (int k = 0; k < Nk; k++) {
      
      for (int d = 0; d < Nd; d++) {
        theta(d,k) = (theta_sums(d,k) / i_diff);
      }
      
      for (int v; v < Nv; v++) {
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


