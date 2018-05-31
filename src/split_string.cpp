#include <Rcpp.h>
using namespace Rcpp;
using std::string;
using std::strtok;

// [[Rcpp::plugins(cpp11)]]
// This is a simple example of splitting a string on spaces

// [[Rcpp::export]]

Rcpp::RawVector split_string(const std::string &s, const char* delim, std::vector<std::string> & v){
  // to avoid modifying original string
  // first duplicate the original string and return a char pointer then free the memory
  char * dup = strdup(s.c_str());
  char * token = strtok(dup, delim);
  while(token != NULL){
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
  
  return token;
}



// std::vector<std::string> split_string(std::string x) {
//   std::string s = x; // don't overwrite inputs
//   
//   char token ;
//   
//   std::vector<std::string> result {
//     std::regex_token_iterator(s.begin(), s.end(), ws_re, -1), {}
//   };  
//   
//   return result;
// }



