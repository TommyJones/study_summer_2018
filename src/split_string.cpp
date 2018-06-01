#include <Rcpp.h>
using namespace Rcpp;
using std::string;

// This is a simple example of splitting a string on spaces

// [[Rcpp::export]]
std::vector<std::string> split_string(const std::string str, const std::string delim = " "){
  
  // initialize the output result
  std::vector<std::string> tokens;
  
  // Skip delim at beginning.
  string::size_type last_pos = str.find_first_not_of(delim, 0);
  
  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delim, last_pos);
  
  while (string::npos != pos || string::npos != last_pos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(last_pos, pos - last_pos));
    
    // Skip delim.  Note the "not_of"
    last_pos = str.find_first_not_of(delim, pos);
    
    // Find next "non-delimiter"
    pos = str.find_first_of(delim, last_pos);
  }
  
  return tokens;
};

