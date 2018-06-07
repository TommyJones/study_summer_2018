#include <Rcpp.h>
using namespace Rcpp;


// This function takes an R character vector as input and returns a list with
// each element a vector tokenized on words

//[[Rcpp::export]]
List split_string_v(const std::vector<std::string> x, 
                    const std::string delim = " "){
  
  // initialize output
  List result(x.size());
  
  for(int j = 0; j < x.size(); j++){
    // don't overwrite inputs
    std::string str = x[j];
    
    // initialize the output result
    std::vector<std::string> tokens;
    
    // Skip delim at beginning.
    std::string::size_type last_pos = str.find_first_not_of(delim, 0);
    
    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delim, last_pos);
    
    while (std::string::npos != pos || std::string::npos != last_pos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(last_pos, pos - last_pos));
      
      // Skip delim.  Note the "not_of"
      last_pos = str.find_first_not_of(delim, pos);
      
      // Find next "non-delimiter"
      pos = str.find_first_of(delim, last_pos);
    }
    result[j] = tokens;
  }
  return result;
}
