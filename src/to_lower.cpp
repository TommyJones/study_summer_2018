#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP to_lower(const std::string s,
                     Function f){
  
  // don't overwrite inputs
  SEXP str = f(s);
  
  return str;
}

// std::string to_lower(const std::string s,
//                      const std::string letters,
//                      const std::string LETTERS){
//   
//   // don't overwrite inputs
//   std::string str = s;
//   
//   // loop over each element and replace if uppercase
//   for(int j = 0; j < str.length(); j++){
//     for(int k = 0; k < letters.length(); k++){
//       if(str[j] == LETTERS[k]){
//         str[j] = letters[k];
//       }
//     }
//   }
//   return str;
// }

// StringVector to_lower(const StringVector s, 
//                                   const StringVector letters, 
//                                   const StringVector LETTERS){
//   
//   // don't overwrite inputs
//   StringVector str = s;
//   
//   // loop over each element of the string and replace if uppercase
//   for(int j = 0; j < str.size(); j++){
//     for(int k = 0; k < letters.size(); k++){
//       std::string a = str[j];
//       std::string b = LETTERS[k];
//       
//       if(str[j] == LETTERS[k]){
//         str[j] = letters[k];
//       }
//     }
//   }
//   
//   return str;
//   
// }

