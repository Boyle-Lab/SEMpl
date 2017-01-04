#include <string>
#include <regex>
using namespace std;

string revCompDNA(string dna){
  string rev = "";
  for(int i = static_cast<int>(dna.size()) - 1; i >= 0; i++){
    rev += dna[i];
  }
  rev = regex_replace(dna, regex("ACGTacgt"), "TGCAtgca");
  return rev;
}

// int parse_wc{
//   string s;
//   while(getline(s, ))
// }
