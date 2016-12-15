#include <string>
#include <vector>
#include "../iterativeSEM.hpp"
using namespace std;

int getLength(Dataset &data);

void alignToGenomeWrapper(Dataset &data, int iteration) {

  string genome = "hg19";
  vector<string> nucleotideStack(4);
  nucleotideStack.push_back("A");
  nucleotideStack.push_back("C");
  nucleotideStack.push_back("T");
  nucleotideStack.push_back("G");

  // get the length of kmer
  //int length = getLength(data);
// UNFINISHED
// UNFINISHED


}


// getlength accesses first(?) element of kmerHash and returns the key
int getLength(Dataset &data){
  return data.kmerHash.begin()->second;
}
