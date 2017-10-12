#include "src/iterativeSEM.hpp"
#ifdef DEBUG
    #include <iostream>
#endif
using namespace std;

// EFFECTS: fills new_kmer_vec, using data.kmerHash
void changeBase(const Dataset &data, int position, const char nucleotide, vector<string> &new_kmer_vec, const string &genome){
  // find new_kmer
  // performs operations on Enumerated_kmer, or kmerHash
  // remove \r from line, not needed here, as stored in memory
    new_kmer_vec.clear();
    new_kmer_vec.reserve(data.kmerHash.size());

    string new_kmer = "";

    for(auto pair : data.kmerHash){
        if(position > static_cast<int>(pair.first.size()) - 1){
            cout << "Position is greater than size!!\n\tEXITING";
            exit(1);
        }
        new_kmer = pair.first;
        new_kmer[position] = nucleotide;
        new_kmer_vec.push_back(new_kmer);
    }
}
