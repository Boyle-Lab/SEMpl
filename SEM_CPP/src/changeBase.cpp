#include "src/iterativeSEM.hpp"
#include "src/common.hpp"
#ifdef DEBUG
    #include <iostream>
#endif
using namespace std;

// EFFECTS: fills new_kmer_vec, using data.kmerHash
void changeBase(const Dataset &data, int position, const char nucleotide, vector<string> &new_kmer_vec){
  // find new_kmer
  // performs operations on Enumerated_kmer, or kmerHash
  // remove \r from line, not needed here, as stored in memory
    new_kmer_vec.clear();
    new_kmer_vec.resize(data.kmerHash.size());

    vector<string> new_kmer_unsort;
    new_kmer_unsort.reserve(data.kmerHash.size());

    string new_kmer = "";

    for(auto pair : data.kmerHash){
        if(position > static_cast<int>(pair.first.size()) - 1){
            cout << "Position is greater than size!!\n\tEXITING";
            exit(1);
        }
        new_kmer = pair.first;
        new_kmer[position] = nucleotide;
        new_kmer_unsort.push_back(new_kmer);
    }

    //Make these unique
    sort(new_kmer_unsort.begin(), new_kmer_unsort.end());
    auto iter2 = unique_copy(new_kmer_unsort.begin(), new_kmer_unsort.end(), new_kmer_vec.begin());
    new_kmer_vec.resize(iter2 - new_kmer_vec.begin());

}
