#include "iterativeSEM.hpp"
#include <algorithm>
using namespace std;

void scramble_kmer(Dataset &data){
    for(auto it = scramble_kmers.begin(); it != scramble_kmers.end(); ++it){
        random_shuffle(it->begin(), it->end());
    }
}
