#include "iterativeSEM.hpp"
#include <algorithm>
using namespace std;

void scramble_kmer(Dataset &data){
    for(auto it = data.scramble_kmers.begin(); it != data.scramble_kmers.end(); ++it){
        random_shuffle(it->begin(), it->end());
    }
}
