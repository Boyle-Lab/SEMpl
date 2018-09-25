#include "iterativeSEM.hpp"
#include "common.hpp"
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>


/*
*
*   Effect: Applies quality control measures to SEM output files.
*
*/

using namespace std;

static int count_kmer(const Dataset &data);
//int findMaximumPerRow(Dataset &data, const string dest);

// REQUIRES: data.signal_scramble_output is filled
//           along with data.signal_enumerate_output
void quality_control(const Dataset &data){

    ofstream quality_output;
    stringstream output_stream;
    output_stream << data.output_dir << "/Quality_control.txt";
    quality_output.open(output_stream.str());

    //Step 1: Total kmer check
    int total_kmers = count_kmer(data);
    quality_output << "Total k-mer count: " << total_kmers <<"\n";

    //Step 2: T-test on signal to background
//    double p_val = ttest(data);
//    quality_output << "Signal to background T-test: " << p_val << '\n';
}

int count_kmer(const Dataset &data){

    //Must use the output from enumerate kmer and count the kmers

#ifdef DEBUG
    if(data.size_of_kmerHash != data.kmerHash.size()){
        cout << "size of kmerHash should not have changed!!!!!!\n"
             << "exiting!!!\n";
             exit(1);
    }
#endif

    return static_cast<int>(data.kmerHash.size());

}

