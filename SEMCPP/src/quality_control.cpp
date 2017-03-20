#include "iterativeSEM.hpp"
#include "common.h"
#include <fstream>
#include <sstream>


/*
*
*   Effect: Applies quality control measures to SEM output files.
*
*/

using namespace std;

int count_kmer(Dataset &data);
int findMaximumPerRow(Dataset &data, const string dest);
int ttest(Dataset &data, int max_scram, int max_enum);

void quality_control(Dataset &data){

    ofstream quality_output;
    stringstream output_stream;
    output_stream << data.output_dir << "/Quality_control.txt";
    quality_output.open(output_stream.str());

    if(!quality_output){

    }

    //Step 1: Total kmer check
    int total_kmers = count_kmer(data);
    quality_output << "Total k-mer count: " << total_kmers <<"\n";

    //Step 2: T-test on signal to background

    Dataset::accumSummaryData::accumSummary_dest dest = 
                            Dataset::accumSummaryData::accumSummary_dest::none;

    //First generate baseline value file
    dest = Dataset::accumSummaryData::accumSummary_dest::scrambled;

    //int max_scram = findMaximumPerRow(data, dest);

    // will operate on final data from creating baseline from initially
    // scrambled kmers

    //Next generate signal value file
    dest = Dataset::accumSummaryData::accumSummary_dest::enumerated;

    //int max_enum = findMaximumPerRow(data, dest);

    // will operate on final data from creating baseline from initially
    // enumerated kmers

    //Then run t-test
    double t_result = ttest(data, max_scram, max_enum);
    quality_output << "Signal to background T-test: "
                   << t_result << "\n";
    quality_output.close();

}

int count_kmer(Dataset &data){

    //Must use the output from enumerate kmer and count the kmers

#ifdef DEBUG
    if(data.size_of_kmerHash != data.kmerHash.size()){
        cout << "size of kmerHash should not have changed!!!!!!\n"
             << "exiting!!!\n";
             exit(1);
    }
#endif

    return static_cast<int>(data.size_of_kmerHash);

}

int findMaximumPerRow(Dataset &data, 
                      Dataset::accumSummaryData::accumSummary_dest dest){

    //Uses output from accum_Summary to do calculations
    //to ensure no large dead regions

        

}

double ttest(Dataset &data, int max_scram, int max_enum){

    //Interacts with R to perform a T-test

}
