#include "iterativeSEM.hpp"
#include "common.hpp"
#include <fstream>
#include <sstream>
#include <cstdio>


/*
*
*   Effect: Applies quality control measures to SEM output files.
*
*/

using namespace std;

static int count_kmer(const Dataset &data);
//int findMaximumPerRow(Dataset &data, const string dest);
double ttest(const Dataset &data);

// REQUIRES: data.signal_scramble_output is filled
//           along with data.signal_enumerate_output
void quality_control(Dataset &data){

    ofstream quality_output;
    stringstream output_stream;
    output_stream << data.output_dir << "/Quality_control.txt";
    quality_output.open(output_stream.str());

    

    //Step 1: Total kmer check
    int total_kmers = count_kmer(data);
    quality_output << "Total k-mer count: " << total_kmers <<"\n";

    //Step 2: T-test on signal to background
    double p_val = ttest(data);
    quality_output << "Signal to background T-test: " << p_val << '\n';


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

// THIS MIGHT BE DEPRECATED BELOW
/*
int findMaximumPerRow(Dataset &data, 
                      Dataset::accumSummary_type::accumSummary_dest dest){

    //Uses output from accum_Summary to do calculations
    //to ensure no large dead regions

        

}
*/

double ttest(const Dataset &data){

    //Interacts with R to perform a T-test

    // writes the three files necessary to do the t-test
    ofstream OUT(data.output_dir + "/runTtest.R");
    ofstream OUT1(data.output_dir + "Enumerated_kmer_filtered.signal");
    ofstream OUT2(data.output_dir + "Scrambled_kmer.filtered.signal");
    for(auto val : data.signal_scramble_output){
        OUT2 << val << '\n';
    }
    for(auto val : data.signal_enumerate_output){
        OUT1 << val << '\n';
    }
    OUT << "signal <- read.table(\"" <<  data.output_dir 
        << "Enumerated_kmer_filtered.signal\")" << '\n'
        << "baseline <- read.table(\"" <<  data.output_dir 
        << "Scrambled_kmer_filtered.signal\")" << "\n\n"
        << "res <- t.test(signal$V6, baseline$V6)\n"
        << "pval <- -log10(res$p.value)\n\n"
        << "cat(pval, sep=\"\\n\")\n";

    string cmd = "R --vanilla --no-save --slave < " 
                 + data.output_dir 
                 + "/runTtest.R";

    FILE * strm = popen(cmd.c_str(), "r");
    double p_val = -1.0;
    int message = fscanf(strm, " %lf", &p_val);
    if(p_val == -1.0){
        cerr << "Failure to read p_val!!!!!!\n\tEXITING";
        exit(1);
    }
    if(message != 1){
        cerr << "incorrect message from fscanf(args)!!!!!\n\tEXITING";
        exit(1);
    }
    return p_val;
}
