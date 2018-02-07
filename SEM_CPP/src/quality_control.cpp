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
void quality_control(const Dataset &data){

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

// This doesn't work - can we keep this in c++?
double ttest(const Dataset &data){

    //Interacts with R to perform a T-test

    // writes the three files necessary to do the t-test
    ofstream OUT(data.output_dir + "/runTtest.R");

    OUT << "signal <- read.table(\"" << data.output_dir 
        << "/BASELINE/Enumerated_kmer.signal\")" << '\n'
        << "baseline <- read.table(\"" <<  data.output_dir 
        << "/BASELINE/Scrambled_kmer.signal\")" << "\n\n"
        << "res <- t.test(signal$V6, baseline$V6)\n"
        << "pval <- -log10(res$p.value)\n\n"
        << "cat(pval, sep=\"\\n\")\n";


    OUT.close();

    string cmd = "R --vanilla --no-save --slave < " 
                 + data.output_dir 
                 + "/runTtest.R > " + data.output_dir + "/ttest.txt";

    if(system(cmd.c_str() ) ){
        cerr << "problem running " << cmd << endl;
//        exit(1);
    }
    double p_val = -1.0;

    // int message = fscanf(strm, " %lf", &p_val);
    if(p_val == -1.0){
        cerr << "Failure to read p_val!!!!!!\n\tEXITING";
//        exit(1);
    }

    return p_val;
}
