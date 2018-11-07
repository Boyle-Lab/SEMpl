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
string ttest(const Dataset &data);
static void generate_input(const Dataset &data);
string run_R(const Dataset &data);

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
    string p_val = ttest(data);
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

string ttest(const Dataset &data){

    generate_input(data);
    return run_R(data);
}

static void generate_input(const Dataset &data){

    stringstream Rout;
    Rout << data.output_dir  << "/ttest.R";
    ofstream Rfile;
    Rfile.open(Rout.str());

    Rfile << "signal <- read.table(\"" << data.output_dir << "/BASELINE/Enumerated_kmer.signal\")\n"
          << "baseline <- read.table(\"" << data.output_dir << "/BASELINE/Scrambled_kmer.signal\")\n"
          << "res <- t.test(signal$V6[which(signal$V6>-255)], baseline$V6[which(baseline$V6>-255)])\n"
          << "pval <- -log10(res$p.value)\n"
//          << "if(is.infinite(pval)) { pval<-310 }\n"
          << "cat(pval, sep=\"\\n\")";
    Rfile.close();
}

string run_R(const Dataset &data){
    string s = "R --vanilla --slave < " + data.output_dir + "/ttest.R";
    std::stringstream IN = exec(s.c_str()); // stringstream instead of writing to a file
    string result = "";
    std::getline(IN, result, '\n');

    string s = "rm " + data.output_dir + "/ttest.R";
    system(s.c_str());

    return result;
}
