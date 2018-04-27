#include "iterativeSEM.hpp"
#include "common.hpp"
#include "t-test.hpp"
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

double ttest_internal(const std::vector<double> &one, const std::vector<double> &two){
    
    double p_val = 0.0;
    double signal_mean = 0.0, baseline_mean = 0.0;
    double signal_stdev = 0.0, baseline_stdev = 0.0;
    unsigned signal_size = 0, baseline_size = 0.0;

    signal_size = one.size();
    baseline_size = two.size();

    for(const auto& value : one){
        signal_mean += value;
    }
    signal_mean /= signal_size;
    for(const auto& value : two){
        baseline_mean += value;
    }
    baseline_mean /= baseline_size;

    double temp = 0.0;
    for(const auto& value : one){
        temp += (value - signal_mean) * (value - signal_mean);
    }
    temp /= signal_size - 1;
    signal_stdev = pow(temp, 0.5);


    temp = 0.0;
    for(const auto& value : two){
        temp += (value - baseline_mean) * (value - baseline_mean);
    }
    temp /= baseline_size - 1;
    baseline_stdev = pow(temp, 0.5);

    
    p_val = two_samples_t_test_unequal_sd(  signal_mean,   // Sm1 = Sample Mean 1.
                                            signal_stdev,   // Sd1 = Sample Standard Deviation 1.
                                            signal_size,   // Sn1 = Sample Size 1.
                                            baseline_mean,   // Sm2 = Sample Mean 2.
                                            baseline_stdev,   // Sd2 = Sample Standard Deviation 2.
                                            baseline_size);   // Sn2 = Sample Size 2.
    return p_val;
}

double ttest(const Dataset &data){

    //Interacts with R to perform a T-test

    // Perform's t-test between signal (enumerated kmers)
    // and baseline (scrambled kmers)
   
    // performs type conversion to double
    std::vector<double> grab_score_vec_one, grab_score_vec_two;
    std::string s = "";
    for(const auto& value : data.signal_scramble_output){
       s = grab_string_last_index(value);
       grab_score_vec_two.push_back(stod(s));
    }
    for(const auto& value : data.signal_enumerate_output){
       s = grab_string_last_index(value);
       grab_score_vec_one.push_back(stod(s));
    }

    return ttest_internal(grab_score_vec_one, grab_score_vec_two);
}


