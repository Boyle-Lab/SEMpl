#include "iterativeSEM.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <unistd.h>
#include <sys/wait.h>
using namespace std;


// Takes fasta file and aligns to genome, filters by DNase, and sorts
//  stores results in final output
// bowtie can take input as - to be stdin (possible direct pipe from kmer enumeration)
void bowtie_genome_map(int length, const string& genome, const string& file, 
                       const string& final_output, const string& dnase_file, int threads, bool verbose){

    string cmd = "./map-filer " + std::to_string(length) + " " + genome + " " + file + " " + final_output + " " + dnase_file + " " + std::to_string(threads) + " " + std::to_string(verbose);

    if(verbose){
        cout << "Running command: " << cmd << "\n\tRunning...." << flush;
    }

    bool bad_run = false;
    pid_t pid;
    int child_status;
    do {
        pid = fork();
        if(pid == 0) { /* Child */
            exec(cmd.c_str());
            exit(0);
        } else { /* Parent */
            wait(&child_status);
            if(child_status != 0) {
                cerr << "Status: " << child_status << " Need to run the process again!\n" << flush;
                bad_run = true;
            }
        }
    } while(bad_run);

    if(verbose){
        cout << "FINISH" << endl;
    }

}

