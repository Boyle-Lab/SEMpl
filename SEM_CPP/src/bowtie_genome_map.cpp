#include "iterativeSEM.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>
#include <sstream>
using namespace std;


// Takes fasta file and aligns to genome, filters by DNase, and sorts
//  stores results in final output
// bowtie can take input as - to be stdin (possible direct pipe from kmer enumeration)
void bowtie_genome_map(int length, const string& genome, const string& file, 
                       const string& final_output, const string& dnase_file, int threads, bool verbose){

    const string temp_file = "temp_bowtie.dat";
    string cmd = "./bin/bowtie --threads " + std::to_string(threads) + " --quiet -a -v 0 ./data/hg19 -r " + file + " > /dev/null";

    if(verbose){
        cout << "Running command: " << cmd << "\n\tRunning...." << flush;
    }

//    ofstream OUT(temp_file);
//    if(!OUT){
//        cerr << "Failure to open \"" << temp_file << "\".\n";
//        exit(1);
//    }

    exec(cmd.c_str());
//    std::stringstream IN = exec(cmd.c_str()); // stringstream instead of writing to a file
//    string map_line = "";
//    vector<string> dat;
//    string DNA = "";
//    int DNAlen;

/*
    // Read output, reformat, and write to temp file
    // This could all be kept in memory in the future.
    while(std::getline(IN, map_line, '\n')){
        split_string(map_line, "\t", dat);

        if(dat[1] == "+"){
            DNA = dat[4];
        }
        else if(dat[1] == "-"){
            DNA = revCompDNA(dat[4]);
        }
        else{
            cerr << "unknown strand!\n";
            cerr << map_line << endl;
            exit(1);
        }

//        if(std::stoi(dat[6]) > 0) {  // bowtie sometimes returns alignments with first base missing
        DNAlen = DNA.length();
        if(DNAlen == length) {
          OUT << dat[2] << '\t' << dat[3] << '\t'
              << std::stoi(dat[3]) + length << '\t' << DNA << '\t' << dat[1] << '\n';
        }

        dat.clear();
    }
*/
    if(verbose){
        cout << "Bowtie complete..." << endl;
    }

//    OUT.close();
//    IN.str(std::string());

    cmd = "./bin/bedtools intersect -a " + temp_file + " -b " + dnase_file + " -wa -u | sort | uniq > " + final_output;

    if(verbose){
        cout << "Running command: " << cmd << "\n\tRunning...." << flush;
    }

    exec(cmd.c_str()); // stringstream instead of writing to a file

    if(verbose){
        cout << "FINISH" << endl;
    }

}

int main(int argc, char* argv[])
{
    
