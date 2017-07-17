#include "iterativeSEM.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>
//#include "../lib/bowtie-1.0.0/ebwt.h"
using namespace std;


// from bowtie
/*static void driver(const char * type, // "DNA", by Cody
                   const string& ebwtFileBase,
                   const string& query,
                   const vector<string>& queries,
                   const vector<string>& qualities,
                   const string& outfile);*/


// if I understand correctly, bowtie places output in the filename specified
// at the end of the command                              // intermediate.dat, as
//                                                        // defined in iterativeSEM.hpp
// genome is "./data/hg19", DNA_FA_FILE  is the file that contains .fa


//int bowtie(int argc, const char **argv);


void bowtie_genome_map(int length, const string& genome, const string& file, 
                       const string& final_output, bool verbose){
    //const char *argvs[9] = {"./bin/bowtie", "--quiet", "-a", "-v 0", 
    //                        genome.c_str(), "-f", file.c_str(), "temp.dat" };
/*
*    argvs[0] = "./bin/bowtie\0";
*    argvs[1] = "--quiet\0";
*    argvs[2] = "-a\0";
*    argvs[3] = "-v\0";
*    argvs[4] = "0\0";
*    argvs[5]  = genome;
*    argvs[6] = "-f\0";
*    argvs[7] = "../data/hg19\0";
*    argvs[8] = "temp.dat\0";
*/
    const string temp_file = "./results/temp.dat";

    string cmd = "./bin/bowtie --quiet -a -v 0 ./data/hg19 -f " 
                 + file + ' ' + temp_file;
    
    
    if(verbose){
        cout << "Running command: " << cmd << "\n\tRunning...." << flush;
    }
    
    if(system(cmd.c_str()) != 0){
        cerr << "problem running " << cmd << endl;
        exit(1);
    }

    if(verbose){
        cout << "FINISH" << endl;
    }


    ofstream OUT(final_output);
    if(!OUT){
        cerr << "Failure to open \"" << final_output << "\".\n";
        exit(1);
    }

    ifstream IN(temp_file);
    if(!IN){
        cerr << "Failure to open \"" << file << "\".\n";
        exit(1);
    }

    string map_line = "";
    vector<string> dat;
    string DNA = "";
    while(getline(IN, map_line)){
        #ifdef DEBUG
            // cout << map_line << endl;
        #endif
        split_string(map_line, "\t", dat);
        if(dat[1] == "+"){
            DNA = dat[4];
        }
        else if(dat[1] == "-"){
            #ifdef DEBUG
            // cout << "kmer:     " << dat[4] << endl;
            #endif
            DNA = revCompDNA(dat[4]);
            #ifdef DEBUG
            // cout << "reversed: " << DNA << endl;
            #endif

        }
        else{
            cerr << "unknown strand!\n";
            cerr << map_line << endl;
            exit(1);
        }

        OUT << dat[2] << '\t' << dat[3] << '\t'
            << atoi(dat[3].c_str()) + length
            << '\t' << DNA << '\t' << dat[1] << '\n';

        dat.clear();
    }
    IN.close();
}
