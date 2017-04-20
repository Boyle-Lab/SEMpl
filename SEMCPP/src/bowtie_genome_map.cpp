#include "iterativeSEM.hpp"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "src/common.hpp"
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


void bowtie_genome_map(int length, const string& genome, const string& file, const string& final_output){
    //const char *argvs[9] = {"./bin/bowtie", "--quiet", "-a", "-v 0", genome.c_str(), "-f", file.c_str(), "temp.dat" };
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

    stringstream cmd;
    cmd << "./bin/bowtie --quiet -a -v 0 " << genome << " -f " << file;
    system(cmd.str().c_str());

//    try{
//        // bowtie throws exceptions
//    }
//    catch(exception e){
//        cerr << e.what() << '\n';
//        exit(1);
//    }
//    catch(...){
//        cerr << "Error with bowtie\n";
//        exit(1);
//    }


    ofstream OUT(final_output);
    if(!OUT){
        cerr << "Failure to open \"" << final_output << "\".\n";
        exit(1);
    }

    ifstream IN(file);
    if(!IN){
        cerr << "Failure to open \"" << file << "\".\n";
        exit(1);
    }

    string map_line = "";
    vector<string> dat;
    string DNA = "";
    while(getline(IN, map_line)){
        split(map_line, "\t", dat);
        if(dat[1] == "+"){
            DNA = dat[4];
        }
        else if(dat[1] == "-"){
            DNA = revCompDNA(dat[4]);
        }
        else{
            cerr << "unknown strand!\n";
            cerr << map_line << '\n';
            exit(1);
        }
        // bowtie_output is a vector of arrays, each array has 5 string spots
        //              convert to int, perform addition, convert back to string
        //data.bowtie_output.push_back({{dat[2], dat[3],
        //                              to_string(atoi(dat[3].c_str() +  length) ),
        //                                  DNA, dat[1]}});

        OUT << dat[2] << '\t' << dat[3] << '\t'
            << atoi(dat[3].c_str()) + length
            << '\t' << DNA << '\t' << dat[1] << '\n';
    }
    IN.close();
}
