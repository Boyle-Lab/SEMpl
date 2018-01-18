#include "iterativeSEM.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>
#include <sstream>
using namespace std;


// if I understand correctly, bowtie places output in the filename specified
// at the end of the command                              // intermediate.dat, as
//                                                        // defined in iterativeSEM.hpp
// genome is "./data/hg19", DNA_FA_FILE  is the file that contains .fa


//int bowtie(int argc, const char **argv);


void bowtie_genome_map(int length, const string& genome, const string& file, 
                       const string& final_output, const string& dnase_file, bool verbose){
    //const char *argvs[9] = {"./bin/bowtie", "--quiet", "-a", "-v 0", 
    //                        genome.c_str(), "-f", file.c_str(), "temp.dat" };

//    string cmd = "./bin/bowtie --threads 10 --quiet -a -v 0 ./data/hg19 -f " + file;
    string cmd = "./bin/bowtie --threads 10 --quiet -a -v 0 ./data/hg19 -f " + file 
       + " | awk '{print $3\"\\t\"$4\"\\t\"$4+" + std::to_string(length) + "\"\\t\"$5\"\\t\"$2}' | " 
       + "./bin/bedtools intersect -a stdin -b " + dnase_file + " -wa -u | sort | uniq";

    if(verbose){
        cout << "Running command: " << cmd << "\n\tRunning...." << flush;
    }


    ofstream OUT(final_output);
    if(!OUT){
        cerr << "Failure to open \"" << final_output << "\".\n";
        exit(1);
    }


    std::stringstream IN = exec(cmd.c_str());
    string map_line = "";
    vector<string> dat;
    string DNA = "";
    while(std::getline(IN, map_line, '\n')){
        #ifdef DEBUG
           //cout << map_line << endl << flush;
        #endif
        split_string(map_line, "\t", dat);

/*
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

        OUT << dat[2] << '\t' << dat[3] << '\t'
            << atoi(dat[3].c_str()) + length
            << '\t' << DNA << '\t' << dat[1] << '\n';
*/
        if(dat[4] == "+"){
            DNA = dat[3];
        }
        else if(dat[4] == "-"){
            DNA = revCompDNA(dat[3]);
        }
        else{
            cerr << "unknown strand!\n";
            cerr << map_line << endl;
            exit(1);
        }

        OUT << dat[0] << '\t' << dat[1] << '\t'
            << dat[2] << '\t' << DNA << '\t' << dat[4] << '\n';

        dat.clear();
    }

    if(verbose){
        cout << "FINISH" << endl;
    }

    OUT.close();

}

