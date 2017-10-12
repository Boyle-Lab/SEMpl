#include "iterativeSEM.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>
using namespace std;




// if I understand correctly, bowtie places output in the filename specified
// at the end of the command                              // intermediate.dat, as
//                                                        // defined in iterativeSEM.hpp
// genome is "./data/hg19", DNA_FA_FILE  is the file that contains .fa


//int bowtie(int argc, const char **argv);


void bowtie_genome_map(int length, const string& genome, const string& file, 
                       const string& final_output, bool verbose){
    //const char *argvs[9] = {"./bin/bowtie", "--quiet", "-a", "-v 0", 
    //                        genome.c_str(), "-f", file.c_str(), "temp.dat" };

    const string temp_file = "temp_bowtie.dat";

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

        dat.clear();
    }
    IN.close();
}
