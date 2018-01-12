#include "iterativeSEM.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>
#include <sstream>
using namespace std;


//From https://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix
std::stringstream exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::stringstream result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result << buffer.data();
    }
    return result;
}


// if I understand correctly, bowtie places output in the filename specified
// at the end of the command                              // intermediate.dat, as
//                                                        // defined in iterativeSEM.hpp
// genome is "./data/hg19", DNA_FA_FILE  is the file that contains .fa


//int bowtie(int argc, const char **argv);


void bowtie_genome_map(int length, const string& genome, const string& file, 
                       const string& final_output, bool verbose){
    //const char *argvs[9] = {"./bin/bowtie", "--quiet", "-a", "-v 0", 
    //                        genome.c_str(), "-f", file.c_str(), "temp.dat" };

    string cmd = "./bin/bowtie --threads 5 --quiet -a -v 0 ./data/hg19 -f " + file;

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

    if(verbose){
        cout << "FINISH" << endl;
    }


}

