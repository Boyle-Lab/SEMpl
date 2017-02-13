#include "iterativeSEM.hpp"
#include <iostream>
#include "../lib/bowtie-1.0.0/ebwt.h"
using namespace std;

// from bowtie
/*static void driver(const char * type, // "DNA", by Cody
                   const string& ebwtFileBase,
                   const string& query,
                   const vector<string>& queries,
                   const vector<string>& qualities,
                   const string& outfile);*/
void split(std::string str, std::string splitBy, std::vector<std::string>& tokens);
string revCompDNA(string dna);
// if I understand correctly, bowtie places output in the filename specified
// at the end of the command                              // intermediate.dat, as
//                                                        // defined in iterativeSEM.hpp
// genome is "./data/hg19", DNA_FA_FILE  is the file that contains .fa
void bowtie_genome_map(Dataset &data, int length, const string& genome){
    char *argvs[9];

    argvs[0] = "./bin/bowtie\0";
    argvs[1] = "--quiet\0";
    argvs[2] = "-a\0";
    argvs[3] = "-v\0";
    argvs[4] = "0\0";
    argvs[5]  = genome;
    argvs[6] = "-f\0";
    argvs[7] = "../data/hg19\0";
    argvs[8] = "temp.dat\0";

    try{
        // bowtie throws exceptions, I believe
        bowtie(9, argvs);
    }
    catch(exception e){
        cerr << e.what() << '\n';
        exit(1);
    }
    catch(...){
        cerr << "Error with bowtie\n";
        exit(1);
    }

    ifstream IN("temp.dat");
    if(!IN){
        cerr << "Failure to open temp.dat, ephemeral file\n";
        exit(1);
    }
    string map_line = "";
    vector<string> dat;
    string DNA = "";
    while(getline(IN, map)){
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
        //                                            convert to int, perform addition, convert back to string
        data.bowtie_output.push_back({dat[2],
                                      dat[3],
                                      to_string(atoi(dat[3].c_str() +  length) ),
                                      DNA,
                                      dat[1]});
    }
    IN.close();
}
