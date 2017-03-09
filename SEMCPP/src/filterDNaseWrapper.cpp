#include "iterativeSEM.hpp"
#include "common.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std;

static void read_in_files(Dataset &data, const string &file);

// REQUIRES: Boost library installed
// MEMORY: incrementally reads in output
void filterDNaseWrapper(Dataset &data){
   /*
    * I have added a member variable for the output of this function
    * The type is std::vector<std::pair<std::string, std::string > >
    * and the name of the variable is filterDNaseWrapper_output
    *
    */


    //Working on the implementation of the bedtools library

    //bedtools(Dataset &data);

    //See bedtools.txt for notes on how it works and an alternative for using
    //the file in order to work with struct data


    /*
    *   I will make a quick implementation of this function, involving
    *   reading in and out, just so that the function is working for now
    *       -- Cody
    */

    //
    vector<string> files;
    string targetDir = "../" + data.output_dir + "/ALIGNMENT/";

    GetFilesInDirectory(files, "../ALIGNMENT/");
    for(auto file : files){
        if(file.find("pos") != string::npos){
            string bedfile = targetDir + file + "_filtered.bed";
            string readfile = targetDir + file + ".bed";
            ifstream IN(file);

            if(fileExists(readfile)){
                // CALL BEDTOOLS AS PRESCRIBED IN ORIGINAL ALGORITHM
#ifdef DEBUG
                assert(!data.DNase_file.empty());
                assert(!targetDir.empty());
                assert(!file.empty());
#endif


                string bed_cmd = "../lib/bedtools intersect -a " + readfile
                    + " -b " + data.DNase_file + " -wa -u | sort | uniq > "
                    + bedfile;
#ifdef DEBUG
                int val = system(bed_cmd.c_str());
                assert(val == 0);
#else
                system(bed_cmd.c_str());
#endif
                //I believe this is the correct implementation of the call
                // to bedtools with the binary file
                // located within lib
            }
            else{
                system(string("touch " + bedfile).c_str());
            }
#ifdef DEBUG
            cout << "\tdeleting readfile\n";
#endif
            system(string("rm -f " + readfile).c_str());

            read_in_files(data, bedfile);
        }

    }
}


static void read_in_files(Dataset &data, const string &file){
    ifstream fin(file);

    if(!fin){
        cerr << "unable to open " << file << " for read in\n";
        exit(1);
    }

    string line = "";

    while(getline(fin, line, '\n')){
        data.filterDNaseWrapper_output.push_back(line);
    }
}
