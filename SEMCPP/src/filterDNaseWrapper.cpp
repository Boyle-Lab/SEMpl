#include "iterativeSEM.hpp"
#include "common.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std;

// static void read_in_files(Dataset &data, const string &file);

// REQUIRES: Boost library installed
// MEMORY: incrementally reads in output
// UNSURE OF ABOVE LINE UPDATE
void filterDNaseWrapper(const Dataset &data){

    vector<string> files;
    string targetDir = "./" + data.output_dir + "/ALIGNMENT/";

    #ifdef DEBUG
        // cout << "targetDir: " << targetDir << endl;
    #endif

    GetFilesInDirectory(files, targetDir);
#ifdef DEBUG
    sort(files.begin(), files.end());
#endif
    for(auto file : files){
        if(file.find("pos") != string::npos 
            && file.find("filtered") == string::npos
            && file.find("fa") == string::npos){

            

            string bedfile = /* targetDir + */ file + "_filtered";
            string &readfile = /* targetDir + */ file /* + ".bed" */;

            #ifdef DEBUG
                // cout << "\tfile: " << file << endl;
                // cout << "\tbedfile: " << bedfile << endl;
                // cout << "\treadfile: " << readfile << endl;
            #endif

            ifstream IN(file);

            if(fileExists(readfile)){
                // CALL BEDTOOLS AS PRESCRIBED IN ORIGINAL ALGORITHM
#ifdef DEBUG
                assert(!data.DNase_file.empty());
                assert(!targetDir.empty());
                assert(!file.empty());
#endif


                string bed_cmd = "./bin/bedtools intersect -a " + readfile
                    + " -b " + data.DNase_file + " -wa -u | sort | uniq > "
                    + bedfile;
#ifdef DEBUG
                cout << "Command: " << bed_cmd 
                     << "\n\tRunning..." << flush;
                int val = system(bed_cmd.c_str());
                assert(val == 0);
                cout << "FINISH" << endl;
#else
                if(system(bed_cmd.c_str()) != 0){
                    cerr << "problem running " << bed_cmd << endl;
                    exit(1);
                }
#endif
                //I believe this is the correct implementation of the call
                // to bedtools with the binary file
                // located within lib
            }
            else{
                if(system(string("touch " + bedfile).c_str()) != 0){
                    cerr << "problem running touch " << bedfile << endl;
                    exit(1);
                }
            }
        }

    }
}
