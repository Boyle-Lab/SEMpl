#include "iterativeSEM.hpp"
#include "common.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std;


// REQUIRES: Boost library installed
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
                string bed_cmd = "../lib/bedtools intersect -a " + readfile+ " -b " + data.DNase_file + " -wa -u | sort |uniq > " + bedfile;
                system(bed_cmd.c_str());
                //I believe this is the correct implementation of the call to bedtools with the binary file
                // located within lib
            }
            else{
                system(string("touch " + bedfile).c_str());
            }
        }
        // if ($da == 1){
        // rm -rf readfile
        // }
    }


}
