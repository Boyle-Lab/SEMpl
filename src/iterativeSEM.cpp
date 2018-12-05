///
//  iterativeSEM.cpp
//  SEMCPP
//
//  Created by Cody Morterud and Colten Williams on 11/9/16.
//  Copyright © 2016 Boyle Lab. All rights reserved.
//

#include "iterativeSEM.hpp"
#include <iostream>
#include <ctime>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <map>
#include <fstream>
#include <getopt.h>
using namespace std;

/*
 example execution from command line
 "./iterativeSEM.pl -PWM examples/MA0114.1.pwm
 -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak
 -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig
 -TF_name HNF4A -output examples/HNF4A/"
*/


int main(int argc, char **argv){

	Dataset data;

	int total_iterations = 250;

	time_t timer;
	time(&timer);

	cout << "Running Iterative SEM building..\n";

    int index = -1;
    const option long_opts[] = {
        {"PWM", required_argument, NULL, 'p'},
        {"merge_file", required_argument, NULL, 'm'},
        {"big_wig", required_argument, NULL, 'b'},
        {"TF_name", required_argument, NULL, 't'},
        {"output", required_argument, NULL, 'o'},
        {"readcache", required_argument, NULL,  'c'},
        {"verbose", no_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };
    char c = '\0';

    while( (c = static_cast<char>(getopt_long_only(argc, argv, "p:m:b:t:o:c:v", long_opts, &index)) ) != -1){
        switch (c) {
            case 'p':
                data.PWM_file = optarg;
#ifdef DEBUG
                cout << "\tPWM: " << optarg << '\n';
#endif
                break;
            case 'm':
                data.DNase_file = optarg;
#ifdef DEBUG
                cout << "\tmerge_file: " << optarg << '\n';
#endif
                break;
            case 'b':
                data.bigwig_file = optarg;
#ifdef DEBUG
                cout << "\tbigwig: " << optarg << '\n';
#endif
                break;
            case 't':
                data.TF_name = optarg;
#ifdef DEBUG
                cout << "\tTF_name: " << optarg << '\n';
#endif
                break;
            case 'o':
                data.base_dir = optarg;
#ifdef DEBUG
                cout << "\t output: " << optarg << '\n';
#endif
                break;
            case 'c':
            if(optarg){
                data.cachefile = optarg;

#ifdef DEBUG
                cout << "\tcachefile flag: " << optarg << '\n';
#endif
            }
            else{
                cout << "readcache, but no argument given, using default\n";
            }
            break;
            case 'v':
#ifdef DEBUG
                cout << "\tverbose\n";
#endif
                data.settings.verbose = true;
                break;
            default:
                cout << "unknown option!" << c << '\n';
                break;
        }
    }

    // must have pwm file
	if(data.PWM_file.empty()){
		cout << "No PWM file given" << endl;
		exit(1);
	}
    // must have output directory
	if(data.base_dir.empty()){
		cout << "No output file given" << endl;
		exit(1);
	}
    // must have a transcription factor name
    if(data.TF_name.empty()){
        cout << "No TF name given" << endl;
    }

    // if no database given, assume default location
    if(data.cachefile.empty())  {
        data.cachefile = data.base_dir + "/CACHE.db";
    }

    vector<double> pvals;
    pvals.reserve(total_iterations + 1);
    pvals.push_back(pow(4, -5));

    for(int iteration =  1; iteration < total_iterations; ++iteration){
        pvals.push_back(pow(4, -5.5));
    }

    //needed to prepare cache
    read_pwm(data, data.PWM_file);

    data.settings.useCache = true;
    connectCache(data, data.cachefile, data.cacheDB);

    data.settings.threads = 20;
    data.settings.maxKmers = 1000000;
    data.settings.minKmers = 32;

    int converge = 0;
    vector<string> line_2;
    int diff = 0;
    int same = 0;
    int total_1 = 0;
    int total_diff = 0;
    string final_run = "";
    string line = "";
    string newPwm = "";
    map <string, double> kmers;

    ofstream outFile(data.base_dir + "/kmer_similarity.out");

    for (int iteration = 0; iteration < total_iterations; iteration++){

        // output the results from comparing kmers
        if(iteration > 1 && converge < 10){
            // compare new and old kmers
            total_1 = same = diff = total_diff = 0;
            // for every key-value pair, find if the kmer is found was found in
            // the other kmer key-value map
            for(const auto &kmer_val_pair : kmers){
                if(data.kmerHash.find(kmer_val_pair.first) != data.kmerHash.end()){
                    // an iterator was found that was not the end iterator
                    // therefore the kmer exists in both
                    ++same;
                }
                else{
                    // an end iterator was returned
                    // kmer doesn't exist
                    ++diff;
                }
                ++total_1;
            }

            if(diff == 0){
                // all kmers are the same
                // converge increases
                ++converge;
            }
            else{
                // reset converge
                converge = 0;
            }
            // assign current kmerHash to kmers
            // so on the next iteration a comparison can be made
            kmers = data.kmerHash;

        }

        if(converge < 10 && iteration < total_iterations){
            outFile << iteration << "\t" << converge << "\t" << same << "\t" << diff << "\n" << flush;
            cout << "---Iteration " << iteration << "---"<< '\n';

            data.output_dir = data.base_dir + "/" + "it" + std::to_string(iteration) + "/";

            if (converge == 9){
                // the folder containing the final iteration data
                final_run = "it" + std::to_string(iteration);
            }

            //get pvalues
            newPwm = data.base_dir + "/" + "it" + std::to_string(iteration-1) + "/" + data.TF_name + ".pwm";
            if(iteration == 0) {
                newPwm = data.PWM_file;
            }
            read_pwm(data, newPwm);

            data.settings.threshold = get_threshold(data, pvals.front());
            if(data.settings.threshold < 0.0){
                data.settings.threshold = 0.0;
            }
            pvals.erase(pvals.begin());

            // generate SEM
            try{
                generateSNPEffectMatrix(data);
            }
            catch(...){
                cerr << "Problem with generateSNPEffectMatrix!!!\n\tEXITING\n";
                exit(1);
            }

            // generate final PWM from final SEM
            try{
                generatePWMfromSEM(data,
                           data.output_dir + "/" + data.TF_name + ".sem",
                           data.output_dir + "/" + data.TF_name + ".pwm");
            }
            catch(...){
                cerr << "Problem with generatePWMfromSEM!!!\n\tEXITING\n";
                exit(1);
            }


            cout << "\n";
        }
        else{
            //stop iterations post-convergence

            break;
        }
    }

    // Everything is done now, so clean up
    closeCache(data.cachefile, data.cacheDB);
    outFile.close();

    //link last iteration
    ostringstream newOutput;
    newOutput << data.base_dir << "/final";
    string cmd = "ln -s " + final_run + " " + newOutput.str();
    if(system(cmd.c_str()) != 0){
        cerr << "problem running " << cmd << endl;
    }
    double diff_t;
    time_t endTime;
    time(&endTime);
    diff_t = difftime(endTime, timer);
    cout << "**************************" << '\n';
    printf( "Job took %f seconds", diff_t );
    cout << "**************************\n" << '\n';

    return 0;
}

// Requires: .pwm file uses '\t' to separate fields
// Effects: fills in PWM_data and returns the second field of the first line
// in the .pwm file specified
string read_pwm(Dataset &data, string file){
    ifstream fin(file);
#ifdef DEBUG
    assert(fin);
#endif
    data.PWM_data.matrix_arr.clear();

    string s = "";
    // ignore first line of pwm
    fin >> s >> s;


    // s now contains second field
    fin.ignore(10000, '\n');
    fin.ignore(10000, '\t');
    int matrix_element = -1;
    for(int column = 0; column < data.PWM_data.NUM_COLUMNS; ++column){
        // adds a vector for each column
        data.PWM_data.matrix_arr.push_back(vector<int>());
    }

    // while fin reads in the first element of a row
    while(fin >> matrix_element){
        #ifdef DEBUG
            // cerr << matrix_element << endl;
        #endif
        // adding a row
        data.PWM_data.matrix_arr[0].push_back(matrix_element);
        for(int index = 1; index < data.PWM_data.NUM_COLUMNS; ++index){
            fin >> matrix_element;
            data.PWM_data.matrix_arr[index].push_back(matrix_element);
        }


        // THE EXAMPLE FILE USES A TAB CHARACTER
        // ignore until tab of next line, which ignores the 1st field
        // of the next line also
        // ignore first character (row number)
#ifdef DEBUG
        assert(data.PWM_data.matrix_arr[0].size() == data.PWM_data.matrix_arr[1].size());
        assert(data.PWM_data.matrix_arr[0].size() == data.PWM_data.matrix_arr[2].size());
        assert(data.PWM_data.matrix_arr[0].size() == data.PWM_data.matrix_arr[3].size());
        // cout << endl;
#endif
        fin.ignore(10000, '\n');
        fin.ignore(10000, '\t');
    } // end while
#ifdef DEBUG
    // verifies to stdout that we are storing the matrix
    cout << "\tPWM\n";
    // matrix_arr's first index is the column
    // matrix_arr's second index is the row
        for(int row = 0; row < (int)data.PWM_data.matrix_arr[0].size() ; ++row){
            cout << row << '\t';
            for(int column = 0; column < (int)data.PWM_data.matrix_arr.size(); ++column){
            // fin >> data.PWM_data.matrix_arr[column][row];


                cout << data.PWM_data.matrix_arr[column][row] << '\t';

            }
            cout << endl;
        }
#endif

    data.settings.length = data.PWM_data.matrix_arr[0].size();

    return s;
}
