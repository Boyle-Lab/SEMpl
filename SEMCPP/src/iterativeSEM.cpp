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

	int total_iterations = 0;

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
                data.output_dir = optarg;
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

	if(data.PWM_file.empty()){
		cout << "No PWM file given" << endl;
		exit(1);
	}

	if(data.output_dir.empty()){
		cout << "No output file given" << endl;
		exit(1);
	}

    if(data.TF_name.empty()){
        cout << "No TF name given" << endl;
    }

// data.cachefile.empty() checks if the string is empty, not the actual file
	if(data.cachefile.empty())  data.cachefile = data.output_dir + "/CACHE.db";

    vector<double> pvals;
    pvals.reserve(total_iterations + 1);
    pvals.push_back(0.0009765625);

    for(int iteration =  1; iteration <= total_iterations; ++iteration){
        pvals.push_back(0.0004882812);
    }

    double pVal = pvals.front();
    pvals.erase(pvals.begin());

    read_pwm(data, data.PWM_file);
    data.settings.threshold = get_threshold(data, pVal);
    if (data.settings.threshold < 0){
        data.settings.threshold = 0;
    }
#ifdef DEBUG
    cout << "pval: " << pVal << endl;
    cout << "threshold: " << data.settings.threshold << endl;
#endif

    cout << "--- Iteration 0 ---\n";

    data.settings.iteration = 0;
    try{
#ifdef DEBUG
        cout << "\tgenerating SNPEffectMatrix\n" << flush;
#endif
        generateSNPEffectMatrix(data);
    }
    catch(...){
        cerr << "Problem with generateSNPEffectMatrix!!!\n\tEXITING\n";
        exit(1);
    }
    try{
        generatePWMfromSEM(data);
    }
    catch(...){
        cerr << "Problem with generatePWMfromSEM!!!\n\tEXITING\n";
        exit(1);
    }

    pvals.erase(pvals.begin());


    data.settings.threshold = get_threshold(data, pvals.front());
    pvals.erase(pvals.begin());
    if (data.settings.threshold < 0){
        data.settings.threshold = 0;
    }


    int converge = 0;
    vector<string> line_2;
    int diff = 0;
    int same = 0;
    int total_1 = 0;
    int total_diff = 0;
    string final_run = "";
    string line = "";
//    int iterID = 0;

    map <string, double> kmers;
    // kmers_2 is new k-mers or data.kmerHash
    // kmers is old k-mer or previous data.kmerHash
    data.settings.fastrun = false;
    for (int iteration = 1; iteration < total_iterations; ++iteration){
	    
        // iterID = rand() % 16777216;
        ofstream outFile(data.output_dir + "/kmer_similarity.out");
        if(iteration > 1 && converge < 10){
            // int j = iteration - 1;
            total_1 = same = diff = total_diff = 0;
            for(const auto &kmer_val_pair : kmers){
                if(data.kmerHash.find(kmer_val_pair.first) != data.kmerHash.end()){
                    ++same;
                }
                else{
                    ++diff;
                }
                ++total_1;
            }

            if(diff == 0){
                ++converge;
            }
            else{
                converge = 0;
            }
            kmers = data.kmerHash;

        }

        if(converge < 10){
            outFile << iteration << "\t" << converge << "\t" << same << "\t" << diff << "\n";
            cout << "---Iteration " << iteration << "---"<< '\n';
            ostringstream newOutput;
            newOutput << data.output_dir << "/" << "it" << iteration << "/";
            if (converge == 9){
		        data.settings.fastrun = true;
                final_run = newOutput.str();
            }
            generateSNPEffectMatrix(data);
            // kmerHash should be filled in after the above line, within data!!!!

            generatePWMfromSEM(data);

            pVal = pvals.front();
            pvals.erase(pvals.begin());
            // newPwm = data.output_dir+ "/" + tf + ".pwm";
            data.settings.threshold = get_threshold(data, pVal);
            if(data.settings.threshold < 0){
                data.settings.threshold = 0;
            }
            cout << "\n";
        }
        else{
            //stop iterations post-convergence

            //link last iteration
            ostringstream newOutput;
            newOutput << data.output_dir << "/final/";
            string cmd = "ln -s " + final_run + " " + data.output_dir + " " + newOutput.str();
            if(system(cmd.c_str()) != 0){
                cerr << "problem running " << cmd << endl;
                exit(1);
            }
            double diff_t;
            time_t endTime;
            time(&endTime);
            diff_t = difftime(endTime, timer);
            cout << "**************************" << '\n';
            printf( "Job took %f seconds", diff_t );
            cout << "**************************\n" << '\n';
            break;
        }
    }

	return 0;
}

// Requires: .pwm file uses '\t' to separate fields
// Effects: fills in PWM_data and returns the second field of the first line
// in the .pwm file specified
string read_pwm(Dataset &data, string file){
    ifstream fin(file);
    // string s = "";
    // fin.ignore(10000, '\n');
    // int i = 0;
    // while(getline(fin, s)){
    //     if(i == 0 || i == 13){
    //         fin.ignore(10000, '\n');
    //     }
    //     cout << s << endl;
    // }
#ifdef DEBUG
    assert(fin);
#endif


    // fin.ignore(10000, '\n');
    string s = "";
    fin >> s >> s >> s;
    for(int i = 0; i < Dataset::PWM::NUM_ROWS; ++i){
        // THE EXAMPLE FILE USES A TAB CHARACTER
        fin.ignore(10000, '\t');
        for(int j = 0; j < Dataset::PWM::NUM_COLUMNS; ++j){
            fin >> data.PWM_data.matrix_arr[i][j];
#ifdef DEBUG
            // cout << data.PWM_data.matrix_arr[i][j] << ' ';
#endif
        }
#ifdef DEBUG
        // cout << endl;
#endif
        fin.ignore(10000, '\n');
    }
    return s;
}
