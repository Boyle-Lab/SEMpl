//
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

void generateSNPEffectMatrix(Dataset &data);

int main(int argc, char **argv){

	string pwm = "", dnase = "", chip = "", tf = "";

	Dataset data;

	int total_iterations = 500;

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
        {"readcache", optional_argument, NULL,  'c'},
        {0, 0, 0, 0}
    };
    char c = '\0';
    c = static_cast<char>(getopt_long_only(argc, argv, "p:m:b:t:o:c::", long_opts, &index));
    do {
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
            if(optarg)
                data.cachefile = optarg;
#ifdef DEBUG
                cout << "\tcachefile: " << optarg << '\n';
#endif
                break;
            default:
                cout << "unknown option!" << c << '\n';
                break;
        }
        c = static_cast<char>(getopt_long_only(argc, argv, "p:m:b:t:o:c::", long_opts, &index));
    } while(c != -1);

	if(pwm.empty()){
		cout << "No PWM file given" << '\n';
		exit(1);
	}

	if(data.output_dir.empty()){
		cout << "No output file given" << '\n';
		exit(1);
	}

// data.cachefile.empty() checks if the string is empty, not the actual file
	if(data.cachefile.empty())  data.cachefile = data.output_dir + "/CACHE.DB";

    vector<double> pvals(total_iterations + 1);
    pvals.push_back(0.0009765625);
    for(int iteration = pvals.size(); iteration <= total_iterations; ++iteration){
        pvals.push_back(0.0004882812);
    }

    pvals.erase(pvals.begin());
    double pVal = pvals.front();

//  ostringstream threshstream;
//  threshstream << "./get_threshold " << pwm << " " << pVal;
//  string threshCmd = threshstream.str();
    data.settings.threshold = get_threshold(data);
    if (data.settings.threshold < 0){
        data.settings.threshold = 0;
    }

    int iterID = 0;
    cout << "--- Iteration 0 ---" << '\n';

    data.settings.iteration = 0;
    generateSNPEffectMatrix(data);
    generatePWMfromSEM(data);

/*
*	will change the Cmd's to functions, once the functions are implemented
*/

    pvals.erase(pvals.begin());
    pVal = pvals.front();
//    string newPwm = data.output_dir + "/" + tf + ".pwm";
//    threshstream.str("");
//    threshstream << "src/get_threshold.cpp " << newPwm << " " << pVal;
//    threshCmd = threshstream.str();
//    threshold = threshCmd.c_str();
//    threshold = get_threshold( something here ) ;

    data.settings.threshold = get_threshold(data);
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

    map <string, int> kmers, kmers_2;

    for (int iteration = 1; iteration < total_iterations; iteration++){
	      data.settings.fastrun = false;
        iterID = rand() % 16777216;
        ofstream outFile(data.output_dir + "/kmer_similarity.out");
        if(iteration > 1 && converge < 10){
            int j = iteration - 1;
            total_1 = same = diff = total_diff = 0;
            ostringstream Ekmerstream;
            Ekmerstream << data.output_dir<< "/it" << j << "/Enumerated_kmer.txt";
            string EkmerFile = Ekmerstream.str();
            ifstream Ekmer(EkmerFile);
            if (!Ekmer){
                cerr << "Problem opening " << EkmerFile << '\n';
                exit(1);
            }
            while (Ekmer){
                for (string temp; getline(Ekmer, temp, '\n'); line_2.push_back(temp));
                kmers_2[line_2[1]] = 1;
                ++total_1;
                if (kmers_2.find(line_2[1]) != kmers_2.end()){
                    ++same;
                }
                else{
                    ++diff;
                }
            }

            if (diff == 0){
                ++converge;
            }
            else{
                converge = 0;
            }
            kmers.clear();
            kmers.insert(kmers_2.begin(), kmers_2.end());
            kmers_2.clear();
        }

        if(converge < 10){
            outFile << iteration << "\t" << converge << "\t" << same << "\t" << diff << "\n";
            cout << "---Iteration " << iteration << "---"<< '\n';
            ostringstream newOutput;
            newOutput <<data.output_dir<< "/" << "it" << iteration << "/";
            if (converge == 9){
		            data.settings.fastrun = true;
            //    wkCmdstream.str("");
            //    wkCmdstream << "./generateSNPEffectMatrix.cpp -PWM " << newPwm << " -merge_file " << dnase << " -big_wig " << chip <<" -TF_name " << tf << " -output " << newOutput.str() << " -threshold " << threshold << " -iteration " << iterID << " -writedata.cachefile-readdata.cachefile"<< data.cachefile<< " -fastrun -verbose";
                final_run = newOutput.str();
            }
            else{
            //    wkCmdstream.str("");
            //    wkCmdstream << "./generateSNPEffectMatrix.cpp -PWM " << newPwm << " -merge_file " << dnase << " -big_wig " << chip <<" -TF_name " << tf << " -output " << newOutput.str() << " -threshold " << threshold << " -iteration " << iterID << " -writedata.cachefile-readdata.cachefile"<< data.cachefile<< " -verbose";
            }
            generateSNPEffectMatrix(data);
            // kmerHash should be filled in after the above line, within data!!!!

            generatePWMfromSEM(data);

            pvals.erase(pvals.begin());
            pVal = pvals.front();
            // newPwm = data.output_dir+ "/" + tf + ".pwm";
            data.settings.threshold = get_threshold(data);
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
            system(cmd.c_str());
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
