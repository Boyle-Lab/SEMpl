//
//  iterativeSEM.cpp
//  SEMCPP
//
//  Created by Cody Morterud on 11/9/16.
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
using namespace std;

/*
 example execution from command line
 "./iterativeSEM.pl -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output examples/HNF4A/"
 */


int main(int argc, char **argv){

	string pwm = "", dnase = "", chip = "", tf = "", output = "", cache = "";

    int total_iterations = 500;

	time_t timer;
	time(&timer);

	cout << "Running Iterative SEM building..\n";

	string parse = "";

	for(int i = 0; i < argc; i++){
		cout << argv[i] << ' ';

		parse = argv[i];

		if(parse == "-PWM"){
			pwm = argv[i+1];
		}
		else if(parse == "-merge_file"){
			dnase = argv[i+1];
		}
		else if(parse == "-big_wig"){
			chip = argv[i+1];
		}
		else if(parse == "-TF_name"){
			tf = argv[i+1];
		}
		else if(parse == "-output"){
			output = argv[i+1];
		}
		else if(parse == "-readcache"){
			cache = argv[i+1];
		}

	}
	cout << endl;

	if(cache.empty())  cache = output + "/CACHE.DB";

    vector<double> pvals;
    pvals.push_back( pow(4, -5));
    double minPval= pow(4,-5.5);
    for(int i = pvals.size(); i <= total_iterations; i++){
        pvals.push_back(minPval);
    }

    pvals.erase(pvals.begin());
    double pVal = pvals.front();
    ostringstream threshstream;
    threshstream << "./get_threshold " << pwm << " " << pVal;
    string threshCmd = threshstream.str();
    double threshold = system(threshCmd.c_str());     // I believe we are supposed make function calls here? 
    if (threshold < 0){				      // as opposed to system(string) calls?
        threshold = 0;
    }

    int iterID = rand() % 16777216 ;
    cout << "--- Iteration 0 ---" << endl;
    ostringstream wkCmdstream;
    wkCmdstream << "./generateSNPEffectMatrix.cpp -PWM " << pwm << " -TF_name " << tf << " -output " << output;
    string wkCmd = wkCmdstream.str();
    system(wkCmd.c_str());
    ostringstream pwmCmdstream;
    pwmCmdstream << "./src/generatePWMfromSEM.cpp -PWM " << pwm << " -TF_name " << tf << " -output " << output;
    string pwmCmd = pwmCmdstream.str();
    system(pwmCmd.c_str());

/*
*	will change the Cmd's to functions, once the functions are implemented
*/

    pvals.erase(pvals.begin());
    pVal = pvals.front();
    string newPwm = output + "/" + tf + ".pwm";
    threshstream.str("");
    threshstream << "src/get_threshold.cpp " << newPwm << " " << pVal;
    threshCmd = threshstream.str();
    threshold = system(threshCmd.c_str());
    if (threshold < 0){
        threshold = 0;
    }

    int converge = 0;
    vector<string> line_2;
    map <string, int> kmers;
    map <string, int> kmers_2;
    int diff = 0;
    int same = 0;
    int total_1 = 0;
    int total_diff = 0;
    string final_run = "";
    string line = "";

    for (int i = 1; i < total_iterations; i++){
        iterID = rand() % 16777216 ;
        ofstream outFile(output+"/kmer_similarity.out");
        if(i > 1 && converge < 10){
            int j = i -1;
            total_1 = same = diff = total_diff = 0;
            ostringstream Ekmerstream;
            Ekmerstream << output << "/it" << j << "/Enumerated_kmer.txt";
            string EkmerFile = Ekmerstream.str();
            ifstream Ekmer(EkmerFile);
            if (!Ekmer){
                cerr << "Problem opening " << EkmerFile << endl;
                exit(1);
            }
            while (Ekmer){
                for (string temp; getline(Ekmer, temp, ' '); line_2.push_back(temp));
                kmers_2[line_2[1]] = 1;
                total_1++;
                if (kmers_2.find(line_2[1]) != kmers_2.end()){
                    same++;
                }
                else{
                    diff++;
                }
            }
            if (diff == 0){
                converge++;
            }
            else{
                converge = 0;
            }
            kmers.clear();
            kmers.insert(kmers_2.begin(),kmers_2.end());
            kmers_2.clear();
        }

        if(converge < 10){
            outFile << i << "\t" << converge << "\t" << same << "\t" << diff << "\n";
            cout << "---Iteration " << i << "---"<< endl;
            ostringstream newOutput;
            newOutput << output << "/" << "it" << i << "/";
            if (converge == 9){
                wkCmdstream.str("");
                wkCmdstream << "./generateSNPEffectMatrix.cpp -PWM " << newPwm << " -merge_file " << dnase << " -big_wig " << chip <<" -TF_name " << tf << " -output " << newOutput.str() << " -threshold " << threshold << " -iteration " << iterID << " -writecache -readcache "<< cache << " -fastrun -verbose";
                final_run = newOutput.str();
            }
            else{
                wkCmdstream.str("");
                wkCmdstream << "./generateSNPEffectMatrix.cpp -PWM " << newPwm << " -merge_file " << dnase << " -big_wig " << chip <<" -TF_name " << tf << " -output " << newOutput.str() << " -threshold " << threshold << " -iteration " << iterID << " -writecache -readcache "<< cache << " -verbose";
            }
            wkCmd = wkCmdstream.str();
            system(wkCmd.c_str());

            pwmCmdstream.str("");
            pwmCmdstream << "./src/generatePWMfromSEM.cpp -PWM " << newPwm << " -TF_name " << tf << " -output " << newOutput.str();
            pwmCmd = pwmCmdstream.str();
            system(pwmCmd.c_str());

            pvals.erase(pvals.begin());
            pVal = pvals.front();
            newPwm = output + "/" + tf + ".pwm";
            threshstream.str("");
            threshstream << "src/get_threshold.cpp" << newPwm << pVal;
            threshCmd = threshstream.str();
            threshold = system(threshCmd.c_str());
            if(threshold < 0){
                threshold = 0;
            }
            cout << "\n";
        }
        else{
            //stop iterations post-convergence

            //link last iteration
            ostringstream newOutput;
            newOutput << output << "/final/";
            string cmd = "ln -s " + final_run + " " +output + " " + newOutput.str();
            system(cmd.c_str());
            double diff_t;
            time_t endTime;
            time(&endTime);
            diff_t = difftime(endTime, timer);
            cout << "**************************" << endl;
            printf( "Job took %f seconds", diff_t );
            cout << "**************************\n" << endl;
            break;
        }
    }



	return 0;
}
