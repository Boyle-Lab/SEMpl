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
using namespace std;

/*
 example execution from command line
 "./iterativeSEM.pl -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output examples/HNF4A/"
 */

void generateSNPEffectMatrix(Dataset &data);
void generatePWMfromSEM(Dataset &data);
double get_threshold(Dataset &data);

int main(int argc, char **argv){

	string pwm = "", dnase = "", chip = "", tf = "", output_dir = "", cache = "";

	Dataset data;

	int total_iterations = 500;

	time_t timer;
	time(&timer);

	cout << "Running Iterative SEM building..\n";

	string parse = "";

	for(int i = 0; i < argc; i++){
		cout << argv[i] << ' ';

		data.command += argv[i];
		data.command += " ";

		parse = argv[i];

		if(parse == "-PWM"){
			data.PWM_file = pwm = argv[i+1];
		}
		else if(parse == "-merge_file"){
			data.DNase_file = dnase = argv[i+1];
		}
		else if(parse == "-big_wig"){
			chip = argv[i+1];
		}
		else if(parse == "-TF_name"){
			data.TF_name = tf = argv[i+1];
		}
		else if(parse == "-output"){
			data.output_dir = output_dir = argv[i+1];
		}
		else if(parse == "-readcache"){
			cache = argv[i+1];
		}



	}
	cout << endl;

	if(pwm.empty()){
		cout << "No PWM file given" << endl;
		exit(EXIT_FAILURE);
	}

	if(output_dir.empty()){
		cout << "No output file given" << endl;
		exit(EXIT_FAILURE);
	}

	if(cache.empty())  cache = output_dir + "/CACHE.DB";

    vector<double> pvals;
    pvals.push_back( pow(4, -5));
    double minPval= pow(4,-5.5);
    for(int i = pvals.size(); i <= total_iterations; i++){
        pvals.push_back(minPval);
    }

    pvals.erase(pvals.begin());
    double pVal = pvals.front();

//  ostringstream threshstream;
//  threshstream << "./get_threshold " << pwm << " " << pVal;
//  string threshCmd = threshstream.str();
    data.settings.threshold = get_threshold(data);
    if (data.settings.threshold < 0){				      // as opposed to system(string) calls?
        data.settings.threshold = 0;
    }

    int iterID = rand() % 16777216 ;
    cout << "--- Iteration 0 ---" << endl;

    generateSNPEffectMatrix(data);
    generatePWMfromSEM(data);

/*
*	will change the Cmd's to functions, once the functions are implemented
*/

    pvals.erase(pvals.begin());
    pVal = pvals.front();
    string newPwm = output_dir + "/" + tf + ".pwm";
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
    map <string, int> kmers, kmers_2;
    int diff = 0;
    int same = 0;
    int total_1 = 0;
    int total_diff = 0;
    string final_run = "";
    string line = "";

    for (int i = 1; i < total_iterations; i++){
	data.settings.fastrun = false;
        iterID = rand() % 16777216;
        ofstream outFile(output_dir+"/kmer_similarity.out");
        if(i > 1 && converge < 10){
            int j = i -1;
            total_1 = same = diff = total_diff = 0;
            ostringstream Ekmerstream;
            Ekmerstream << output_dir << "/it" << j << "/Enumerated_kmer.txt";
            string EkmerFile = Ekmerstream.str();
            ifstream Ekmer(EkmerFile);
            if (!Ekmer){
                cerr << "Problem opening " << EkmerFile << endl;
                exit(EXIT_FAILURE);
            }
            while (Ekmer){
                for (string temp; getline(Ekmer, temp, '\n'); line_2.push_back(temp));
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
            newOutput << output_dir << "/" << "it" << i << "/";
            if (converge == 9){
		data.settings.fastrun = true;
            //    wkCmdstream.str("");
            //    wkCmdstream << "./generateSNPEffectMatrix.cpp -PWM " << newPwm << " -merge_file " << dnase << " -big_wig " << chip <<" -TF_name " << tf << " -output " << newOutput.str() << " -threshold " << threshold << " -iteration " << iterID << " -writecache -readcache "<< cache << " -fastrun -verbose";
                final_run = newOutput.str();
            }
            else{
            //    wkCmdstream.str("");
            //    wkCmdstream << "./generateSNPEffectMatrix.cpp -PWM " << newPwm << " -merge_file " << dnase << " -big_wig " << chip <<" -TF_name " << tf << " -output " << newOutput.str() << " -threshold " << threshold << " -iteration " << iterID << " -writecache -readcache "<< cache << " -verbose";
            }
            generateSNPEffectMatrix(data);

            generatePWMfromSEM(data);

            pvals.erase(pvals.begin());
            pVal = pvals.front();
            newPwm = output_dir + "/" + tf + ".pwm";
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
            newOutput << output_dir << "/final/";
            string cmd = "ln -s " + final_run + " " +output_dir + " " + newOutput.str();
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
