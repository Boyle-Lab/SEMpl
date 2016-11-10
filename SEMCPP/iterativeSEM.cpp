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
using namespace std;

/*
 example execution from command line
 "./iterativeSEM.pl -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output examples/HNF4A/"
 */


int main(int argc, char **argv){

	string pwm = "", dnase = "", chip = "", tf = "", output = "", cache = "";

	time_t timer;
	time(&timer);

	cout << "Running Iterative SEM building..\n";

	string parse;

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
	


	return 0;
}
