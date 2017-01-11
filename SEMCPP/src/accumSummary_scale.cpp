#include "../lib/libBigWig-master/BigWig.h"
#include "iterativeSEM.hpp"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
using namespace std;


// taken from internet
void split(string str, string splitBy, vector<string>& tokens);

			// contains all data, contains bigwig filename, region file, scale
//REQUIRES: data is valid Dataset, receives bigwig file, file containing regions to center, and scale size
//MODIFIES: data, specifically accumsummary data
//EFFECTS: gathers raw bigwig data (?)

//Extra parameters should be passed in the struct Dataset

void accumSummary_scale(Dataset &data, string hfile, string cfile, int scale){

	// open file using library, below code is necessary
	// because of C++ type system regarding const
	char *fname = new char[hfile.length() + 1];
	strcpy(fname, hfile.c_str());
	const char* const mode = "r";
	bigWigFile_t *bwFile = bwOpen(fname, NULL, mode);

	int dist = 500;
	int total_size = dist * 2 + scale;
	double max = 0.0;
	int hitcount = 0;

	// create a vector with correct size, reduce memory use as opposed to repeatedly creating a new vector in the loop below

	vector<string> output(total_size);

	//////////////////////////////////
	// Read each peak location and add signal values
	/////////////////////////////////

	ifstream input(cfile);
	if(!input){
		cout << "Failure to open cfile in accumSummary_scale.cpp\n";
		exit(1);
	}
	const string splitBy = "\t";
	char *chrom = nullptr;
	string line = "", seqid = "", direction = "";
	vector<string> temp, signal_array(total_size);
	int start = 0, end = 0, counter = 0, upstart = 0, upend = 0;
	// pointer to hold double values from library function;
	double *values = nullptr;
	bwStatsType type = mean;

	while(getline(input, line)){
		// initialize vairables
		start = 0;
		end = 0;
		counter = 0;
		direction = '\0';
		chrom = nullptr;
		split(line, splitBy, temp);
		seqid = temp[0];

		// if lines begins with chr
		if(temp[0][0] == 'c'){
			if(temp[0][1] == 'h')
				if(temp[0][2] == 'r')
					seqid = temp[0];
		}
		// if line doesn't begin with chr
		else{
			seqid = "chr" + temp[0];
		}
		start = stoi(temp[1]) - 1;
		end = stoi(temp[2]);
		direction = temp[4];

		upstart = start - dist;
		upend = end + dist;

		//char *fname = new char[hfile.length() + 1];
		//strcpy(fname, hfile.c_str());
		chrom = new char[seqid.length() + 1];
		strcpy(chrom, seqid.c_str());

		// stats function used from library
		//double *bwStats(bigWigFile_t *fp, char *chrom, uint32_t start, uint32_t end, uint32_t nBins, enum bwStatsType type);
		values = bwStats(bwFile, chrom, static_cast<uint32_t>(0), static_cast<uint32_t>(total_size), static_cast<uint32_t>(total_size), type);
    if(values == NULL){
      cerr << "Failure to use bwStats!\n\tEXITING\n";
      exit(1);
    }
		counter = 0;

		for(counter = 0; counter < total_size; counter++){
			values[counter] = roundf(values[counter] * 1000) / 1000;
			signal_array[counter] = values[counter];
		}

    free(values);
		delete [] chrom;

		//output results
		if(direction.find('+') != string::npos)
			for(int k = 0; k < total_size; k++)
				// need to determine how to check for definition of signal_array[k]
				output[k] = signal_array[k];
		else
			for(int k = total_size - 1; k >= 0; k--)
				output[k] = signal_array[k];


		max = 0;
		hitcount = 0;
		for(int l = 0; l < static_cast<int>(output.size()); l++){
			if(stod(output[l]) > max)
				max = stod(output[l]);
			if(output[l] != "N")
				hitcount++;
		}

		if(hitcount / static_cast<double>(output.size()) < 0.9)
			max = numeric_limits<double>::max();
		// if max is maximum possible double value, then it is not applicable

		data.accumSummary_data.accum_lines.push_back(line);
		data.accumSummary_data.accum_max.push_back(max);
		}
	bwClose(bwFile);
	delete [] fname;
}
