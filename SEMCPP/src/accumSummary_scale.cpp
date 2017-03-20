#include "lib/libBigWig-master/bigWig.h"
#include "iterativeSEM.hpp"
#include "common.h"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
using namespace std;



			// contains all data, contains bigwig filename, region file, scale
//REQUIRES: data is valid Dataset, receives bigwig file, file containing regions to center, and scale size
//MODIFIES: data, specifically accumsummary data
//EFFECTS: gathers raw bigwig data (?)

//Extra parameters should be passed in the struct Dataset

void accumSummary_scale(Dataset &data, const string &hfile,
                        const string &cfile, int scale,
                        Dataset::accumSummaryData::accumSummary_dest dest){

	// open file using library, below code is necessary
	// because of C++ type system regarding const
	char *fname = new char[hfile.length() + 1];
	strcpy(fname, hfile.c_str());
	const char* const mode = "r";
	bigWigFile_t *bwFile = bwOpen(fname, NULL, mode);
    	if(bwFile == NULL){
        	cerr << "Failed to open hfile: " << hfile << '\n';
        	exit(1);
    	}

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
	int start = 0, end = 0, counter = 0;
    int upstart = 0, upend = 0;
	// pointer to hold double values from library function;
	double *values = nullptr;


	while(getline(input, line)){
		// initialize vairables
		start = 0;
		end = 0;
		counter = 0;
		direction = '\0';
		chrom = nullptr;
		split(line, splitBy, temp);
		seqid = temp[0];
        upstart = 0, upend = 0;

		// if lines begins with chr
		if(temp[0][0] == 'c'){
			if(temp[0][1] == 'h')
				if(temp[0][2] == 'r'){
					seqid = temp[0];
#ifdef DEBUG
                    cout << "temp[0] begins with chr: " << temp[0] << '\n';
#endif
                }
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
		//double *bwStats(bigWigFile_t *fp, char *chrom, uint32_t start,
        //                uint32_t end, uint32_t nBins, enum bwStatsType type);
        // should return each individual value, where the mean of a single value
        // is that same value
		values = bwStats(bwFile, chrom, static_cast<uint32_t>(upstart),
                         static_cast<uint32_t>(upend),
                         static_cast<uint32_t>(upend - upstart), bwStatsType::mean);
#ifdef DEBUG
        // check that this is the correect value
        double *test1 = nullptr;
        test1 = bwStats(bwFile, chrom, static_cast<uint32_t>(upstart),
                         static_cast<uint32_t>(upend),
                         static_cast<uint32_t>(upend - upstart), bwStatsType::min);
        double *test2 = nullptr;
        test2 = bwStats(bwFile, chrom, static_cast<uint32_t>(upstart),
                         static_cast<uint32_t>(upend),
                         static_cast<uint32_t>(upend - upstart), bwStatsType::max);
        for(int i = 0; i < upend - upstart; ++i){
            if(test1[i] != test2[i]){
                cerr << "test1 AND test2 should match!!!!!!\nIf I understand it"
                     << " correctly\n";
                cerr << "\ttest1 val: " << test1[i] << "\n\ttest2 val: "
                          << test2[i] << "\n\tEXITING";
                exit(1);
            }
            if(test1[i] != values[i]){
                cerr << "test1 AND values should match!!!!!!\nIf I understand it"
                     << " correctly\n";
                cerr << "\ttest1 val: " << test1[i] << "\n\tvalues val: "
                     << values[i] << "\n\tEXITING";
                exit(1);
            }
        }
        free(test2);
        free(test1);
#endif

        if(values == NULL){
            cerr << "Failure to use bwStats!\n\tEXITING\n";
            exit(1);
        }

		counter = 0;
        //                      I DON'T KNOW IF THIS
        //                      IS THE RIGHT END OF INTERVAL
        //                      UPEND - UPSTART
		for(counter = 0; counter < upend - upstart; counter++){
			values[counter] = roundf(values[counter] * 1000) / 1000;
			signal_array[counter] = to_string(values[counter]);
		}

        free(values);
		delete [] chrom;

		//output results
        // I AM GOING TO ASSUME THAT A VALUE OF 0.0 MEANS UNDEFINED.
        // ACCESSING SOMETHING THAT IS UNDEFINED IN C++ IS A SEGMENTATION FAULT,
        // A RUNTIME ERROR. I DON'T KNOW IF values DOESN'T POINT TO ANY
        // UNDEFINED DATA (0.0 as I am assuming)
		if(direction.find('+') != string::npos){
			for(int k = 0; k < total_size; ++k){
				// need to determine how to check for definition of signal_array[k]

                // checks that signal_array[k] is non-zero
                if(stod(signal_array[k]) != 0.0){
                    output[k] = signal_array[k];
                }
                else{
                    output[k] = "N";
                }
            }
        }
		else{
            // use reserve so I can keep the same iteration direction
            reverse(signal_array.begin(), signal_array.end());
			for(int k = 0; k < total_size; ++k){
                // need to determine how to check for definition of signal_array[k]

                // checks that signal_array[k] is non-zero
                if(stod(signal_array[k]) != 0.0){
                    output[k] = signal_array[k];
                }
                else{
                    output[k] = "N";
                }
            }
        }


		max = 0;
		hitcount = 0;
		for(int l = 0; l < static_cast<int>(output.size()); ++l){
			if(stod(output[l]) > max) max = stod(output[l]);
                                    // string to double
			if(output[l] != string("N") ) ++hitcount;
		}
        // if max is maximum possible double value, then it is not applicable
		if(hitcount / static_cast<double>(output.size()) < 0.9){
			max = numeric_limits<double>::max();
        }
        switch (dest) {
            case Dataset::accumSummaryData::accumSummary_dest::none:
                cerr << "dest shouldn't be none!!!!\n";
                exit(1);
            break;
            case Dataset::accumSummaryData::accumSummary_dest::enumerated:
#ifdef DEBUG
                if(!data.accumSummary_data.enum_accum_lines.empty() 
                    || !data.accumSummary_data.enum_accum_max.empty()){
                    cout << "enum accum data should be empty!!! I think!!\n"
                         << "\tEXITING\n";
                         exit(1);
                }
#endif
                data.accumSummary_data.enum_accum_lines.push_back(line);
                data.accumSummary_data.enum_accum_max.push_back(max);
            break;
            case Dataset::accumSummaryData::accumSummary_dest::scrambled:
#ifdef DEBUG
                if(!data.accumSummary_data.scramble_accum_lines.empty() 
                    || !data.accumSummary_data.scramble_accum_max.empty()){
                    cout << "scramble accum data should be empty!!! I think!!\n"
                         << "\tEXITING\n";
                         exit(1);
                }
#endif
                data.accumSummary_data.scramble_accum_lines.push_back(line);
                data.accumSummary_data.scramble_accum_max.push_back(max);
            break;
            case Dataset::accumSummaryData::accumSummary_dest::alignment:
#ifdef DEBUG
                if(!data.accumSummary_data.align_accum_lines.empty() 
                    || !data.accumSummary_data.align_accum_max.empty()){
                    cout << "align accum data should be empty!!! I think!!\n"
                         << "\tEXITING\n";
                         exit(1);
                }
#endif            
                data.accumSummary_data.align_accum_lines.push_back(line);
                data.accumSummary_data.align_accum_max.push_back(max);
            break;
            default:
                cerr << "there is no default for dest's switch statement!!!\n";
                exit(1);
            break;
        }
	}
	bwClose(bwFile);
	delete [] fname;
}
