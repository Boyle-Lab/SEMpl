extern "C" {
    #include "./lib/libBigWig-master/bigWig.h"
}
#include "iterativeSEM.hpp"
#include "common.hpp"
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

void accumSummary_scale(Dataset &data, const string &hfile,
                        const string &cfile, int scale,
                        Dataset::accumSummaryData::accumSummary_dest dest){

	// open file using library, below code is necessary
	// because of C++ type system regarding const
	char *fname = new char[hfile.length() + 1];
	strcpy(fname, hfile.c_str());
	bigWigFile_t *bwFile = bwOpen(fname, NULL, "r");
    	
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
        temp.clear();
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

		chrom = new char[seqid.length() + 1];
		strcpy(chrom, seqid.c_str());

        cout << "chrom: " << chrom << endl
        << "upstart: " << static_cast<uint32_t>(upstart) << endl
        << "upend: " << static_cast<uint32_t>(upend) << endl
        << "nBins: " << static_cast<uint32_t>(upend - upstart) << endl;

		values = bwStats(bwFile, chrom, static_cast<uint32_t>(upstart),
                         static_cast<uint32_t>(upend),
                         static_cast<uint32_t>(upend - upstart), bwStatsType::mean);


        if(values == NULL){
            cerr << "Failure to use bwStats!\n\tEXITING\n";
            exit(1);
        }

		
		for(counter = 0; counter < upend - upstart; counter++){
			values[counter] = roundf(values[counter] * 1000) / 1000;
			signal_array[counter] = to_string(values[counter]);
		}

        free(values);
		delete [] chrom;


        // UPDATE:
        // nan IS AN ACTUAL POSSIBLE DOUBLE VALUE
        // where nan can be found by nan("") or strtod("nan")
		if(direction.find('+') != string::npos){
			for(int k = 0; k < total_size; ++k){
				// need to determine how to check for definition of signal_array[k]

                // checks that signal_array[k] is non-zero
                if(stod(signal_array[k]) != nan("")){
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
                if(stod(signal_array[k]) != nan("")){
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
                    cerr << "enum accum data should be empty!!! I think!!\n"
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
                    cerr << "scramble accum data should be empty!!! I think!!\n"
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
                    cerr << "align accum data should be empty!!! I think!!\n"
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
