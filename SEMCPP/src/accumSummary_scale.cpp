#include "../lib/libBigWig-master/BigWig.h"
#include "../iterativeSEM.hpp"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// taken from internet
void split(string str, string splitBy, vector<string>& tokens);

			// contains all data, contains bigwig filename, region file, scale
void accumySummary_scale(Dataset &data, string hfile, string cfile, int scale){

	// open file using library, below code is necessary
	// because of C++ type system regarding const
	char *fname = new char[hfile.length() + 1];
	strcpy(fname, hfile.c_str());
	const char* const mode = "r";
	bigWigFile_t *bwFile = bwOpen(fname, NULL, mode);

	int dist = 500;
	int total_size = dist * 2 + scale;


	//////////////////////////////////
	// Read each peak location and add signal values
	/////////////////////////////////

	ifstream input(cfile);
	if(!input){
		cout << "Failure to open cfile in accumSummary_scale.cpp\n";
		exit(EXIT_FAILURE);
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
		
		for(counter = 0; counter < total_size; counter++){
			values[counter] = roundf(values[counter] * 1000) / 1000;
			signal_array[counter] = values[counter];
		}
		
		delete [] chrom;
		
	}




	//void bwClose(bigWigFile_t *fp);

	bwClose(bwFile);
	delete [] fname;
}

// function taken from internet that splits  a string by a character
void split(std::string str, std::string splitBy, std::vector<std::string>& tokens)
{
    /* Store the original string in the array, so we can loop the rest
     * of the algorithm. */
    tokens.push_back(str);

    // Store the split index in a 'size_t' (unsigned integer) type.
    size_t splitAt;
    // Store the size of what we're splicing out.
    size_t splitLen = splitBy.size();
    // Create a string for temporarily storing the fragment we're processing.
    std::string frag;
    // Loop infinitely - break is internal.
    while(true)
    {
        /* Store the last string in the vector, which is the only logical
         * candidate for processing. */
        frag = tokens.back();
        /* The index where the split is. */
        splitAt = frag.find(splitBy);
        // If we didn't find a new split point...
        if(splitAt == string::npos)
        {
            // Break the loop and (implicitly) return.
            break;
        }
        /* Put everything from the left side of the split where the string
         * being processed used to be. */
        tokens.back() = frag.substr(0, splitAt);
        /* Push everything from the right side of the split to the next empty
         * index in the vector. */
        tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
    }
}
