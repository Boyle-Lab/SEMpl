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
//REQUIRES: data is valid Dataset, receives bigwig file,
//          file containing regions to center, and scale size
//MODIFIES: data, specifically accumsummary data
//EFFECTS: fills appropriate accumSummary_data vectors

void accumSummary_scale(Dataset &data, const string &hfile,
                        const string &cfile, int scale,
                        Dataset::accumSummary_type::accumSummary_dest dest){

    // clear accum data of corresponding type
    switch (dest) {
            case Dataset::accumSummary_type::accumSummary_dest::enumerated:
                // data.accumSummary_data.enum_accum_max.clear();
                data.accumSummary_data.enum_accum_lines.clear();
            break;
            case Dataset::accumSummary_type::accumSummary_dest::scrambled:
                // data.accumSummary_data.scramble_accum_max.clear();
                data.accumSummary_data.scramble_accum_lines.clear();
            break;
            case Dataset::accumSummary_type::accumSummary_dest::alignment:
                // data.accumSummary_data.align_accum_max.clear();
                data.accumSummary_data.align_accum_lines.clear();
            break;
            case Dataset::accumSummary_type::accumSummary_dest::none:
                cerr << "dest shouldn't be none!!!!" << endl;;
                exit(1);
            break;
            default:
                cerr << "there is no default for dest's switch statement!!!" << endl;
                exit(1);
            break;
    }


    // open file using library, below code is necessary
    // because of C++ type system regarding const
    char *fname = new char[hfile.length() + 1];
    strcpy(fname, hfile.c_str());
    bigWigFile_t *bwFile = bwOpen(fname, NULL, "r");

    if(bwFile == NULL){
    	cerr << "Failed to open hfile: " << hfile << endl;
    	exit(1);
    }

    int dist = 50;
    float max = 0.0;


    //////////////////////////////////
    // Read each peak location and add signal values
    /////////////////////////////////

    ifstream input(cfile);
    if(!input){
        cout << "Failure to open cfile in accumSummary_scale.cpp" << endl;
        exit(1);
    }

    const string splitBy = "\t";
    char *chrom = nullptr;
    string line = "", seqid = "", direction = "";

    vector<string> temp;

    int start = 0, end = 0;
    int upstart = 0, upend = 0;

    #ifdef DEBUG
        //string signal = cfile + ".signal";
        //ofstream sFile(signal);
    #endif

    while( getline(input, line) ){
        // initialize variables
        max = 0.0;
        temp.clear();
        chrom = nullptr;
        split_string(line, splitBy, temp);
        seqid = temp[0];
       	upstart = 0, upend = 0;

        // if lines begins with chr
        if(temp[0][0] == 'c'){
            if(temp[0][1] == 'h'){
                if(temp[0][2] == 'r'){
                    seqid = temp[0];
#ifdef DEBUG
                    // cout << "temp[0] begins with chr: " << temp[0] << '\n';
#endif
                }
            }
        }
        // if line doesn't begin with chr
        else{
            seqid = "chr" + temp[0];
        }
        start = stoi(temp[1]) - 1;

#ifdef DEBUG
        // cerr << "temp[1]: #" << temp[1] << "# stoi: #" << start + 1 << '#' << endl;
#endif
	if(scale != end-start) {
//		cerr << "Incorrect size region!" << endl;
	}

        end = stoi(temp[2]);
        direction = temp[4];
        upstart = start - dist;
        upend = end + dist - 1;

        chrom = new char[seqid.length() + 1];
        strcpy(chrom, seqid.c_str());

#ifdef DEBUG
        //cerr << "Upstart: " << upstart << " Upend: " << upend << " Chr: " << chrom << endl;
#endif

       	bwOverlappingIntervals_t *ptr = bwGetValues(bwFile, chrom, 
                                        static_cast<uint32_t>(upstart),
                                        static_cast<uint32_t>(upend),
                                        1);
        if(!ptr){
            cerr << "problem with bwGetValues!!! " << chrom << " " << upstart << " " << upend << endl 
                 << "\tEXITING" << endl;
            exit(1);
        }

        // cout << "\t for chrom: " << seqid << endl;
        // cout << "\tdeleted chrom" << endl;
        delete [] chrom;

        int hitcount = 0;

        try{
            for(int k = 0; k < (int)(ptr->l); ++k){
                if(!isnan( ptr->value[k] )){
                    ptr->value[k] = roundf(ptr->value[k] * 10000.0) / 10000.0;
                    if(ptr->value[k] > max){
                        max = ptr->value[k];
                    }
                    ++hitcount;
                }
            }
        }
        catch(...){
            cerr << "nan exception thrown" << endl;
            exit(1);
        }


        free(ptr->start);
        free(ptr->end);
        free(ptr->value);
        free(ptr);

//	cerr << "Hitcount: " << hitcount << " Total size: " << total_size << " Max: " << max;

        // Bug in BigWig reading -- 0s seem to not appear
        // Setting this only for all NAs for now
//        if( hitcount < 1){
        if(max == 0) {
            max = NAN_VALUE;
        }

        line += '\t' + to_string(max);

#ifdef DEBUG
        //sFile << line << endl;
#endif

        switch (dest) {
            case Dataset::accumSummary_type::accumSummary_dest::none:
                cerr << "dest shouldn't be none!!!!" << endl;
                exit(1);
            break;
            case Dataset::accumSummary_type::accumSummary_dest::enumerated:
                data.accumSummary_data.enum_accum_lines.push_back(line);
                // data.accumSummary_data.enum_accum_max.push_back(max);
            break;
            case Dataset::accumSummary_type::accumSummary_dest::scrambled:
                data.accumSummary_data.scramble_accum_lines.push_back(line);
                // data.accumSummary_data.scramble_accum_max.push_back(max);
            break;
            case Dataset::accumSummary_type::accumSummary_dest::alignment:
                data.accumSummary_data.align_accum_lines.push_back(line);
                // data.accumSummary_data.align_accum_max.push_back(max);
            break;
            default:
                cerr << "there is no default for dest's switch statement!!!" << endl;
                exit(1);
            break;
        }
	}
        bwClose(bwFile);
        delete [] fname;
}
