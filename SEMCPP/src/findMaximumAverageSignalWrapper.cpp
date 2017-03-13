#include "iterativeSEM.hpp"
#include "common.h"
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

    /* Effects: Reads in the data from filterDnaseWrapper
    *           will output integer values into the data struct
    *            containing the maximum, count, stdev, and sterr
    */

//vector<string> split(const string &s, char delim);
                                                    // PASS scramble FOR BASELINE
void findMaximumAverageSignalWrapper(Dataset &data,
                                     Dataset::accumSummaryData::accumSummary_dest dest){

        vector<string> *max_ptr = nullptr;

        switch (dest) {
            case Dataset::accumSummaryData::accumSummary_dest::none:
                cerr << "dest shouldn't be none!!!!\n";
                exit(1);
            break;
            case Dataset::accumSummaryData::accumSummary_dest::enumerated:
                max_ptr = &data.accumSummary_data.enum_accum_max;
            break;
            case Dataset::accumSummaryData::accumSummary_dest::scrambled:
                max_ptr = &data.accumSummary_data.scramble_accum_max;
            break;
            case Dataset::accumSummaryData::accumSummary_dest::alignment:
                max_ptr = &data.accumSummary_data.align_accum_max;
            break;
            default:
                cerr << "there is no default for dest's switch statement!!!\n";
                exit(1);
            break;
        }

        if(max_ptr->empty()){
            cerr << "corresponding accumSummary max data is missing!!!!\n";
            exit(1);
        }

        double sum = 0.0;
        int counter = 0;
        double mean = 0.0;
        double stdev = 0.0;
        double sterr = 0.0;

        //Finds maximums of each line and stores into a vector called
        // maximums.

        for(int i = 0; i < max_ptr->size(); ++i){
            // checks if the value corresponds to "NA"
            if(max_ptr->at(i) != numeric_limits<double>::max()){
                sum += max_ptr->at(i);
                ++counter;
            }
        }

        if(counter > 0){
            mean = sum / counter;
        }

        // The following is the calculations for the standard deviation
        // and standard error.
        double sqtotal = 0.0;
        for (int i = 0; i < data.max_ptr->size(); i++){
            sqtotal += pow((mean - max_ptr.at(i), 2.0);
        }
        stdev = pow( (sqtotal / (max_ptr->size()-1) ), 0.5);
        sterr = stdev / pow(counter, 0.5);

        switch (dest) {
            case Dataset::accumSummaryData::accumSummary_dest::none:
                cerr << "dest shouldn't be none!!!!\n";
                exit(1);
            break;
            case Dataset::accumSummaryData::accumSummary_dest::enumerated:
            case Dataset::accumSummaryData::accumSummary_dest::scrambled:
                data.signal_Data.base_counter = counter;
                data.signal_Data.base_maximum = mean;
                data.signal_Data.base_stdev = stdev;
                data.signal_Data.base_sterr = sterr;
            break;
            case Dataset::accumSummaryData::accumSummary_dest::alignment:
                data.signal_Data.alignment_counter = counter;
                data.signal_Data.alignment_maximum = mean;
                data.signal_Data.alignment_stdev = stdev;
                data.signal_Data.alignment_sterr = sterr;
            break;
            default:
                cerr << "there is no default for dest's switch statement!!!\n";
                exit(1);
            break;
        }


}

//Function used for parsing lines by a delimiter.
/*
vector<string> split(const string &str, char delim) {
    stringstream stream(str);
    string temp;
    vector<string> contents;
    while (getline(stream, temp, delim)) {
        contents.push_back(temp);
    }
    return contents;
}
*/
