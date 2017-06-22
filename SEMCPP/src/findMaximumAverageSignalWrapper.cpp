#include "iterativeSEM.hpp"
#include "common.hpp"
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

// REQUIRES: appropriate accumSummary_data is filled in
//           corresponding to type
// EFFECTS: Reads in the data from filterDnaseWrapper
//          will output integer values into the data struct
//          containing the maximum, count, stdev, and sterr

// NOTE: PASS scramble FOR BASELINE
void findMaximumAverageSignalWrapper(Dataset &data,
                                     Dataset::accumSummary_type::accumSummary_dest dest){

        vector<double> *max_ptr = nullptr;
        //vector<string> *line_ptr = nullptr;

        switch (dest) {
            case Dataset::accumSummary_type::accumSummary_dest::none:
                cerr << "dest shouldn't be none!!!!\n";
                exit(1);
            break;
            case Dataset::accumSummary_type::accumSummary_dest::enumerated:
                max_ptr = &data.accumSummary_data.enum_accum_max;
                //line_ptr = &data.accumSummary_data.enum_accum_lines;
            break;
            case Dataset::accumSummary_type::accumSummary_dest::scrambled:
                max_ptr = &data.accumSummary_data.scramble_accum_max;
                //line_ptr = &data.accumSummary_data.scramble_accum_lines;
            break;
            case Dataset::accumSummary_type::accumSummary_dest::alignment:
                max_ptr = &data.accumSummary_data.align_accum_max;
                //line_ptr = &data.accumSummary_data.align_accum_lines;
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

        for(size_t i = 0; i < max_ptr->size(); ++i){
            // checks if the value corresponds to "NA"
            // if(max_ptr->at(i) != numeric_limits<double>::max()){
            // if( !isnan( max_ptr->at(i) ) ) {
            if( max_ptr->at(i) != NAN_VALUE ) {
                sum += max_ptr->at(i);
                ++counter;
            }

            // cout << "value at " << i << ": " << max_ptr->at(i) << endl;
        }

        if(counter > 0){
            mean = sum / counter;
        }
        // cout << "max_ptr->size(): " << max_ptr->size() << endl;
        // cout << "mean: " << mean << endl << "sum: " << sum << endl 
        //      << "counter: " << counter << endl;

        // The following is the calculations for the standard deviation
        // and standard error.
        double sqtotal = 0.0;
        for (size_t i = 0; i < max_ptr->size(); i++){
            sqtotal += pow((mean - max_ptr->at(i)), 2.0);
        }
        stdev = pow( (sqtotal / ( max_ptr->size() - 1) ), 0.5);
        sterr = stdev / pow(counter, 0.5);

        // print OUT_HANDLE "$file\t$maximum\t$count\t$stdev\t$sterr\n";

        switch (dest) {
            case Dataset::accumSummary_type::accumSummary_dest::enumerated:
                data.Signal_data.enumerate_maximum = mean;
                data.Signal_data.enumerate_counter = counter;
                data.Signal_data.enumerate_stdev = stdev;
                data.Signal_data.enumerate_sterr = sterr;
            break;
            case Dataset::accumSummary_type::accumSummary_dest::scrambled:
                data.Signal_data.scramble_maximum = mean;
                data.Signal_data.scramble_counter = counter;
                data.Signal_data.scramble_stdev = stdev;
                data.Signal_data.scramble_sterr = sterr;
            break;
            case Dataset::accumSummary_type::accumSummary_dest::alignment:
                data.Signal_data.alignment_maximum = mean;
                data.Signal_data.alignment_counter = counter;
                data.Signal_data.alignment_stdev = stdev;
                data.Signal_data.alignment_sterr = sterr;
            break;
            case Dataset::accumSummary_type::accumSummary_dest::none:
                cerr << "dest shouldn't be none!!!!\n";
                exit(1);
            break;
            default:
                cerr << "there is no default for dest's switch statement!!!\n";
                exit(1);
            break;
        }


}
