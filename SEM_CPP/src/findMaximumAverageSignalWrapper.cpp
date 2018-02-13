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
//              MAYBE
void findMaximumAverageSignalWrapper(const std::vector<std::string> &alignments,
                                     double &mean_out, int &counter_out,
                                     double &stdev_out, double &sterr_out){

    mean_out = 0.0;
    counter_out = 0;
    stdev_out = 0.0;
    sterr_out = 0.0;

    if(alignments.empty()){
        cerr << "corresponding findMaximumAverageSignalWrapper"
             << " data is missing!!!!"
             << endl;
        exit(1);
    }


    double sum = 0.0;
    string line = "", val = "";
    vector<double> values;
    double max = 0.0;

    for(size_t i = 0; i < alignments.size(); ++i){
        //Finds maximums of each line and stores into a vector called
        // maximums.
        line = alignments[i];

        val = grab_string_last_index(line);

        max = stod(val);

        if(max != NAN_VALUE){
            values.push_back(max);
            sum += max;
            ++counter_out;
        }
    }

    if(counter_out > 0){
        mean_out = sum / counter_out;
    }

    double sqtotal = 0.0;
    for(size_t i = 0; i < values.size(); ++i){
        sqtotal += pow( (mean_out - values[i] ), 2.0);
    }

    stdev_out = pow( (sqtotal / (values.size() - 1) ), 0.5);
    sterr_out = stdev_out / pow(counter_out, 0.5);

    cout << "Mean: " << mean_out << " Total: " << alignments.size() << " Count: " << counter_out << " Stdev: " << stdev_out << " Sterr: " << sterr_out << endl;

}
