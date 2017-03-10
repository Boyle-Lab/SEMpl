#include "iterativeSEM.hpp"
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

vector<string> split(const string &s, char delim);

void findMaximumAverageSignalWrapper(Dataset &data, const string &file_dir){

    if (data.accumSummary_data.accum_max.size() != 0){
            //Assures that there is data to be used
//            vector<double> maximums;
            double sum = 0.0;
            int counter = 0;
            double mean = 0.0;
            double stdev = 0.0;
            double sterr = 0.0;

        //Finds maximums of each line and stores into a vector called
        // maximums.

        for(int i =0; i < data.accumSummary_data.accum_max.size(); i++){
            sum += data.accumSummary_data.accum_max[i];
            counter += 1;
        }

        mean = sum/counter;

        //The following is the calculations for the standard deviation
        //and standard error.
        double sqtotal = 0.0;
        for (int i = 0; i < data.accumSummary_data.accum_max.size(); i++){
            sqtotal += pow((mean - data.accumSummary_data.accum_max[i]), 2.0);
        }
        stdev = pow((sqtotal/(data.accumSummary_data.accum_max.size()-1)), 0.5);
        sterr = stdev/pow(counter, 0.5);

        //Stores the data according to where it would previously be
        //held in file directory
        if (file_dir == "baseline"){
            data.signal_Data.base_counter = counter;
            data.signal_Data.base_maximum = mean;
            data.signal_Data.base_stdev = stdev;
            data.signal_Data.base_sterr = sterr;
        }
        else if (file_dir == "alignment"){
            data.signal_Data.alignment_counter = counter;
            data.signal_Data.alignment_maximum = mean;
            data.signal_Data.alignment_stdev = stdev;
            data.signal_Data.alignment_sterr = sterr;
        }
    }

}

//Function used for parsing lines by a delimiter.
vector<string> split(const string &str, char delim) {
    stringstream stream(str);
    string temp;
    vector<string> contents;
    while (getline(stream, temp, delim)) {
        contents.push_back(temp);
    }
    return contents;
}
