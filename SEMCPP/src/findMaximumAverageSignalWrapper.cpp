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

#ifdef DEBUG
    ofstream debug("fMASW.txt");
#endif

    vector<double> values(alignments.size(), 0.0);
#ifdef DEBUG
    debug << "align size " << alignments.size() << endl;
    debug << "values size " << values.size() << endl;
#endif
    
    
        for(size_t i = 0; i < alignments.size(); ++i){
            //Finds maximums of each line and stores into a vector called
            // maximums.
    #ifdef DEBUG
            debug << "line: #" << line << '#' << " i: " << i << endl;
    #endif
            line = alignments[i];

            // removed from last try block
            val = grab_string_last_index(line);
    #ifdef DEBUG
            debug << "val:#" << val << '#' << endl;
    #endif
            try{
                values.at(i) = stod(val);
            }
            catch(invalid_argument e){
                cerr << "hee" << e.what() << endl;
                cerr << "val: #" << val << '#' << endl;
                exit(1);
            }
            catch(out_of_range e){
                cerr << "from .at: " << !(i < values.size() ) << endl;
                cerr << e.what() << endl;
                exit(1);
            }
            catch(...){
                cerr << "problem iterating findMaximumAverageSignalWrapper" << endl;
                exit(1);
            }
    #ifdef DEBUG
            debug << "\tval: #" << values.at(i) << '#' << endl;
    #endif 
    
    #ifdef DEBUG 
            debug << "here1" << endl;
    #endif

            if(values.at(i) != NAN_VALUE){
                #ifdef DEBUG
                    debug << "here2 val: #" << values.at(i) << '#' << endl
                          << '#' << sum << '#' << endl;
                #endif
                sum += values.at(i);
                #ifdef DEBUG
                    debug << "here3" << endl;
                #endif
                ++counter_out;
                #ifdef DEBUG
                    debug << "here4: #" << counter_out << '#' << endl;
                #endif
            }
        }
    
    
            // cout << "value at " << i << ": " << max_ptr->at(i) << endl;
        

    if(counter_out > 0){
        mean_out = sum / counter_out;
    }

    double sqtotal = 0.0;
    for(size_t i = 0; i < values.size(); ++i){
        sqtotal += pow( (mean_out - values[i] ), 2.0);
    }
    stdev_out = pow( (sqtotal / (values.size() - 1) ), 0.5);
    sterr_out = stdev_out / pow(counter_out, 0.5);




        // print OUT_HANDLE "$file\t$maximum\t$count\t$stdev\t$sterr\n";

        
    

}
