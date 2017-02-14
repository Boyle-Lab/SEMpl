#include "iterativeSEM.hpp"
#include "../lib/sqlite3/sqlite3.h"
#include "common.h"
using namespace std;

// REQUIRES: accumSummary_scale is filled with the correct data
// EFFECTS: writes output of accumSummary_scale to cache
void writeCache(const Dataset &data, const string &cache){
    /*
    *   takes the output from accumSummary_scale
    *
    */
    // take infile, grab 3rd field, indexed from 0, and then grab entire line
    // store in %kmers, which is a hash

    if(data.settings.verbose){
        cout << "Building cache for processed kmers.\n";
    }

    


    if(!fileExists(cache)){

    }
    else{

    }
}
