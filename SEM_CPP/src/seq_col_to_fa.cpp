#include "iterativeSEM.hpp"
using namespace std;

int seq_col_to_fa(const vector<string> &column, const string &file){

    // takes input from checkCache, made a mark in checkCache at line 150
    // UPDATE: made changes

    // will need to discuss the "intersect" function, or the bed tools function
    // and adding them onto the program with speed and time in mind

    // writes to a file sucessfully

    if(column.empty()) {
        return 0;
    }

    ofstream OUTF(file);

    int counter = 0;
    for(auto val : column){
        OUTF << '>' << val << '\n' << val << '\n';
        ++counter;
    }

    OUTF.close();

    return counter;
}
