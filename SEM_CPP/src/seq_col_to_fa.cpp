#include "iterativeSEM.hpp"
using namespace std;

int seq_col_to_fa(const vector<string> &column, const string &file){

    if(column.empty()) {
        return 0;
    }

    ofstream OUTF(file);

    int counter = 0;
    for(auto val : column){
//        OUTF << '>' << val << '\n' << val << '\n';
        OUTF << val << '\n';
        ++counter;
    }

    OUTF.close();

    return counter;
}
