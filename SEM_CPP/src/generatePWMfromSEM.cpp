
#include "iterativeSEM.hpp"
#include <iostream>
#include <ctime>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <map>
#include <fstream>
#include <string>
#include "src/common.hpp"

using namespace std;

// parses pwm, places appropriate data in map
void static parse_pwm(const string &pwm, map<char, vector<double> > &map);

//REQUIRES: requires data is a valid Dataset with pwm filled in
//MODIFIES: modifies data, specifically pwm
//EFFECTS: Creates a pwm matrix from the previous sem matrix in same
//         directory
void generatePWMfromSEM(const Dataset & data, 
                        string input, // input SEM file
                        string output){ // output PWM file

    // raw baseline is the output from running findMaximumAverageSignal... on
    // "enumerate" data

    // sem file for TF
    // string rawInput = data.output_dir + "/" + data.TF_name + ".sem";
    // pwm file for TF
    // string pwmOutput = data.output_dir + "/" + data.TF_name + ".pwm";


    double avgScore = data.Signal_data.enumerate_maximum;
    avgScore *= avgScore;
    // second power
    if(avgScore == 0){
        cerr << "avgScore shouldn't be 0!!!!!!\n\tEXITING\n";
        exit(1);
    }

    vector<double> AA;
    vector<double> TT;
    vector<double> CC;
    vector<double> GG;
    vector<int> A;
    vector<int> T;
    vector<int> C;
    vector<int> G;

    ifstream INF(input);

    if(INF){
        cerr << "INF is not good!!!!!\n\tEXITING\n";
        exit(1);
    }

    string line = "";
    vector<string> fields_str;
    vector<double> fields_dbl;

    double minimumScore = 1000.0;

    map<char, vector<double> > pwm;
    parse_pwm(data.PWM_file, pwm);

    while(getline(INF, line, '\n')){

        // split line into fields_str by tab
        split_string(line, "\t", fields_str);
        // resize the double vector to match number of fields in line
        fields_dbl.resize(fields_str.size());

        for(size_t i = 1; i < 5; ++i){
            try{
                // copied from line 62 of perl version
                fields_dbl.at(i) = avgScore * pow(2, stod(fields_str.at(i) ) );
            }
            catch(...){
                cerr << "Something has been caught!!\n\tEXITING\n";
                exit(1);
            }
            // if a field is below minimum score, adjust minimum score
            if(fields_dbl[i] < minimumScore){
                minimumScore = fields_dbl[i];
            }
        }
        // copied directly from line 69 of perl version
        AA.push_back(fields_dbl[1]);
        CC.push_back(fields_dbl[2]);
        GG.push_back(fields_dbl[3]);
        TT.push_back(fields_dbl[4]);
        fields_str.clear();
    }
    INF.close();

    double rowmin = 0, denom = 0.0;
    for(size_t i = 0; i < AA.size(); ++i){
        rowmin = minimumScore;
        denom = 0.0;
        denom += (AA[i] - rowmin) / (avgScore - rowmin);
        denom += (CC[i] - rowmin) / (avgScore - rowmin);
        denom += (GG[i] - rowmin) / (avgScore - rowmin);
        denom += (TT[i] - rowmin) / (avgScore - rowmin);

        A.push_back(static_cast<int>( ( AA[i] - rowmin )
                                    / ( avgScore - rowmin )
                                    * ( 1000.0 / denom ) + 0.5 ) );
        C.push_back(static_cast<int>( ( CC[i] - rowmin )
                                    / ( avgScore - rowmin )
                                    * ( 1000.0 / denom ) + 0.5 ) );
        G.push_back(static_cast<int>( ( GG[i] - rowmin )
                                    / ( avgScore - rowmin )
                                    * ( 1000.0 / denom ) + 0.5 ) );
        T.push_back(static_cast<int>( ( TT[i] - rowmin )
                                    / ( avgScore - rowmin )
                                    * ( 1000.0 / denom ) + 0.5 ) );
    }

    // p-value
    double alpha = 0.05;

    ofstream OUTF(output);

    if(!OUTF){
        cerr << "OUTF is not good!!!!!\n\tEXITING\n";
        exit(1);
    }

    // formatting
    OUTF << "DE\t" + data.TF_name + '\n';

    int rowsum = 0;
    #ifdef DEBUG
        assert(A.size() == C.size());
        assert(A.size() == G.size());
        assert(A.size() == T.size());
    #endif
    for(size_t i = 0; i < A.size(); ++i){
        rowsum = A[i] + C[i] + G[i] + T[i];

        OUTF << i + '\t';
        OUTF << static_cast<int>( (static_cast<double>(A[i]) * alpha)
                                    + (static_cast<double>(pwm['A'][i] * rowsum * (1 - alpha) ))
                                    + .5 );
        OUTF << '\t';
        OUTF << static_cast<int>( (static_cast<double>(C[i]) * alpha)
                                    + (static_cast<double>(pwm['C'][i] * rowsum * (1 - alpha) ))
                                    + .5 );
        OUTF << '\t';
        OUTF << static_cast<int>( (static_cast<double>(G[i]) * alpha)
                                    + (static_cast<double>(pwm['G'][i] * rowsum * (1 - alpha) ))
                                    + .5 );
        OUTF << '\t';
        OUTF << static_cast<int>( (static_cast<double>(T[i]) * alpha)
                                    + (static_cast<double>(pwm['T'][i] * rowsum * (1 - alpha) ))
                                    + .5 );
        OUTF << "\tX\n";
    }

    OUTF << "XX\n";
    OUTF.close();


}

// REQUIRES: PWM within pwm file is formatted according to example
static void parse_pwm(const string &pwm, map<char, vector<double> > &motif){
    // clear motif for sanity
    motif.clear();
    // pwm is the filename containing a position weight matrix
    ifstream IN_HANDLE(pwm);

    if(!IN_HANDLE){
        cerr << "IN_HANDLE is not good!!!!\n\tEXITING\n";
        exit(1);
    }

    string line = "";
    vector<string> fields;
    vector<int> fields_int;
    // while getline works successfully
    while(getline(IN_HANDLE, line, '\n')){
        // getline removes the newline character
        // if neither of those strings are present in the current line
        if( ( line.find("DE") == string::npos )
         && ( line.find("XX") == string::npos ) ){
            split_string(line, "\t", fields);

            for(size_t i = 0; i < fields.size(); ++i){
                fields_int.push_back(stoi(fields[i]));
            }

            int rowsum = fields_int[1] + fields_int[2]
                       + fields_int[3] + fields_int[4];
            // fill motif, which is a map from char to vector,
            // stores integers corresponding to a char
            motif['A'].push_back(static_cast<double>(fields_int[1]) /
                                 static_cast<double>(rowsum) );
            motif['C'].push_back(static_cast<double>(fields_int[2]) /
                                 static_cast<double>(rowsum) );
            motif['G'].push_back(static_cast<double>(fields_int[3]) /
                                 static_cast<double>(rowsum) );
            motif['T'].push_back(static_cast<double>(fields_int[4]) /
                                 static_cast<double>(rowsum) );

            // clear for next line of pwm
            fields.clear();
            fields_int.clear();
        }
    }
}

