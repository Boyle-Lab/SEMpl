
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
//EFFECTS: Creates a pwm matrix from the previous sem matrix
void generatePWMfromSEM(Dataset & data){

    // raw baseline is the output from running findMaximumAverageSignal... on
    // "enumerate" data

    string rawInput = data.output_dir + "/" + data.TF_name + ".sem";
    string pwmOutput = data.output_dir + "/" + data.TF_name + ".pwm";

    double avgScore = data.signal_Data.enumerate_maximum;
    avgScore = avgScore * avgScore;
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

    ifstream INF(rawInput);

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

        split_string(line, "\t", fields_str);
        fields_dbl.resize(fields_str.size());

        for(size_t i = 1; i < 5; ++i){
            try{
                fields_dbl[i] = avgScore * pow(2, stod(fields_str[i]));
            }
            catch(...){
                cerr << "Something has been caught!!\n\tEXITING\n";
                exit(1);
            }

            if(fields_dbl[i] < minimumScore){
                minimumScore = fields_dbl[i];
            }
        }

        AA.push_back(fields_dbl[1]);
        CC.push_back(fields_dbl[2]);
        GG.push_back(fields_dbl[3]);
        TT.push_back(fields_dbl[4]);
        fields_str.clear();
    }
    INF.close();


    for(size_t i = 0; i < AA.size(); ++i){
        double rowmin = minimumScore;
        double denom = 0.0;
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

    double alpha = 0.05;

    ofstream OUTF(pwmOutput);

    if(!OUTF){
        cerr << "OUTF is not good!!!!!\n\tEXITING\n";
        exit(1);
    }

    OUTF << "DE\t" + data.TF_name + '\n';

    for(size_t i = 0; i < A.size(); ++i){
        int rowsum = A[i] + C[i] + G[i] + T[i];

        OUTF << i + '\t';
        OUTF << static_cast<int>( (static_cast<double>(A[i]) * alpha)
                                    + (static_cast<double>(pwm['A'][i] * rowsum * (1- alpha) ))
                                    + .5 );
        OUTF << i + '\t';
        OUTF << static_cast<int>( (static_cast<double>(C[i]) * alpha)
                                    + (static_cast<double>(pwm['C'][i] * rowsum * (1- alpha) ))
                                    + .5 );
        OUTF << i + '\t';
        OUTF << static_cast<int>( (static_cast<double>(G[i]) * alpha)
                                    + (static_cast<double>(pwm['G'][i] * rowsum * (1- alpha) ))
                                    + .5 );
        OUTF << i + '\t';
        OUTF << static_cast<int>( (static_cast<double>(T[i]) * alpha)
                                    + (static_cast<double>(pwm['T'][i] * rowsum * (1- alpha) ))
                                    + .5 );
        OUTF << "\tX\n";
    }

    OUTF << "XX\n";
    OUTF.close();


}

static void parse_pwm(const string &pwm, map<char, vector<double> > &map){
    ifstream IN_HANDLE(pwm);

    if(!IN_HANDLE){
        cerr << "IN_HANDLE is not good!!!!\n\tEXITING\n";
        exit(1);
    }

    string line = "";
    vector<string> fields;
    vector<int> fields_int;
    while(getline(IN_HANDLE, line, '\n')){
        // getline removes the newline character
        // if neither of those strings are present in the current line
        if( ( line.find("DE") == string::npos )
         && ( line.find("XX") == string::npos ) ){
            split_string(line, "\t", fields);

            fields_int.push_back(-1);

            for(size_t i = 0; i < fields.size(); ++i){
                if(i == 0) continue;
                fields_int.push_back(stoi(fields[i]));
            }

            int rowsum = fields_int[1] + fields_int[2]
                       + fields_int[3] + fields_int[4];

            map['A'].push_back(static_cast<double>(fields_int[1]) /
                               static_cast<double>(rowsum) );
            map['C'].push_back(static_cast<double>(fields_int[2]) /
                               static_cast<double>(rowsum) );
            map['G'].push_back(static_cast<double>(fields_int[3]) /
                               static_cast<double>(rowsum) );
            map['T'].push_back(static_cast<double>(fields_int[4]) /
                               static_cast<double>(rowsum) );

            fields.clear();
            fields_int.clear();
        }
    }
}

