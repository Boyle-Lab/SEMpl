
#include "iterativeSEM.hpp"
#include <iostream>
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

#ifdef DEBUG
    assert(INF);
#endif

    string line = "";
    vector<string> fields_str;
    vector<double> fields_dbl;

    double minimumScore = 1000.0;

    //map<char, vector<double> > pwm;
    //parse_pwm(data.PWM_file, pwm); // we already have this in memory

    //remove first line from file
    getline(INF, line, '\n');

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

    double rowsum = 0.0;
    double rowsum_pwm = 0.0;
    #ifdef DEBUG
        assert(A.size() == C.size());
        assert(A.size() == G.size());
        assert(A.size() == T.size());
    #endif
    for(size_t i = 0; i < A.size(); ++i){
        rowsum = A[i] + C[i] + G[i] + T[i];
        rowsum_pwm = data.PWM_data.matrix_arr[0][i] + data.PWM_data.matrix_arr[1][i] + data.PWM_data.matrix_arr[2][i] + data.PWM_data.matrix_arr[3][i];


//	OUTF << "A: " << A[i] << " MatA " << data.PWM_data.matrix_arr[0][i] << " rowsump: " << rowsum_pwm << " rowsum: " << rowsum << endl;
        OUTF << i << '\t';
        OUTF << static_cast<int>( (static_cast<double>(A[i]) * alpha)
                                    + (static_cast<double>(data.PWM_data.matrix_arr[0][i]) / rowsum_pwm * rowsum * (1 - alpha) )
                                    + .5 );
        OUTF << '\t';
        OUTF << static_cast<int>( (static_cast<double>(C[i]) * alpha)
                                    + (static_cast<double>(data.PWM_data.matrix_arr[1][i]) / rowsum_pwm * rowsum * (1 - alpha) )
                                    + .5 );
        OUTF << '\t';
        OUTF << static_cast<int>( (static_cast<double>(G[i]) * alpha)
                                    + (static_cast<double>(data.PWM_data.matrix_arr[2][i]) / rowsum_pwm * rowsum * (1 - alpha) )
                                    + .5 );
        OUTF << '\t';
        OUTF << static_cast<int>( (static_cast<double>(T[i]) * alpha)
                                    + (static_cast<double>(data.PWM_data.matrix_arr[3][i]) / rowsum_pwm * rowsum * (1 - alpha) )
                                    + .5 );
        OUTF << "\tX\n";
    }

    OUTF << "XX\n";
    OUTF.close();


}

