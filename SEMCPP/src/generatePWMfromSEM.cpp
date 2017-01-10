
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

using namespace std;


//REQUIRES: requires data is a valid Dataset with pwm filled in
//MODIFIES: modifies data, specifically pwm
//EFFECTS: Creates a pwm matrix from the previous sem matrix
void generatePWMfromSEM(Dataset & data){
    int minimumScore = 1000;
    int avgScore = 0;
    string rawBaseline = data.output_dir + "/BASELINE/baseline.maximums";
    ifstream inf(rawBaseline);
    string temp;
    vector<string> line;
    while (inf >> temp){
           line.push_back(temp);
    }
    if (line[0] == "Enumerated_kmer_filtered.signal"){
        int tempint = atoi(line[1].c_str());
        avgScore = pow(tempint, 2);
    }
    inf.close();

    vector<int> AA;
    vector<int> TT;
    vector<int> CC;
    vector<int> GG;
    vector<int> A;
    vector<int> T;
    vector<int> C;
    vector<int> G;


    for(int i=0; i < data.PWM_data.NUM_ROWS; i++){
            AA.push_back(data.PWM_data.matrix_arr[1][i]);
    }
    for(int i=0; i < data.PWM_data.NUM_ROWS; i++){
            CC.push_back(data.PWM_data.matrix_arr[1][i]);
    }
    for(int i=0; i < data.PWM_data.NUM_ROWS; i++){
            GG.push_back(data.PWM_data.matrix_arr[1][i]);
    }
    for(int i=0; i < data.PWM_data.NUM_ROWS; i++){
            TT.push_back(data.PWM_data.matrix_arr[1][i]);
    }

    for(int i = 0; i < data.PWM_data.NUM_ROWS; i++){
        int rowmin = minimumScore;
        int denom =0;
        denom += (AA[i] - rowmin)/(avgScore-rowmin);
        denom += (CC[i] - rowmin)/(avgScore-rowmin);
        denom += (GG[i] - rowmin)/(avgScore-rowmin);
        denom += (TT[i] - rowmin)/(avgScore-rowmin);

        A.push_back(round((AA[i]-rowmin)/(avgScore-rowmin) * (1000/denom) + 0.5));
        C.push_back(round((CC[i]-rowmin)/(avgScore-rowmin) * (1000/denom) + 0.5));
        G.push_back(round((GG[i]-rowmin)/(avgScore-rowmin) * (1000/denom) + 0.5));
        T.push_back(round((TT[i]-rowmin)/(avgScore-rowmin) * (1000/denom) + 0.5));
    }

    double alpha = 0.5;
    for (int i = 0; i < A.size(); i++){
        int rowsum = A[i] + C[i] + G[i] + T[i];
        data.PWM_data.sem_arr[1][i] = ((A[i] * alpha) + (data.PWM_data.matrix_arr[1][i] * rowsum * (1-alpha)) + .5);
        data.PWM_data.sem_arr[2][i] = ((A[i] * alpha) + (data.PWM_data.matrix_arr[2][i] * rowsum * (1-alpha)) + .5);
        data.PWM_data.sem_arr[3][i] = ((A[i] * alpha) + (data.PWM_data.matrix_arr[3][i] * rowsum * (1-alpha)) + .5);
        data.PWM_data.sem_arr[4][i] = ((A[i] * alpha) + (data.PWM_data.matrix_arr[4][i] * rowsum * (1-alpha)) + .5);
    }
}
