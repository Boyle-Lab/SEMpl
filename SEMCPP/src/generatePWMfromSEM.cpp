
#include "../iterativeSEM.hpp"
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

//REQUIRES: (please fill in what the function requires, before we can call the function, such as certain data being filled in
//MODIFIES: modifies data
//EFFECTS: (brief statement of what the function does)
map <char, int> parse_pwm(Dataset & data){
    map <char, int> motif;
    int rowsum = 0;
    int counter = 0;
    while (counter < data.PWM_data.MATRIX_SIZE){
        for (int i = 1; i < data.PWM_data.ROW_SIZE; i++){
            rowsum += data.PWM_data.matrix_arr[i][counter];
        }
        motif['A',counter ] = data.PWM_data.matrix_arr[1][counter]/rowsum;
        motif['C',counter ] = data.PWM_data.matrix_arr[2][counter]/rowsum;
        motif['G',counter ] = data.PWM_data.matrix_arr[3][counter]/rowsum;
        motif['T',counter ] = data.PWM_data.matrix_arr[4][counter]/rowsum;
        counter++;
    }
    return motif;
}

//REQUIRES: (please fill in what the function requires, before we can call the function, such as certain data being filled in
//MODIFIES: modifies data
//EFFECTS: (brief statement of what the function does)
void generatePWMfromSEM(Dataset & data, string output){
    int minimumScore = 1000;
    int avgScore = 0;
    string rawBaseline = output + "/BASELINE/baseline.maximums";
    map <char, int> pwm = parse_pwm(data);
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
}

    //Incomplete and builds with error in reference to the array in parse_pwm



