
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

using namespace std;



map <char, int> parse_pwm(Dataset & data){
    map <char, int> motif;
    int rowsum;
    int counter = 0;
    while (counter < data.PWM_data.MATRIX_SIZE){
        for (int i = 1; i < data.PWM_data.ROW_SIZE; i++){
            rowsum += data.PWM_data.matrix_arr(i, counter);
        }
        motif['A',counter ] = data.PWM_data.matrix_arr((1), counter)/rowsum;
        motif['C',counter ] = data.PWM_data.matrix_arr((2), counter)/rowsum;
        motif['G',counter ] = data.PWM_data.matrix_arr((3), counter)/rowsum;
        motif['T',counter ] = data.PWM_data.matrix_arr((4), counter)/rowsum;
        counter++;
    }
    return motif;
}

void generatePWMfromSEM(Dataset & data, string output){
    int minimumScore = 1000;
    int avgScore = 0;
    map <char, int> pwm = parse_pwm(data);
    ifstream inf(output);
    string line;
    while (inf >> line){

        }
    }
    //Incomplete and builds with error in reference to the array in parse_pwm


