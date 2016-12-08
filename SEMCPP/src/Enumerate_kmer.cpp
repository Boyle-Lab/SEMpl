#include "../iterativeSEM.hpp"
#include <map>
#include <cmath>
#include <utility>
#include <sstream>
#include <iostream>
#include <stdlib.h>
using namespace std;

/*
* Does this file use the original PWM? Is new_PWM something else????
*/


class hash_comp{
public:
  bool operator()(pair<char, int> one, pair<char, int> two){
    if(one.first == two.first){
      return one.second < two.second;
    }
    return one.first < two.second;
  }
};

int findMax(const map<int, double> &v){
  int max = 0;
  for(auto i : v){
    if(i.second > max)
      max = i.second;
  }
  return max;
}

static double get_cutoff(Dataset &data,
                      map<pair<char, int>, int, hash_comp> &pwmHash,
                      vector<char> &nucleotideStack,
                      vector<int> &bestCase, map<string, int> &kmerHash){
  if(data.settings.verbose){
    cout << "Searching for pre-calculated cutoff" << '\n';
  }

  string curr_line = "";
  ifstream IN_HANDLE("src/PWM_SCORES_FINAL.txt");

  if(!IN_HANDLE){
    cerr << "Failure to open src/PWM_SCORES_FINAL.txt\n\tEXITING" << '\n';
    exit(EXIT_FAILURE);
  }

  while(getline(IN_HANDLE, curr_line)){
    if(curr_line.find(data.TF_name) != string::npos){
      stringstream parse; // NEEDS TO BE TESTED
      parse << curr_line; // NEEDS TO BE TESTED
      parse >> curr_line; // NEEDS TO BE TESTED
      double a;
      parse >> a;
      return a;
    }
  }
  cerr << "Unable to find pre-caluclated cutoff in file\n\tEXITING\n";
  exit(EXIT_FAILURE);
  return 0.0;
}




static void parse_pwm(Dataset &data,
                      map<pair<char, int>, int, hash_comp> &pwmHash,
                      vector<char> &nucleotideStack,
                      vector<int> &bestCase, map<string, int> &kmerHash){
  nucleotideStack.push_back('A');
  nucleotideStack.push_back('C');
  nucleotideStack.push_back('G');
  nucleotideStack.push_back('T');

  {
    int sum = 0;

    for(int i = 0; i < Dataset::PWM::NUM_ROWS; i++){
      map<int, double> modifiedFields;
      for(int j = 0; j < Dataset::PWM::NUM_COLUMNS; j++){
        sum += data.PWM_data.matrix_arr[i][j];
        modifiedFields[log((data.PWM_data.matrix_arr[i][j] + 0.25)/ sum + 1) / log(2) - log(0.25) / log(2)];
      }
      int max = findMax(modifiedFields);
      bestCase.push_back(max);
      for(int k = 1; k < Dataset::PWM::NUM_COLUMNS; k++){

      }
    }
  }

}

// REQUIRES:
void Enumerate_kmer(Dataset &data){
  map<pair<char, int>, int, hash_comp> pwmHash;
  vector<char> nucleotideStack;
  vector<int> bestCase;
  map<string, int> kmerHash;

  parse_pwm(data, pwmHash, nucleotideStack, bestCase, kmerHash);
  get_cutoff(data, pwmHash, nucleotideStack, bestCase, kmerHash);

  if(data.settings.verbose){
    cout << "Using user defined cutoff." << '\n';
  }



}
