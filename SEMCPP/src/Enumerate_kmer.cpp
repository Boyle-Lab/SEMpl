#include "../iterativeSEM.hpp"
#include <map>
#include <cmath>
#include <utility>
#include <sstream>
#include <iostream>
#include <exception>
#include <stdlib.h>
using namespace std;

/*
* Does this file use the original PWM? Is new_PWM something else????
*   By Cody
*/

// custom comparator necessary for map
// is a function object
class hash_comp{
public:
  bool operator()(const pair<int, string> &one, const pair<int, string> &two) const {
    if(one.first == two.first){
      return one.second < two.second;
    }
    return one.first < two.first;
  }
};

// EFFECTS: calculates log base 2 of x
double log2(double x){
  return log(x) / log(2);
}

// EFFECTS: finds maximum mapped value
int findMax(const map<int, double> &v){
  int max = 0;
  for(auto i : v){
    if(i.second > max)
      max = i.second;
  }
  return max;
}

static void create_kmer(const Dataset &data,
                        map<pair<int, string>, int, hash_comp> &pwmHash,
                        vector<string> &nucleotideStack,
                        vector<int> &bestCase, map<string, int> &retHash,
                        const double cutoff){
                        if(data.settings.verbose){
                          cout << "Generating kmer list.\n";
                        }
                        { // test block begin
                          size_t check = pwmHash.size();
                          for(size_t i = 0; i < nucleotideStack.size(); i++){
                            retHash[nucleotideStack[i]] = pwmHash[{1, nucleotideStack[i]}];
                          }
                          if(check != pwmHash.size()){
                            cerr << "size of pwmHash was modified within create_kmer\n\tEXITING\n";
                            exit(EXIT_FAILURE);
                          }
                        } // test block end

                        vector<int> maxScores(bestCase.size() + 1);
                        maxScores[bestCase.size()] = 0;
                        for(int i = static_cast<int>(bestCase.size()) - 1; i >= 0; i--){
                          maxScores[i] = maxScores[i+1] + bestCase[i];
                        }

                          // watch to debug below here
                          for(size_t i = 1; i < bestCase.size(); i++){
                            for(auto pair : retHash){
                              int score = retHash[pair.first];
                              int length = static_cast<int>(pair.first.size());
                              int maxscore = maxScores[length] + score;
                              if(maxscore < cutoff){
                                retHash.erase(pair.first);
                              }
                              else{
                                retHash.erase(pair.first);
                                for(size_t j = 0; j < nucleotideStack.size(); j++){
                                  string newkey = pair.first + nucleotideStack[j];
                                  int newscore = score + pwmHash[{i + 1, nucleotideStack[j]}];
                                  retHash[newkey] = newscore;
                                }
                              }
                            }
                          // debug above here
                          }
                      }

static double get_cutoff(const Dataset &data,
                      map<pair<int, string>, int, hash_comp> &pwmHash,
                      vector<string> &nucleotideStack,
                      vector<int> &bestCase, map<string, int> &retHash){
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
      // assumes that there is whitespace between the name and the value
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

static void parse_pwm(const Dataset &data,
                      map<pair<int, string>, int, hash_comp> &pwmHash,
                      vector<string> &nucleotideStack,
                      vector<int> &bestCase, map<string, int> &retHash){


  nucleotideStack.push_back("A");
  nucleotideStack.push_back("C");
  nucleotideStack.push_back("G");
  nucleotideStack.push_back("T");

  {
    // modifiedFields is indexed from 1 in original implementation
    // in original implementation, fields[0] contains number of row
    // first row with numbers is indexed from 0, thus the last row is 12
    // but there are 13 total rows
    int sum = 0;
    for(int i = 0; i < Dataset::PWM::NUM_ROWS; i++){
      map<int, double> modifiedFields;
      for(int j = 0; j < Dataset::PWM::NUM_COLUMNS; j++){
        sum += data.PWM_data.matrix_arr[i][j];
        modifiedFields[j+1] = log2((data.PWM_data.matrix_arr[i][j] + 0.25) / (sum + 1)) - log2(0.25);
      }
      // TEST ME!!!
      // TEST ME!!!
      // TEST ME!!!
      int max = findMax(modifiedFields);
      bestCase.push_back(max);
      for(int k = 1; k < Dataset::PWM::NUM_COLUMNS; k++){
        pwmHash[{data.PWM_data.matrix_arr[i][0] + 1, nucleotideStack[k - 1]}] = modifiedFields[k];
      }
    }
  }
}

// REQUIRES: within data, PWM_data is filled in
// EFFECTS: created enumerated kmers using a cutoff and a PWM matrix, returns the map
// note: for now, searches for a pre-calculated cutoff
void Enumerate_kmer(Dataset &data){
  map<pair<int, string>, int, hash_comp> pwmHash;
  vector<string> nucleotideStack;
  vector<int> bestCase;
  map<string, int> retHash;
  double cutoff = 0.0;

  parse_pwm(data, pwmHash, nucleotideStack, bestCase, retHash);
  // pwmHash is now filled in, along with bestCase
  if(data.settings.verbose){
    cout << "No cutoff defined, so earching for pre-calculated cutoff." << '\n';
  }
  cutoff = get_cutoff(data, pwmHash, nucleotideStack, bestCase, retHash);

  if(cutoff == 0.0){
    cerr << "cutoff value unchanged within Enumerate_kmer.cpp\n\tEXITING\n";
    exit(EXIT_FAILURE);
  }

  create_kmer(data, pwmHash, nucleotideStack, bestCase, retHash, cutoff);

  // print_kmer is replaced by a simple object assignment
  for(auto pair : retHash){
    if(pair.second > cutoff){
      retHash.erase(pair.first);
    }
  }
  data.kmerHash = retHash;
}
