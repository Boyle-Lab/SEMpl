#include "iterativeSEM.hpp"
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

static double findMax(const map<int, double> &v);
static double get_cutoff(const Dataset &data);
static void create_kmer(const Dataset &data,
                        const map<pair<int, string>, int, hash_comp> &pwmHash,
                        const vector<string> &nucleotideStack,
                        const vector<int> &bestCase,
                        map<string, int> &retHash,
                        const double cutoff);
static void parse_pwm(const Dataset &data,
                      map<pair<int, string>, int, hash_comp> &pwmHash,
                      const vector<string> &nucleotideStack,
                      vector<int> &bestCase);

// REQUIRES: within data, PWM_data is filled in
// EFFECTS: created enumerated kmers using a cutoff and a PWM matrix, returns the map
// note: for now, searches for a pre-calculated cutoff
void Enumerate_kmer(Dataset &data){
    map<pair<int, string>, int, hash_comp> pwmHash;
    vector<string> nucleotideStack {"A", "C", "G", "T"};
    vector<int> bestCase;
    double cutoff = 0.0;

    parse_pwm(data, pwmHash, nucleotideStack, bestCase);
  // pwmHash is now filled in, along with bestCase
    if(data.settings.verbose){
        cout << "No cutoff defined, so earching for pre-calculated cutoff." << '\n';
    }
    cutoff = get_cutoff(data);

    if(cutoff == 0.0){
        cerr << "cutoff value unchanged within Enumerate_kmer.cpp\n\tEXITING\n";
        exit(1);
    }
    try{
        create_kmer(data, pwmHash, nucleotideStack, bestCase, data.kmerHash, cutoff);
    }
    catch(...){
        cerr << "Problem with create_kmer!\n\tEXITING" << endl;
        exit(1);
    }
    // print_kmer is replaced by a simple object assignment
    for(auto pair : data.kmerHash){
        if(pair.second <= cutoff){
            data.kmerHash.erase(pair.first);
        }
    }
    // data.kmerHash is now ready for use

#ifdef DEBUG
    data.size_of_kmerHash = data.kmerHash.size();
#endif
}

// EFFECTS: finds maximum mapped value
static double findMax(const map<int, double> &v){
    int max = 0;
    for(auto i : v){
        if(i.second > max){
            max = i.second;
        }
    }
    return max;
}

static void parse_pwm(const Dataset &data,
                      map<pair<int, string>, int, hash_comp> &pwmHash,
                      const vector<string> &nucleotideStack,
                      vector<int> &bestCase){

  // modifiedFields is indexed from 1 in original implementation
  // in original implementation, fields[0] contains number of row
  // first row with numbers is indexed from 0, thus the last row is 12
  // but there are 13 total rows

#ifdef DEBUG
    cout << "\tparse_pwm!" << endl;
#endif

    // cout << Dataset::PWM::NUM_ROWS << ' ' << Dataset::PWM::NUM_COLUMNS << endl;

    for(int i = 0; i < Dataset::PWM::NUM_ROWS; i++){
        int sum = 0;
        map<int, double> modifiedFields;
        for(int j = 0; j < Dataset::PWM::NUM_COLUMNS; j++){
            sum += data.PWM_data.matrix_arr[i][j];
#ifdef DEBUG
            // cout << data.PWM_data.matrix_arr[i][j] << endl;;
#endif
            modifiedFields[j+1] = log2((static_cast<double>(data.PWM_data.matrix_arr[i][j]) + 0.25) / (sum + 1)) - log2(0.25);
        }
#ifdef DEBUG
        for(auto val : modifiedFields) cout << val.first << ' ' << val.second << endl;
            cout << endl;

        size_t sz = modifiedFields.size();
#endif
        // TEST ME!!!
        // TEST ME!!!
        // TEST ME!!!
        bestCase.push_back(findMax(modifiedFields));
        for(int k = 1; k < Dataset::PWM::NUM_COLUMNS + 1; k++){
            pwmHash[{data.PWM_data.matrix_arr[i][0] + 1, nucleotideStack[k - 1]}] = modifiedFields[k];
        }
#ifdef DEBUG
        if(sz != modifiedFields.size()){
            cerr << "modifiedFields has changed in size!!\n\tEXITING\n";
            exit(1);
        }
        if(sz == 0){
            cerr << "sz shouldn't be 0!!\n\tEXITING";
            exit(1);
        }
#endif
    }
    #ifdef DEBUG

        cout << endl;
        for(auto val : pwmHash) cout << val.first.first << ' ' << val.first.second << ' ' << val.second << endl;
    #endif

}

static void create_kmer(const Dataset &data,
                        const map<pair<int, string>, int, hash_comp> &pwmHash,
                        const vector<string> &nucleotideStack,
                        const vector<int> &bestCase,
                        map<string, int> &retHash,
                        const double cutoff){
    if(data.settings.verbose){
        cout << "Generating kmer list.\n";
    }
    for(size_t i = 0; i < nucleotideStack.size(); i++){
        try{
            retHash[nucleotideStack[i]] = pwmHash.at({1, nucleotideStack[i]});
        }
        catch(...){
            cerr << "line 103 Enumerate_kmer.cpp" << endl;
            exit(1);
        }
    }

    vector<int> maxScores(bestCase.size() + 1);
    maxScores.at(bestCase.size()) = 0;
    for(int i = static_cast<int>(bestCase.size()) - 1; i >= 0; i--){
      maxScores[i] = maxScores[i+1] + bestCase[i];
    }

  // watch to debug below here

    for(size_t i = 1; i < bestCase.size(); ++i){
        for(auto pair : retHash){
            int score = pair.second;
            int length = static_cast<int>(pair.first.size());
            int maxscore = maxScores[length] + score;
            try{
                retHash.erase(pair.first);
            }
            catch(...){
                cerr << "line 125 Enumerate_kmer.cpp" << endl;
                exit(1);
            }
            if(maxscore < cutoff){
                
            }
            else{
                for(size_t j = 0; j < nucleotideStack.size(); j++){
                    string newkey = pair.first + nucleotideStack[j];
                    int newscore = score + pwmHash.at({i + 1, nucleotideStack[j]});
                    retHash[newkey] = newscore;
                }
            }
        }
  // debug above here
    }
}

static double get_cutoff(const Dataset &data){
    if(data.settings.verbose){
        cout << "Searching for pre-calculated cutoff" << '\n';
    }

      string curr_line = "";
      // will need to change below line eventually to be specified
      ifstream IN_HANDLE("src/PWM_SCORES_FINAL.txt");

    if(!IN_HANDLE){
        cerr << "Failure to open src/PWM_SCORES_FINAL.txt\n\tEXITING" << '\n';
        exit(1);
    }

    while(getline(IN_HANDLE, curr_line)){
        if(curr_line.find(data.TF_name) != string::npos){
          // assumes that there is whitespace between the name and the value
              stringstream parse; // NEEDS TO BE TESTED
              parse << curr_line; // NEEDS TO BE TESTED
              parse >> curr_line; // NEEDS TO BE TESTED
              double a;
              parse >> a;
              if(data.settings.verbose){
                  cout << "Pre-calculated threshold identified as " << a << '\n';
              }
          IN_HANDLE.close();
          return a;
        }
    }
    cerr << "Unable to find pre-caluclated cutoff in file\n\tEXITING\n";
    exit(1);
    return 0.0;
}
