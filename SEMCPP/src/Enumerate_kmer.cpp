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
                        const vector<double> &bestCase,
                        map<string, double> &retHash,
                        const double cutoff);
static void parse_pwm(const Dataset &data,
                      map<pair<int, string>, int, hash_comp> &pwmHash,
                      const vector<string> &nucleotideStack,
                      vector<double> &bestCase);

// REQUIRES: within data, PWM_data is filled in
// EFFECTS: created enumerated kmers using a cutoff and a PWM matrix, returns the map
// note: for now, searches for a pre-calculated cutoff
void Enumerate_kmer(Dataset &data){
    map<pair<int, string>, int, hash_comp> pwmHash;
    const vector<string> nucleotideStack {"A", "C", "G", "T"};
    vector<double> bestCase;
    double cutoff = 0.0;

    try{
        parse_pwm(data, pwmHash, nucleotideStack, bestCase);
    }
    catch(...){
        cerr << "exception thrown from parse_pwm" << endl;
        exit(1);
    }
  // pwmHash is now filled in, along with bestCase
    if(data.settings.verbose){
        cout << "\tNo cutoff defined, so searching for pre-calculated cutoff." << '\n';
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
    std::vector<std::string> to_erase;
    for(auto pair : data.kmerHash){
        // cout << pair.first << ' ' << pair.second << ' ';
        if(pair.second <= cutoff){
            // data.kmerHash.erase(pair.first);
            to_erase.push_back(pair.first);
            // cout << "erased!";
        }
        else{
            // cout << "not erased!";
        }
        // cout << endl;
    }
    for(auto val : to_erase){
        data.kmerHash.erase(val);
    }
    // data.kmerHash is now ready for use

#ifdef DEBUG
    data.size_of_kmerHash = data.kmerHash.size();
#endif
}

static void parse_pwm(const Dataset &data,
                      map<pair<int, string>, int, hash_comp> &pwmHash,
                      const vector<string> &nucleotideStack,
                      vector<double> &bestCase){

  // modifiedFields is indexed from 1 in original implementation
  // in original implementation, fields[0] contains number of row
  // first row with numbers is indexed from 0, thus the last row is 12
  // but there are 13 total rows

#ifdef DEBUG
    cout << "\tparse_pwm!" << endl;
#endif

    // cout << Dataset::PWM::NUM_ROWS << ' ' << Dataset::PWM::NUM_COLUMNS << endl;

    for(int i = 0; i < Dataset::PWM::NUM_ROWS; ++i){
        int sum = 0;
        map<int, double> modifiedFields;
        for(int j = 0; j < Dataset::PWM::NUM_COLUMNS; ++j){
            sum += data.PWM_data.matrix_arr[i][j];
#ifdef DEBUG
            // cout << data.PWM_data.matrix_arr[i][j] << endl;;
#endif
            modifiedFields[j+1] = log2((static_cast<double>(data.PWM_data.matrix_arr[i][j]) + 0.25)
                                     / (sum + 1)) - log2(0.25);
        }
#ifdef DEBUG
        // for(auto val : modifiedFields) cout << val.first << ' ' << val.second << endl;
        //     cout << endl;

        size_t sz = modifiedFields.size();
#endif
        // TEST ME!!!
        // TEST ME!!!
        // TEST ME!!!
        bestCase.push_back(findMax(modifiedFields));
        for(int k = 1; k < Dataset::PWM::NUM_COLUMNS + 1; ++k){
            pwmHash[{i + 1, nucleotideStack[k - 1]}] = modifiedFields.at(k);
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
        //  << endl;
        // for(auto val : pwmHash) cout << val.first.first << ' ' 
        //                              << val.first.second << ' ' << val.second 
        //                              << endl;
    #endif

}

static void create_kmer(const Dataset &data,
                        const map<pair<int, string>, int, hash_comp> &pwmHash,
                        const vector<string> &nucleotideStack,
                        const vector<double> &bestCase,
                        map<string, double> &retHash,
                        const double cutoff){
    if(data.settings.verbose){
        cout << "Generating kmer list.\n";
    }
    for(size_t i = 0; i < nucleotideStack.size(); ++i){
        try{
            retHash[nucleotideStack[i]] = pwmHash.at(std::make_pair(1, nucleotideStack[i]));
        }
        catch(...){
            cerr << "line 166 Enumerate_kmer.cpp " << i << nucleotideStack[i] << endl;
            exit(1);
        }
    }

    vector<double> maxScores(bestCase.size() + 1);
    maxScores.at(bestCase.size()) = 0;
    for(int i = static_cast<int>(bestCase.size()) - 1; i >= 0; --i){
        maxScores.at(i) = maxScores.at(i+1) + bestCase.at(i);
#ifdef DEBUG
      // cout << bestCase[i] << endl;
#endif
    }
#ifdef DEBUG
    // cout << bestCase.size() << endl;
    // for(auto val : maxScores) cout << val << endl;
#endif

  // watch to debug below here

#ifdef DEBUG
    // for(auto val : retHash) cout << val.first << ' ' << val.second << endl;
#endif

// modifying object while iterating through
// will invalidate iterators

    std::vector<std::string> to_erase;
    std::vector<std::pair<std::string, int> > to_add;
    for(size_t i = 1; i < bestCase.size(); ++i){
        for(auto pair : retHash){
            double score = pair.second;
            int length = static_cast<int>(pair.first.size());
            double maxscore = maxScores.at(length) + score;
	    //cout << pair.first << endl;
            try{
                //retHash.erase(pair.first);
                to_erase.push_back(pair.first);
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
                    double newscore = score + pwmHash.at({i + 1, nucleotideStack[j]});
                    // retHash[newkey] = newscore;
                    to_add.push_back({newkey, newscore});
                }
            }
        }
  // debug above here
    	for(auto val : to_erase){
        	retHash.erase(val);
    	}
    	for(auto pair : to_add){
        	retHash[pair.first] = pair.second;
    	}
    }
}

static double get_cutoff(const Dataset &data){
    if(data.settings.verbose){
        cout << "\tSearching for pre-calculated cutoff" << '\n';
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
              stringstream parse; // seems good
              parse << curr_line; 
              parse >> curr_line; 
              double a = 0.0;
              parse >> a;
              if(data.settings.verbose){
                  cout << "\tPre-calculated threshold identified as " << a << '\n';
              }
          IN_HANDLE.close();
          return a;
        }
    }
    cerr << "Unable to find pre-caluclated cutoff in file\n\tEXITING\n";
    exit(1);
    return 0.0;
}

// EFFECTS: finds maximum mapped value
static double findMax(const map<int, double> &v){
    double max = -1000000000.0;
    // cout << "here" << endl;
    for(auto i : v){
        // cout << i.first << ' ' << i.second << endl;
        if(i.second > max){
            max = i.second;
        }
    }
    // cout << "max: " << max << endl;
    return max;
}

