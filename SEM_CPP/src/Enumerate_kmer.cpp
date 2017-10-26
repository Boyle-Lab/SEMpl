#include "iterativeSEM.hpp"
#include <map>
#include <cmath>
#include <utility>
#include <sstream>
#include <iostream>
#include <exception>
#include <stdlib.h>
using namespace std;

// custom comparator necessary for map
// is a function object
// class hash_comp{
// public:
//     bool operator()(const pair<int, string> &one, const pair<int, string> &two) const {
//         if(one.first == two.first){
//             return one.second < two.second;
//         }
//         return one.first < two.first;
//     }
// };

static double findMax(const map<int, double> &v);
static double get_cutoff(const Dataset &data);
static void create_kmer(const Dataset &data,
                        const map<pair<int, char>, double> &pwmHash,
                        const vector<char> &nucleotideStack,
                        const vector<double> &bestCase,
                        map<std::string, double> &retHash,
                        const double cutoff);
static void parse_pwm(const Dataset &data,
                      map<pair<int, char>, double> &pwmHash,
                      vector<char> &nucleotideStack,
                      vector<double> &bestCase);

// REQUIRES: within data, PWM_data is filled in
// EFFECTS: created enumerated kmers using a cutoff and a PWM matrix, returns the map
//          the result is data.kmerHash
// note: for now, searches for a pre-calculated cutoff
void Enumerate_kmer(Dataset &data){
    map< pair<int, char>, double> pwmHash;
    //                      default nucleotide ordering,
    //                      will check in parse_pwm
    //                      for nucleotide ordering (columns) within pwm
    vector<char> nucleotideStack {'A', 'C', 'G', 'T'};
    vector<double> bestCase;
    double cutoff = 0.0;

    // clear data.kmerHash
    data.kmerHash.clear();

    try{
        #ifdef DEBUG
        cout << "\tparsing pwm..." << flush;
        #endif
        parse_pwm(data, pwmHash, nucleotideStack, bestCase);
        #ifdef DEBUG
        cout << "FINISH" << endl;
        #endif
    }
    catch(...){
        cerr << "exception thrown from parse_pwm" << endl;
        exit(1);
    }
  // pwmHash is now filled in, along with bestCase
    if(data.settings.threshold < 0.0){
        // indicated in iterativeSEM.hpp that if
        // data.settings.threshold < 0.0, then the value
        // should be treated as not defined, for purposes of
        // replicating the perl version
        if(data.settings.verbose){
            cout << "\tNo cutoff defined, so searching for pre-calculated cutoff.\n";
        }
        cutoff = get_cutoff(data);
    }
    else{
        if(data.settings.verbose){
            cout << "\tUsing user defined cutoff.\n";
        }
        cutoff = data.settings.threshold;
    }

    if(cutoff == 0.0){
        cerr << "cutoff value unchanged within Enumerate_kmer.cpp\n\tEXITING" << endl;
        exit(1);
    }
    try{
        create_kmer(data, pwmHash, nucleotideStack, bestCase, data.kmerHash, cutoff);
    }
    catch(...){
        cerr << "Problem with create_kmer!\n\tEXITING" << endl;
        exit(1);
    }

    for(auto iter = data.kmerHash.begin(); iter != data.kmerHash.end(); ){
        if(iter->second <= cutoff){
            data.kmerHash.erase(iter++);
        }
        else{
            ++iter;
        }
    }
    // data.kmerHash is now ready for use


#ifdef DEBUG
    ofstream OUT("Enumerated_kmers.txt");
    for(auto val : data.kmerHash){
        OUT << val.first << '\t' << val.second << endl;
    }
    data.size_of_kmerHash = data.kmerHash.size();
#endif
}

static void parse_pwm(const Dataset &data,
                      map<pair<int, char>, double> &pwmHash,
                      vector<char> &nucleotideStack,
                      vector<double> &bestCase){

  // modifiedFields is indexed from 1 in original implementation
  // in original implementation, fields[0] contains number of row
  // first row with numbers is indexed from 0, thus the last row is 12
  // but there are 13 total rows
    bestCase.clear();
    pwmHash.clear();


    assert(nucleotideStack.size() == data.PWM_data.matrix_arr.size());

    int sum = 0;

    // matrix_arr's first index is the column
    // matrix_arr's second index is the row

    for(int row = 0; row < (int)data.PWM_data.matrix_arr[0].size(); ++row){
        // cout << row << endl;
        sum = 0;
        map<int, double> modifiedFields;

        for(int column = 0; column < (int)data.PWM_data.matrix_arr.size(); ++column){
            // sums each row
            sum += data.PWM_data.matrix_arr[column][row];

        }
        // cout << "sum: " << sum << endl;
        for(int column = 0; column < (int)data.PWM_data.matrix_arr.size(); ++column){
            modifiedFields[column+1] = log2((static_cast<double>(data.PWM_data.matrix_arr[column][row]) + 0.25)
                                     / (static_cast<double>(sum) + 1.0)) - log2(0.25);
        }
#ifdef DEBUG
        size_t sz = modifiedFields.size();
#endif
        // TEST ME!!!
        // TEST ME!!!
        // TEST ME!!!
        bestCase.push_back(findMax(modifiedFields));
        try{
            for(int column = 0; column < (int)data.PWM_data.matrix_arr.size(); ++column){
                pwmHash[ { row + 1, nucleotideStack[column] } ] = modifiedFields.at(column + 1);
            }
        }
        catch(...){
            cerr << "out of range line 157" << endl;
            exit(1);
        }
#ifdef DEBUG
        if(sz != modifiedFields.size()){
            cerr << "modifiedFields has changed in size!!\n\tEXITING" << endl;
            exit(1);
        }
        if(sz == 0){
            cerr << "sz shouldn't be 0!!\n\tEXITING" << endl;
            exit(1);
        }
#endif
    }
}

static void create_kmer(const Dataset &data,
                        const map<pair<int, char>, double> &pwmHash,
                        const vector<char> &nucleotideStack,
                        const vector<double> &bestCase,
                        map<std::string, double> &retHash,
                        const double cutoff){
    if(data.settings.verbose){
        cout << "Generating kmer list.\n";
    }
    string temp = "";
    for(size_t nucleotide = 0; nucleotide < nucleotideStack.size(); ++nucleotide){
        try{
            // DO NOT USE to_string(args) on char!!!!!!
            temp = nucleotideStack[nucleotide];
            // use of temp to construct string from char of nucleotideStack
            // then pass as key to retHash, as retHash takes strings
            retHash[temp] = pwmHash.at( {1, nucleotideStack[nucleotide] } );
        }
        catch(...){
            cerr << "line 184 Enumerate_kmer.cpp " << nucleotide << ' '
                 << nucleotideStack[nucleotide] << endl;
            exit(1);
        }
    } // end for
    vector<double> maxScores(bestCase.size() + 1);
    try{
        maxScores.at(bestCase.size()) = 0;
        for(int i = static_cast<int>(bestCase.size()) - 1; i >= 0; --i){
            maxScores.at(i) = maxScores.at(i+1) + bestCase.at(i);
        }
    }
    catch(...){
        cerr << "line 207 or 205 out of range error" << endl;
        exit(1);
    }

    double maxscore = 0.0, score = 0.0;
    int length = 0;
    string key = "";

    for(size_t i = 1; i < bestCase.size(); ++i){
        for(auto iter = retHash.begin(); iter != retHash.end(); /*++iter*/){
            key = iter->first;
            score = iter->second;
            length = static_cast<int>(key.size());
            try{
                maxscore = maxScores.at(length) + score;
            }
            catch(...){
                cerr << "out of range line 229 Enumerate_kmer.cpp" << endl
                     << "i: " << i << ' ' << key << ' '
                     << score << endl;
                exit(1);
            }
            try{
                retHash.erase(iter++);
            }
            catch(...){
                cerr << "line 239 Enumerate_kmer.cpp" << endl;
                exit(1);
            }
            if(maxscore < cutoff){

            }
            else{
                for(size_t j = 0; j < nucleotideStack.size(); j++){
                    string newkey = key + nucleotideStack[j];
                    double newscore = score + pwmHash.at({i + 1, nucleotideStack[j]});
                    retHash[newkey] = newscore;
                }
            }
        } // for retHash loop
    } // for bestCase loop
}

static double get_cutoff(const Dataset &data){
    if(data.settings.verbose){
        cout << "\tSearching for pre-calculated cutoff" << '\n';
    }

      string curr_line = "";
      // will need to change below line eventually to be specified
      ifstream IN_HANDLE("src/PWM_SCORES_FINAL.txt");

    if(!IN_HANDLE){
        cerr << "Failure to open src/PWM_SCORES_FINAL.txt\n\tEXITING" << endl;
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
    cerr << "Unable to find pre-calculated cutoff in file\n\tEXITING" << endl;
    exit(1);
    return 0.0;
}

// EFFECTS: finds maximum mapped value
static double findMax(const map<int, double> &v){
    double max = v.begin()->second;
    for(auto i : v){
        if(i.second > max){
            max = i.second;
        }
    }
    return max;
}
