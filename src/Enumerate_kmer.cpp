#include "iterativeSEM.hpp"
#include <map>
#include <cmath>
#include <utility>
#include <sstream>
#include <iostream>
#include <exception>
#include <stdlib.h>
#include <algorithm>
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
static void create_kmer(const Dataset &data,
                        const map<pair<int, char>, double> &pwmHash,
                        const vector<char> &nucleotideStack,
                        const vector<double> &bestCase,
                        map<std::string, double> &retHash,
                        const double cutoff);
static void parse_pwm(const Dataset &data,
                      map<pair<int, char>, double> &pwmHash,
                      const vector<char> &nucleotideStack,
                      vector<double> &bestCase);

// REQUIRES: within data, PWM_data is filled in
// EFFECTS: created enumerated kmers using a cutoff and a PWM matrix, returns the map
//          the result is data.kmerHash
void Enumerate_kmer(Dataset &data){
    map< pair<int, char>, double> pwmHash;
    //                      default nucleotide ordering,
    //                      will check in parse_pwm
    //                      for nucleotide ordering (columns) within pwm
    vector<char> nucleotideStack {'A', 'C', 'G', 'T'};
    vector<double> bestCase;
    bool cutoffValid = 0;
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
        //this is no longer reachable
        if(data.settings.verbose){
            cout << "\tNo cutoff defined, so searching for pre-calculated cutoff." << endl;
        }
        cutoff = 0.0;
        #ifdef DEBUG
            cout << "\tusing cutoff: " << cutoff << endl;
        #endif
    }
    else{
        if(data.settings.verbose){
            cout << "\tUsing user defined cutoff." << endl;
        }
        cutoff = data.settings.threshold;
        #ifdef DEBUG
            cout << "\t\tusing cutoff: " << cutoff << endl;
        #endif
    }

    // here we can adjust the cutoff for count thresholds as well
    // now enforcing min and max Kmer counts
    while(!cutoffValid) {
        try{
            create_kmer(data, pwmHash, nucleotideStack, bestCase, data.kmerHash, cutoff);
        }
        catch(...){
            cerr << "Problem with create_kmer!\n\tEXITING" << endl;
            exit(1);
        }

        // trim here for last base added to kmer
        for(auto iter = data.kmerHash.begin(); iter != data.kmerHash.end(); ){
            if(iter->second <= cutoff){
                data.kmerHash.erase(iter++);
            }
            else{
                ++iter;
            }
        }
        if(data.kmerHash.size() > data.settings.maxKmers) {
            // above max we can just trim the kmers to meet our threshold
            cout << "\tKmer count is " << data.kmerHash.size() << ", which is above maximum threshold. Trimming kmer list." << endl;

            std::vector<double> kmer_scores;
            for(auto val : data.kmerHash){
                kmer_scores.push_back(val.second);
            }
            std::sort(kmer_scores.begin(), kmer_scores.end());

            cutoff = *(kmer_scores.end()-data.settings.maxKmers-1);
            cout << "\tEffective cutoff: " << cutoff << endl;
            for(auto iter = data.kmerHash.begin(); iter != data.kmerHash.end(); ){
                if(iter->second <= cutoff){
                    data.kmerHash.erase(iter++);
                }
                else{
                    ++iter;
                }
            }

            cutoffValid = 1;
        }
        else if(data.kmerHash.size() < data.settings.minKmers) {
            // below min, we need to enumerate more
            cutoff = cutoff-1;
            cout << "\tKmer count is " << data.kmerHash.size() << ", which is below minimum threshold. Expanding kmer list." << endl;
            cout << "\tEffective cutoff: " << cutoff << endl;
            // clear data.kmerHash
            data.kmerHash.clear();
        }
        else {
            cerr << "Kmer count is " << data.kmerHash.size() << endl;
            cutoffValid = 1; // All good
        }
    }
    // data.kmerHash is now ready for use


#ifdef DEBUG
    ofstream OUT(data.output_dir + "Enumerated_kmers.txt");
    for(auto val : data.kmerHash){
        OUT << val.first << '\t' << val.second << endl;
    }
    OUT.close();

    data.size_of_kmerHash = data.kmerHash.size();
    // cout << "exit " << __LINE__ << endl;
    // exit(1);
    // cout << data.size_of_kmerHash << endl;
#endif
}

static void parse_pwm(const Dataset &data,
                      map<pair<int, char>, double> &pwmHash,
                      const vector<char> &nucleotideStack,
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

    map<int, double> modifiedFields;

    for(int row = 0; row < (int)data.PWM_data.matrix_arr[0].size(); ++row){
        // cout << row << endl;
        sum = 0;
        modifiedFields.clear();

        for(int column = 0; column < (int)data.PWM_data.matrix_arr.size(); ++column){
            // sums each row
            sum += data.PWM_data.matrix_arr[column][row];
        }

        for(int column = 0; column < (int)data.PWM_data.matrix_arr.size(); ++column){
            modifiedFields[column] = log2((static_cast<double>(data.PWM_data.matrix_arr[column][row]) + 0.25)
                                     / (static_cast<double>(sum) + 1.0)) - log2(0.25);
        }
#ifdef DEBUG
        size_t sz = modifiedFields.size();
#endif

        bestCase.push_back(findMax(modifiedFields));
        try{
            for(int column = 0; column < (int)data.PWM_data.matrix_arr.size(); ++column){
                pwmHash[ { row + 1, nucleotideStack[column] } ] = modifiedFields.at(column);
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
            // then pass as key to retHash, as retHash takes strings for keys
            retHash[temp] = pwmHash.at( {1, nucleotideStack[nucleotide] } );
        }
        catch(...){
            cerr << "line " << __LINE__ << "Enumerate_kmer.cpp " << nucleotide << ' '
                 << nucleotideStack[nucleotide] << endl;
            exit(1);
        }
    } // end for
    #ifdef DEBUG
    // verified that retHash and kmerHash from original algorithm
    // has the same contents at this step in the program
        // for(auto i : retHash){
        //     cerr << "rethash: " << i.first << ' ' << i.second << endl;
        // }
        // exit(1);
    #endif
    vector<double> maxScores(bestCase.size() + 1);
    // will contain data necessary to prune searchpaths (kmer permutations)
    try{
        maxScores.at(bestCase.size()) = 0;
        // if kmer is of full length (number of rows in pwm)
        // then maximum score that can be added is 0
        for(int i = static_cast<int>(bestCase.size()) - 1; i >= 0; --i){
            // constructs maximum possible score at each position of any kmer
            maxScores.at(i) = maxScores.at(i+1) + bestCase.at(i);
            #ifdef DEBUG
            // cout << "maxscores[" << i << "] " << maxScores[i] << endl;
            #endif
        }
    }
    catch(...){
        cerr << "line " << __LINE__ << " out of range error" << endl;
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
                cerr << "out of range line " << __LINE__ << " Enumerate_kmer.cpp" << endl
                     << "i: " << i << ' ' << key << ' '
                     << score << endl;
                exit(1);
            }
            try{
                retHash.erase(iter++);
            }
            catch(...){
                cerr << "line " << __LINE__ << "Enumerate_kmer.cpp" << endl;
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

// EFFECTS: finds maximum mapped value
static double findMax(const map<int, double> &v){
    double max = v.begin()->second;
    for(const auto &i : v){
        if(i.second > max){
            max = i.second;
        }
    }
    return max;
}
