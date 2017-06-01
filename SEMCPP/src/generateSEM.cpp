#include "./src/iterativeSEM.hpp"
#include <cmath>
#include "./src/common.hpp"
#include <vector>
#include <cstring>
#include <map>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>

using namespace std;


/*  Effect: To generate a SNP effect matrix using maximum signal
*           signals of the generated SNPS and baseline.
*
*/

// TF_name and output_dir are contained in data!!!!
void generateSEM(const Dataset &data){

    // SNPSignal = data.output_dir + "/SIGNAL/signal.maximums"
    // output from running findMaximumAverageSignalWrapper on alignment
    // BaseSignal = data.output_dir + "/BASELINE/baseline.maximums"
    // output from running findMaximumAverageSignalWrapper on enumerated

#ifdef DEBUG
    // if(data.accumSummary_data.align_accum_lines.empty()){
    //     cout << "\tdata.accumSummary_data.align_accum_lines shouldn't be\n"
    //          << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    //     exit(1);
    // }

    // if(data.accumSummary_data.enum_accum_lines.empty()){
    //     cout << "\tdata.accumSummary_data.enum_accum_lines shouldn't be\n"
    //          << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    //     exit(1);
    // }

    // if(data.accumSummary_data.enum_accum_max.empty()){
    //     cout << "\tdata.accumSummary_data.enum_accum_max shouldn't be\n"
    //          << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    //     exit(1);
    // }

    // if(data.accumSummary_data.align_accum_max.empty()){
    //     cout << "\tdata.accumSummary_data.align_accum_max shouldn't be\n"
    //          << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    //     exit(1);
    // }

#endif

// I don't think it's necessary to check that Enumerated is in $line
// as I distinctly separated the data from calling accumSummary_scale on
// the enumerated and alignment data separately.

// in findMaximumAverageSignalWrapper, the file grabs the output from
// accumSummary_scale when ran on enumerated, as well as scrambled,
// if they are present. Two cases, enumerated data, or scrambled and enumerated
// data.


// print OUT_HANDLE "$file\t$maximum\t$count\t$stdev\t$sterr\n";
// use to translate indexes in this file's Perl equivalent

    double enum_ = data.Signal_data.enumerate_maximum;
    //double enum_err = data.Signal_data.enumerate_sterr;

    double score = data.Signal_data.alignment_maximum;
    double st_err = data.Signal_data.alignment_sterr;

    double maximum = 0.0;

    map<pair<int, string>, int> SNPEffect;
    map<pair<int, string>, int> STDErr;

    // will need to address the splitting of @fields, which would be the first
    // thing printed on a given line in the output from findMaximumAverageSignalWrapper.

    // I believe the first first line of output would 
    // normally be the line from accum_summary
    // which is held in the struct and should be 
    // utilized as seen below.


    vector<string> fields1;
    const string splitBy =  "_";
    for(size_t i = 0; i < data.accumSummary_data.align_accum_lines.size(); ++i){
        split_string(data.accumSummary_data.align_accum_lines[i], splitBy, fields1);
        string BP = fields1[0];
        string loc = fields1[1];

        while ( loc.find("pos") != loc.npos){
            loc.replace(loc.find("pos"), 3, "");
        }
        if (strtol(loc.c_str(), 0, 10) > maximum ){
            maximum = strtol(loc.c_str(), 0, 10);
        }

        // Hash with a pair key of locations and base pair
        SNPEffect [ make_pair(strtol(loc.c_str(), 0, 10), BP)]= log(score/enum_)/log(2.0);
        STDErr [make_pair(strtol(loc.c_str(), 0, 10), BP)]= st_err/enum_;
        fields1.clear();
    }

    //Constructing Output File
    stringstream output1;
    stringstream output2;
    output1 << data.output_dir << "/" << data.TF_name << ".sem";
    output2 << data.output_dir << "/" << data.TF_name << ".sterr";
    ofstream sem_file;
    ofstream sterr_file;
    sem_file.open(output1.str());
    sterr_file.open(output2.str());

    //Making output files for the .sem and .sterr
    sem_file << data.TF_name << "\tA\tC\tG\tT\n";
    sterr_file << data.TF_name << "\tA\tC\tG\tT\n";
    for (int j = 0; j <= maximum; j++){
        int p = j+1;
        sem_file << p << "\t" << SNPEffect[make_pair(j, "A")]
                 << "\t" << SNPEffect[make_pair(j, "C")]
                 << "\t" << SNPEffect[make_pair(j, "G")]
                 << "\t" << SNPEffect[make_pair(j, "T")]
                 << "\n";
        sterr_file << p << "\t" << STDErr[make_pair(j, "A")]
                   << "\t" << STDErr[make_pair(j, "C")]
                   << "\t" << STDErr[make_pair(j, "G")]
                   << "\t" << STDErr[make_pair(j, "T")]
                   << "\n";;
    }

    sem_file.close();
    sterr_file.close();

}
