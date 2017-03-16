#include "iterativeSEM.hpp"
#include <cmath>
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
    if(data.accumSummary_data.align_accum_lines.empty()){
        cout << "\tdata.accumSummary_data.align_accum_lines shouldn't be\n"
             << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        exit(1);
    }

    if(data.accumSummary_data.enum_accum_lines.empty()){
        cout << "\tdata.accumSummary_data.enum_accum_lines shouldn't be\n"
             << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        exit(1);
    }

    if(data.accumSummary_data.enum_accum_max.empty()){
        cout << "\tdata.accumSummary_data.enum_accum_max shouldn't be\n"
             << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        exit(1);
    }

    if(data.accumSummary_data.align_accum_max.empty()){
        cout << "\tdata.accumSummary_data.align_accum_max shouldn't be\n"
             << "\tempty!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        exit(1);
    }

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

    double enum_ = data.accumSummary_data.enumerate_maximum;
    double enum_err = data.accumSummary_data.sterr;

    double score = data.accumSummary_data.alignment_maximum;
    double st_err = data.accumSummary_data.alignment_sterr;

    // will need to address the splitting of @fields, which would be the first
    // thing printed on a given line in the output from findMaximumAverageSignalWrapper.

    


}
