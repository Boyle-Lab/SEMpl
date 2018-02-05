#include "src/iterativeSEM.hpp"
#include "src/common.hpp"
#include <iostream>
using namespace std;

static void align_SNPs(Dataset &data, string CWD, vector<string> &new_kmer, int length, int position, char bp);
void find_signal(Dataset &data, int length, int position, char bp);

void alignToGenomeWrapper(Dataset &data, int iteration, const string genome) {

    const vector<char> nucleotideStack{'A', 'C', 'G', 'T'};
    vector<string> new_kmer;

    // step 1: get the length of kmer
    cout << data.kmerHash.size() << endl;
    int length = getLength(data);
    if(data.settings.verbose){
        cout << "\tAligning..." << length << endl << flush;
    }

    // Build ALIGNMENT directory if it doesn't exist
    const string CWD =  "./" + data.output_dir + "ALIGNMENT/";
    if(system( string("mkdir -p " + CWD).c_str() ) != 0){
        cerr << "problem running mkdir -p " << CWD << endl;
        exit(1);
    }

    // step 2: iterate through all the positions and nucleotides,
    //            creating SNP files and aligning to genome

    // clear data.sig_deets*
    data.sig_deets_maximum.clear();
    data.sig_deets_counter.clear();
    data.sig_deets_stdev.clear();
    data.sig_deets_sterr.clear();
    data.accumSummary_data.align_accum_lines.clear();

    // Iterate through each position and generate all possible nucleotide changes
    for(int position  = 0; position < length; ++position){
        for(int j = 0; j < static_cast<int>(nucleotideStack.size()); ++j){
            new_kmer.clear();

            try{
                // creates new_kmer vector from copying over data.kmerHash
                // and changing a single nucleotide
                changeBase(data, position, nucleotideStack[j], new_kmer,
                           genome);
            }
            catch(...){
                cerr << "exception thrown from changeBase" << endl;
                exit(1);
            }

            // Align to the genome
            align_SNPs(data, CWD, new_kmer, length, position, nucleotideStack[j]);

            // Get signal for all alignments
            //  and Write alignments to cache to prevent duplicate processing
            find_signal(data, length, poition, nucleotideStack[j]);
        }
    }

    if(data.settings.verbose){
        cout << "FINISH" << endl;
    }
}

// Align all pseudo-SNP kmers to genome and filter with DNase peaks
//  in one combined function.
//  output is sorted and unique
static void align_SNPs(Dataset &data, string CWD, vector<string> &new_kmer,
                       int length, int position, char bp) {

    string name = "";
    const string genome = "./data/hg19";
    string fa_file = "";
    string bowtie_output = "./";
    int non_zero_file_size;
    vector<string> cache_to_align;

    cache_to_align.clear();

    name = bp + string("_pos") + to_string(position);

    // Filter our kmers by existing kmers in the cache so that we don't
    //   need to re-process them.
    try{
        checkCache(data, new_kmer, cache_to_align, data.cacheDB,
                Dataset::accumSummary_type::accumSummary_dest::alignment,
                position, bp);
        // new_kmer and cache to align would be swapped if there is no cache
        // rather than assigned, much faster
    }
    catch(...){
        cerr << "exception thrown from checkCache" << endl;
        exit(1);
    }

    #ifdef DEBUG
    cerr << position << bp << " to_align: " << cache_to_align.size() << endl;
    #endif

    // align these files to the genome through bowtie_genome_map function
    fa_file = CWD + name + ".fa";
    bowtie_output = CWD + name + ".bed";

    non_zero_file_size = seq_col_to_fa(cache_to_align, fa_file);
    if(non_zero_file_size > 500){
        bowtie_genome_map(length, "./data/hg19", fa_file, bowtie_output,
                          data.DNase_file, data.settings.threads, data.settings.verbose);
    } else if (non_zero_file_size > 0) {
        bowtie_genome_map(length, "./data/hg19", fa_file, bowtie_output,
                          data.DNase_file, 1, data.settings.verbose);
    } else {
        //create the empty file - we use the files to creat the map later on
        // there is probably a better way to do this
        // i thin kthat we can get rid of this now (2/5)
        ofstream OUT(bowtie_output);
        OUT.close();
    }

}

// assumes alignment is filled from previous function
// and that the output is sorted, and contains only unique string values

// NOTE: iteratively constructs the output from
//       findMaximumAverageSignalWrapper(args),
//       as opposed to aggregately, as done in the original algorithm
void find_signal(Dataset &data, int length, int position, char bp){
    if(data.settings.verbose){
        cout << "Finding the average signal" << endl;
    }

    string file;

    // Calculate summary score for all alignments
    try{
        if(data.settings.verbose){
            cout << "\taccumSummary_scale(args) is running..." << flush;
        }
        accumSummary_scale(data, data.bigwig_file, file, length,
                           Dataset::accumSummary_type::accumSummary_dest::alignment);
        if(data.settings.verbose){
            cout << "FINISH" << endl;
        }
    }
    catch(...){
        cerr << "Problem with accumSummary_scale\n\tEXITING" << endl;
        exit(1);
    }

    // Write NEW computed scores to our cache so that we don't need to do previous steps again
    try{
        if(data.settings.verbose){
            cout << "\twriteCache(args) is running..." << flush;
        }
        writeCache(data, data.cacheDB,
                   Dataset::accumSummary_type::accumSummary_dest::alignment);
        if(data.settings.verbose){
            cout << "FINISH" << endl;
        }
    }
    catch(...){
        cerr << "problem with writeCache\n\tEXITING" << endl;
        exit(1);
    }
// do we need past here?

    // SHOULD THERE BE AN ERROR CHECK IF signal_cache_enumerate IS EMPTY????
    try{
        if(data.settings.verbose){
            cout << "\tsorting..." << flush;
        }
        sort(data.signal_cache[ {position, bp} ].begin(),
             data.signal_cache[ {position, bp} ].end() );
        if(data.settings.verbose){
            cout << "FINISH" << endl;
        }
        data.signal_output.resize(data.signal_cache[ {position, bp} ].size()
                                + data.accumSummary_data.align_accum_lines.size());
        // returns iterator to one past the location of the last copy
        if(data.settings.verbose){
            cout << "\tcopying..." << flush;
        }
        auto iter = data.signal_output.begin();
        if(data.accumSummary_data.align_accum_lines.size() != 0 ){
            iter = copy(data.accumSummary_data.align_accum_lines.begin(),
                             data.accumSummary_data.align_accum_lines.end(),
                             data.signal_output.begin() );
        }
        if(data.settings.verbose){
            cout << "FINISH (calculated)" << data.accumSummary_data.align_accum_lines.size() << endl;
        }

        //  FILLS data.signal_enumerate_output
        if(data.settings.verbose){
            cout << "\tunique copying..." << flush;
        }
        auto end_iter3 = unique_copy(data.signal_cache[ {position, bp} ].begin(),
                    data.signal_cache[ {position, bp} ].end(),
                    iter);
        if(data.settings.verbose){
            cout << "FINISH (cache)" << data.signal_cache[ {position, bp} ].size() << endl;
        }
        data.signal_output.resize(end_iter3 - data.signal_output.begin() );
    }
    catch(...){
        cerr << "problem with algorithm usage" << endl;
        exit(1);
    }

    // Now summarize the scores from all alignments (including those from cache)
    // into a single maximum, counter, stdev, and sterr
    try{
        if(data.settings.verbose){
            cout << "\tfinding findMaximumAverageSignal..." << flush;
        }
        findMaximumAverageSignalWrapper(data.signal_output,
                                        data.Signal_data.alignment_maximum,
                                        data.Signal_data.alignment_counter,
                                        data.Signal_data.alignment_stdev,
                                        data.Signal_data.alignment_sterr);
        if(data.settings.verbose){
            cout << "FINISH" << endl;
        }
    }
    catch(...){
        cerr << "problem with findMaximumAverageSignalWrapper(args)" << endl;
        exit(1);
    }
        // iteratively fills sig_deets


    // filled Signal_data data for appropriate type,
    // (alignment, scrambled, or enumerated)

    // first pair type

#ifdef DEBUG
    if(data.settings.verbose){
        cout << "\tpos: " << position << "  bp: " << bp << endl;
    }
#endif
    try{
        auto iter = data.sig_deets_maximum.insert( { {position, bp},
                                            data.Signal_data.alignment_maximum} );

            // cout << '\t' << data.Signal_data.alignment_maximum << "  max inserted" << endl;

        if(!iter.second){
            cerr << "duplicate key inserted into sig_deets_maximum" << endl;
            exit(1);
        }
        // second pair type
        auto iter1 = data.sig_deets_counter.insert( { {position, bp},
                                            data.Signal_data.alignment_counter} );
            // cout << '\t' << data.Signal_data.alignment_counter << "  counter inserted" << endl;
        if(!iter1.second){
            cerr << "duplicate key inserted into sig_deets_counter" << endl;
            exit(1);
        }
        iter = data.sig_deets_stdev.insert( { {position, bp},
                                            data.Signal_data.alignment_stdev} );
            // cout << '\t' << data.Signal_data.alignment_stdev << "  stdev inserted" << endl;
        if(!iter.second){
            cerr << "duplicate key inserted into sig_deets_stdev" << endl;
            exit(1);
        }
        iter = data.sig_deets_sterr.insert( { {position, bp},
                                            data.Signal_data.alignment_sterr} );
        if(!iter.second){
            cerr << "duplicate key inserted into sig_deets_sterr" << endl;
            exit(1);
        }
    }
    catch(...){
        cerr << "exception thrown when inserting alignment data" << endl;
        exit(1);
    }

}
