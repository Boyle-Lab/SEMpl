#include "src/iterativeSEM.hpp"
#include "src/common.hpp"
#include <iostream>
using namespace std;

static void align_SNPs(Dataset &data, string CWD, vector<string> &new_kmer, int length, int position, char bp);

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

    // lets loop here instead of in the alignment step.
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

            align_SNPs(data, CWD, new_kmer, length, position, nucleotideStack[j]);
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
