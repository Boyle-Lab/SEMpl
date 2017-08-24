#include "src/iterativeSEM.hpp"
#include "src/common.hpp"
#include <iostream>
using namespace std;

static void align_SNPs(Dataset &data, int length, const vector<char> &nucleotideStack);

// default genome is "hg19"
// INFILE FROM ORIGINAL ALGORITHM IS ENUMERATED_KMER
void alignToGenomeWrapper(Dataset &data, int iteration, const string genome) {

    const vector<char> nucleotideStack{'A', 'C', 'G', 'T'};

    // step 1: get the length of kmer
    int length = getLength(data);
    if(data.settings.verbose){
        cout << "\tAligning\n";
    }

  // step 2: iterate through all the positions and nucleotides,
  //            creating SNP files and aligning to genome
    align_SNPs(data, length, nucleotideStack);
}

// INFILE FROM ORIGINAL ALGORITHM IS ENUMERATED_KMER
static void align_SNPs(Dataset &data, int length,
                       const vector<char> &nucleotideStack){
    string name = "";

    const string CWD =  "./" + data.output_dir + "ALIGNMENT/";
    const string genome = "./data/hg19";

    if(system( string("mkdir -p " + CWD).c_str() ) != 0){
        cerr << "problem running mkdir -p " << CWD << endl;
        exit(1);
    }
    

    string fa_file = "";

    string bowtie_output = "./";

    vector<string> new_kmer;

    bool non_zero_file_size = false;

    vector<string> cache_to_align;

    for(int position  = 0; position < length; ++position){
        for(int j = 0; j < static_cast<int>(nucleotideStack.size()); ++j){

            cache_to_align.clear();
            new_kmer.clear();

            name = nucleotideStack[j] + string("_pos") + to_string(position);
            
                                      // nucleotide
	        
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
                        // CHECK THAT THIS IS CORRECT RELATIVE LOCATION
            #ifdef DEBUG
            // if(position == 4){
            //     cerr << "new_kmer:\n";
            //     for(auto val : new_kmer){
            //         cerr << '#' << val << "#\n";
            //     }
            // }
            #endif

            try{
                checkCache(data, new_kmer, cache_to_align, data.cachefile,
                        Dataset::accumSummary_type::accumSummary_dest::alignment,
                        position, nucleotideStack[j]);
            }
            catch(...){
                cerr << "exception thrown from checkCache" << endl;
                exit(1);
            }

            #ifdef DEBUG
            cerr << position << nucleotideStack[j] << " to_align: " << cache_to_align.size() << endl;
            cerr << '\t' << " amount already aligned: " << new_kmer.size() - cache_to_align.size() << endl;
            cerr << '\t' << " new_kmer.size() " << new_kmer.size() << endl;
            // if(position == 4){
            //     cerr << "cache_to_align:\n";
            //     for(auto val : cache_to_align){
            //         cerr << '#' << val << "#\n";
            //     }
            //     cerr << "aligned:\n";
            //     for(auto val : data.signal_cache[ { position, nucleotideStack[j] } ] ){
            //         cerr << '#' << val << "#\n";
            //     }
            //     cerr << endl;
            // }
            #endif
            // pass in a sequence column, which is from output of checkCache
            // cachefile in Dataset is $cache from original algorithm!!!

            fa_file = CWD + name + ".fa";
            bowtie_output = CWD + name + ".bed";

            non_zero_file_size = seq_col_to_fa(cache_to_align, fa_file);
            if(non_zero_file_size){
                bowtie_genome_map(length, "./data/hg19", fa_file, bowtie_output,
                                  data.settings.verbose);
            }
            #ifdef DEBUG
                
            #endif
        }
    }
}
