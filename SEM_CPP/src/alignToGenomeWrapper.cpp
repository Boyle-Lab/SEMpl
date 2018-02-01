#include "src/iterativeSEM.hpp"
#include "src/common.hpp"
#include <iostream>
using namespace std;

static void align_SNPs(Dataset &data, int length, const vector<char> &nucleotideStack);

void alignToGenomeWrapper(Dataset &data, int iteration, const string genome) {

    const vector<char> nucleotideStack{'A', 'C', 'G', 'T'};

    // step 1: get the length of kmer
    cout << data.kmerHash.size() << endl;
    int length = getLength(data);
    if(data.settings.verbose){
        cout << "\tAligning..." << length << endl << flush;
    }

    // step 2: iterate through all the positions and nucleotides,
    //            creating SNP files and aligning to genome
    align_SNPs(data, length, nucleotideStack);

    if(data.settings.verbose){
        cout << "FINISH" << endl;
    }
}

// Align all pseudo-SNP kmers to genome and filter with DNase peaks
//  in one combined function.
//  output is sorted and unique
static void align_SNPs(Dataset &data, int length,
                       const vector<char> &nucleotideStack) {

    string name = "";
    const string CWD =  "./" + data.output_dir + "ALIGNMENT/";
    const string genome = "./data/hg19";
    string fa_file = "";
    string bowtie_output = "./";
    vector<string> new_kmer;
    bool non_zero_file_size = false;
    vector<string> cache_to_align;

    // Build ALIGNMENT directory if it doesn't exist
    if(system( string("mkdir -p " + CWD).c_str() ) != 0){
        cerr << "problem running mkdir -p " << CWD << endl;
        exit(1);
    }

    // Iterate through each position and generate all possible nucleotide changes
    for(int position  = 0; position < length; ++position){
        for(int j = 0; j < static_cast<int>(nucleotideStack.size()); ++j){

            cache_to_align.clear();
            new_kmer.clear();

            name = nucleotideStack[j] + string("_pos") + to_string(position);

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

            // Filter our kmers by existing kmers in the cache so that we don't
            //   need to re-process them.
            try{
                checkCache(data, new_kmer, cache_to_align, data.cacheDB,
                        Dataset::accumSummary_type::accumSummary_dest::alignment,
                        position, nucleotideStack[j]);
                // new_kmer and cache to align would be swapped if there is no cache
                // rather than assigned, much faster
            }
            catch(...){
                cerr << "exception thrown from checkCache" << endl;
                exit(1);
            }

            #ifdef DEBUG
            cerr << position << nucleotideStack[j] << " to_align: " << cache_to_align.size() << endl;
            // cerr << '\t' << " amount already aligned: " << new_kmer.size() - cache_to_align.size() << endl;
            // cerr << '\t' << " new_kmer.size() " << new_kmer.size() << endl;
//            if(position == 0 && nucleotideStack[j] == 'A'){
//                cerr << "cache_to_align:\n";
//                for(auto val : cache_to_align){
//                    cerr << '#' << val << "#\n";
//                }
//                cerr << "aligned:\n";
//                for(auto val : data.signal_cache[ { position, nucleotideStack[j] } ] ){
//                    cerr << '#' << val << "#\n";
//                }
//                cerr << endl;
                // ofstream out("A_pos1.txt");
                // out << "cache_to_align:\n";
                // for(auto val : cache_to_align){
                //     out << val << "\n";
                // }
                // out << "aligned:\n";
                // for(auto val : data.signal_cache[ { position, nucleotideStack[j] } ] ){
                //     cerr << val << "\n";
                // }

//            }
//            cout << "exit " << __LINE__ << endl;
//            exit(1);
            #endif
            // pass in a sequence column, which is from output of checkCache


            // align these files to the genome through bowtie_genome_map function
            fa_file = CWD + name + ".fa";
            bowtie_output = CWD + name + ".bed";

            non_zero_file_size = seq_col_to_fa(cache_to_align, fa_file);
            if(non_zero_file_size){
                bowtie_genome_map(length, "./data/hg19", fa_file, bowtie_output,
                                  data.DNase_file, data.settings.verbose);
            } else {
                //create the empty file - we use the files to creat the map later on
                // there is probably a better way to do this
                ofstream OUT(bowtie_output);
                OUT.close();
            }
            #ifdef DEBUG

            #endif
        }
    }
}
