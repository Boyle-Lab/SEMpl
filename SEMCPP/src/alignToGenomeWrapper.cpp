#include "iterativeSEM.hpp"
#include <iostream>
using namespace std;

static int getLength(Dataset &data);
static void align_SNPs(Dataset &data, int length, const vector<string> &nucleotideStack);
void changeBase(const Dataset &data, int position, string nucleotide,
                vector<string> &new_kmer_vec);
bool seq_col_to_fa(Dataset &data);
//void bowtie_genome_map(Dataset &data, int length, const string& genome,
//                                                     const string& file);

// default genome is "hg19"
// INFILE FROM ORIGINAL ALGORITHM IS ENUMERATED_KMER
void alignToGenomeWrapper(Dataset &data, int iteration, const string &genome) {

    vector<string> nucleotideStack(4);
    nucleotideStack.push_back("A");
    nucleotideStack.push_back("C");
    nucleotideStack.push_back("T");
    nucleotideStack.push_back("G");

    // step 1: get the length of kmer
    int length = getLength(data);
    if(data.settings.verbose){
        cout << "Aligning" << '\n';
    }

  // step 2: iterate through all the positions and nucleotides,
  //            creating SNP files and aligning to genome
    align_SNPs(data, length, nucleotideStack);
}

// INFILE FROM ORIGINAL ALGORITHM IS ENUMERATED_KMER
static void align_SNPs(Dataset &data, int length,
                       const vector<string> &nucleotideStack){

    vector<string> cache_output(nucleotideStack.size());
    string name = "";

    string CWD =  "../" + data.output_dir + "/ALIGNMENT/";

    system( string("mkdir -p " + CWD).c_str() );

    string genome = "";

    string fa_file = "";

    string bowtie_output = "../";

    vector<string> new_kmer;

    bool non_zero_file_size = false;

    for(int position  = 0; position < length; ++position){
        for(int j = 0; j < static_cast<int>(nucleotideStack.size()); ++j){

            name = nucleotideStack[j] +  "_pos" + to_string(position);

            new_kmer.clear();
                                      // nucleotide
            changeBase(data, position, nucleotideStack[j], new_kmer);

// void checkCache(Dataset &data, vector<string> &in_file, vector<string> &out_cache,
//                 const string &cachefile);
            checkCache(data, new_kmer, cache_output, data.cachefile);
            // pass in a sequence column, which is from output of checkCache
            // cachefile in Dataset is $cache from original algorithm!!!

            fa_file = CWD + name + ".fa";
            bowtie_output = CWD + name + ".bed";

            non_zero_file_size = seq_col_to_fa(cache_output, fa_file);
            if(non_zero_file_size){
                bowtie_genome_map(data, length, "../data/hg19", fa_file, bowtie_output);
            }
            cache_output.clear();
        }
    }
}

// getlength accesses first(?) element of kmerHash and returns the key length
static int getLength(Dataset &data){
    return static_cast<int>(data.kmerHash.begin()->first.size());
}
