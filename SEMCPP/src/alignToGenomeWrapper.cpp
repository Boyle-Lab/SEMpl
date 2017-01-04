#include "../iterativeSEM.hpp"
#include <iostream>
using namespace std;

static int getLength(Dataset &data);
static void align_SNPs(Dataset &data, int length, const vector<string> &nucleotideStack);
void changeBase(Dataset &data, int position, string nucleotide, vector<string> &new_kmer_vec);
void checkCache(Dataset &data, string outfile);
bool seq_col_to_fa(Dataset &data, int col);
void bowtie_genome_map(Dataset &data, int length);

// default genome is "hg19"
void alignToGenomeWrapper(Dataset &data, int iteration, string genome = "hg19") {

  vector<string> nucleotideStack(4);
  nucleotideStack.push_back("A");
  nucleotideStack.push_back("C");
  nucleotideStack.push_back("T");
  nucleotideStack.push_back("G");

  // step 1: get the length of kmer
  int length = getLength(data);
  if(data.settings.verbose){
    cout << "Aligning" << endl;
  }

  // step 2: iterate through all the positions and nucleotides, creating SNP files and aligning to genome
  align_SNPs(data, length, nucleotideStack);
}

static void align_SNPs(Dataset &data, int length, const vector<string> &nucleotideStack){
  for(int position  = 0; position  < length; position++){
    for(int j = 0; j < static_cast<int>(nucleotideStack.size()); j++){
      string name =  "fa/" + nucleotideStack[j] + "_pos" + to_string(position);

      vector<string> new_kmer;
      bool zero_file_size = false;
                                  // nucleotide
      changeBase(data, position, nucleotideStack[j], new_kmer);
      checkCache(data);
      // pass in a sequence column, which is from output of checkCache
      zero_file_size = seq_col_to_fa(data, 0);
      zero_file_size = !zero_file_size;
      if(!zero_file_size){
        bowtie_genome_map(data, length);
      }
    }
  }
}

// getlength accesses first(?) element of kmerHash and returns the key length
static int getLength(Dataset &data){
  return static_cast<int>(data.kmerHash.begin()->first.size());
}
