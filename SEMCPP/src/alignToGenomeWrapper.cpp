#include "../iterativeSEM.hpp"
#include <iostream>
using namespace std;

int getLength(Dataset &data);
static void align_SNPs(Dataset &data, int length, const vector<string> &nucleotideStack);
void changeBase(Dataset &data, int position, string nucleotide);
void checkCache(Dataset &data);
void seq_col_to_fa(Dataset &data);
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
  for(int i = 0; i < length; i++){
    for(int j = 0; j < static_cast<int>(nucleotideStack.size()); j++){
      string name = nucleotideStack[j] + "_pos" + to_string(i);
      changeBase(data, i, nucleotideStack[j]);
      checkCache(data);
      seq_col_to_fa(data);
      bowtie_genome_map(data, length);
    }
  }
}

// getlength accesses first(?) element of kmerHash and returns the key length
int getLength(Dataset &data){
  return static_cast<int>(data.kmerHash.begin()->first.size());
}
