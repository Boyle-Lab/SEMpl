#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <iomanip>
#include <vector>
extern "C"{
  #include "./lib/libBigWig/bigWig.h"
}

using namespace std;

float findMethyl(string, string, int, float, string, string);

int main(){
  string nstrand = "+";
  string chr = "chr1";
  int loci = 10059;
  float sig = 100.1;
  string nuc = "C";
  string bw_file = "examples/WGBS/ENCFF073DUG.bigWig";
  float meth = findMethyl(nstrand, chr, loci, sig, nuc, bw_file);
  return 0;
}

float findMethyl(string nstrand,string chr,int loci,float sig,string nuc,string bw_file){
  string str_loci = to_string(loci);
  
  if(nstrand == "+"){
    //    loci = loci + 1;
    //string loci_str = to_string(loci);
    //string pos = chr + ":" + loci_str;
    
  }
  if(nstrand == "-"){
    //string loci_str = to_string(loci);
    //string pos = chr + ":" + loci_str;
  }
  int loci_end = loci + 1;
  //  if(sig.compare("NA") !=0 && sig.compare("-256.000000")){
  if(nuc.compare("C") == 0){
    //bwOverlappingIntervals_t *ptr

    //     bigWigFile_t *fp = NULL;
    //fp = bwOpen(bw_file, NULL, "r");
    //if(!fp) {
    //  fprintf(stderr, "An error occured while opening methylation file\n", bw_file);
    //  return 0;
    //}
    
    //int methyl  = bwStats(bw_file, chr, loci, loci_end, 1, dev);
    //cout << methyl << endl;


    //6/6/20
    //    char *fname = new char[bw_file.length() + 1];
    // strcpy(fname, bw_file.c_str());
    //bigWigFile_t *bwFile = bwOpen(fname, NULL, "r");

    //if(bwFile == NULL){
    //	cerr << "Failed to open bw_file: " << bw_file << endl;
    //	exit(1);
    //}

    //6/8/20
    bwOverlappingIntervals_t *ptr = bwGetValues(bwFile, chr, 
				    static_cast<uint32_t>(loci),
                                    static_cast<uint32_t>(loci_end),
				    1);
    
    intervals = bwGetOverlappingIntervals(bwFile, "chr1", loci, loci_end);
    bwDestroyOverlappingIntervals(intervals);
    
  }
  else{
    cout << "bad sig" << sig << endl;
  }
  return 0;
}
