// # This program will use cached computations to save us some compute time
// #prints previously cached items to $cachefile and non-cached to $outfile
// # also filters out kmers based on previous DNase filtering

// DBI library

#include "../iterativeSEM.hpp"
#include <iostream>
#include <cstdlib>
using namespace std;
                              // defined as "none" in iterativeSEM.hpp
  // takes in_file, out_file, out_cache, cache(location of cache), current iteration,
  // all of them are file names, and are specific to nucleotide and position

  // cachefile is the output cache
  // cache is the input cache
void checkCache(Dataset &data, string outfile){
  bool newcache = false;
  if(data.cachefile != ""){
    string cmd = "rm -f " + data.cachefile;
    system(cmd.c_str());
    newcache = true;
  }
  if(outfile != "none"){
    string cmd = "rm -f " + outfile;
    system(cmd.c_str());
  }

  if(data.settings.verbose){
    cout << "Querying cache for processed kmers.\n";
  }

  if(newcache){
    // open infile
    // open outfile
  }
  else{
    if(data.settings.verbose){
      cout << "No existing cache!\n";
    }

  }

  // ASK ABOUT THIS
  // ASK ABOUT THIS
  // ASK ABOUT THIS
  // ASK ABOUT THIS
  // ASK ABOUT THIS

}
