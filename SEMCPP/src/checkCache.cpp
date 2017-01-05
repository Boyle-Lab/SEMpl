// # This program will use cached computations to save us some compute time
// #prints previously cached items to $cachefile and non-cached to $outfile
// # also filters out kmers based on previous DNase filtering

// DBI library

#include "../iterativeSEM.hpp"
extern "C"{
  #include "../lib/sqlite3/sqlite3.h"
}
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>
using namespace std;

  // takes in_file, out_file, out_cache, cache(location of cache), current iteration,
  // all of them are file names, and are specific to nucleotide and position

  // cachefile is the output cache
  // cache is the input cache
                               // default is "none", defined in iterativeSEM.hpp

// EFFECTS: returns true if file already exists in current directory,
//          false otherwise
static bool fileExists(const string &filename);
static void problemEncountered(int message, const string &what);

void checkCache(Dataset &data, string outfile){
  vector<string> output;
  bool newcache = fileExists(data.cachefile);
    // data.cachefile is Cache in original algorithm!!
    // data.cachefile was defined in generateSNPEffectMatrix!!

  if(data.settings.verbose){
    cout << "Querying cache for processed kmers.\n";
  }

  sqlite3 *cacheDB;
  int message;
  string msg = "SELECT count(*) FROM seen_cache WHERE kmer=? AND iter!=?";
  const char* tail_msg;
  
  message = sqlite3_open(data.cachefile.c_str(), &cacheDB);
  problemEncountered(message, "open");

  sqlite3_stmt* seen_query;
  message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &seen_query, &tail_msg);
  delete[] tail_msg;
  problemEncountered(message, "seen query");

  msg = "SELECT kmer, alignment FROM kmer_cache WHERE kmer=?";
  sqlite3_stmt* data_query;
  message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &data_query, &tail_msg);
  delete[] tail_msg;
  problemEncountered(message, "data query");

  msg = "INSERT OR IGNORE INTO seen_cache VALUES(?, ?)";
  sqlite3_stmt* staged_query;
  message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &staged_query, &tail_msg);
  delete[] tail_msg;
  problemEncountered(message, "staged query");


  

  if(outfile != "none"){
    string cmd = "rm -f " + outfile;
    system(cmd.c_str());
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
}


// taken from internet
// returns true if file exists
static bool fileExists(const string &filename){
  struct stat buffer;
  return (stat (filename.c_str(), &buffer) == 0);
}

static void problemEncountered(int message, const string &what){
  if(message != SQLITE_OK){  
    cerr << "Problem encountered opening " << what << "!\n\tEXITING\n";
    exit(EXIT_FAILURE);
  }
}
