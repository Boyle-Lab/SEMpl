// # This program will use cached computations to save us some compute time
// #prints previously cached items to $cachefile and non-cached to $outfile
// # also filters out kmers based on previous DNase filtering

// DBI library

#include "iterativeSEM.hpp"
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
static bool isRowReady(int message);


void checkCache(Dataset &data, string outfile){
  vector<string> output;
  bool newcache = fileExists(data.cachefile);
    // data_local.cachefile is Cache in original algorithm!!
    // data_local.cachefile was defined in generateSNPEffectMatrix!!

  if(data.settings.verbose){
    cout << "Querying cache for processed kmers.\n";
  }

  sqlite3 *cacheDB;
  int message;
  string msg = "SELECT count(*) FROM seen_cache WHERE kmer=? AND iter!=?";
  const char* tail_msg;

  message = sqlite3_open(data.cachefile.c_str(), &cacheDB);
  problemEncountered(message, "open");

  if(newcache){
    sqlite3_stmt* seen_query;
    message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &seen_query, &tail_msg);
    if(*tail_msg) delete[] tail_msg;
    problemEncountered(message, "seen query");

    msg = "SELECT kmer, alignment FROM kmer_cache WHERE kmer=?";
    sqlite3_stmt* data_query;
    message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &data_query, &tail_msg);
    if(*tail_msg) delete[] tail_msg;
    problemEncountered(message, "data query");

    msg = "INSERT OR IGNORE INTO seen_cache VALUES(?, ?)";
    sqlite3_stmt* staged_query;
    message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &staged_query, &tail_msg);
    if(*tail_msg) delete[] tail_msg;
    problemEncountered(message, "staged query");
    // range for loop, ranges over every pair within the "map"
    for(auto kmer_pair : data.kmerHash){
      message = sqlite3_bind_text(data_query, 1, kmer_pair.first.c_str(), static_cast<int>(kmer_pair.first.size()), NULL);
      problemEncountered(message, "bind_text for data_query");

      message = sqlite3_step(data_query);
      problemEncountered(message, "step data_query");
      if(!isRowReady(message)){
        cerr << "Row isn't ready!!\n\tEXITING\n";
        exit(EXIT_FAILURE);
      }

      // seems like I will have to extract data with multiple calls, as it appears
      // that fetchrow_array is not a function within sqlite3.h, or has an equivalent
      // thus, I will do that myself

      // will need to discuss exactly what checkCache.pl does within the DBI interface

      int num_col = sqlite3_column_count(data_query);
      #ifdef DEBUG
        cout << "There are " << num_col << " columns in data_query" << endl;
      #endif

      if(num_col < 0) {
        cerr << "Number of columns from data_query is less than 0!!\n\tEXITING" << endl;
        exit(EXIT_FAILURE);
      }
      vector<string> data_local;
      // const unsigned char *sqlite3_column_text(sqlite3_stmt*, int iCol);
      for(int i = 0; i < num_col; i++){

        const unsigned char* store = sqlite3_column_text(data_query, i);
        // conversion needs to be made so string can receive, will investigate
        //data_local.push_back(store);
      }
      if(data_local.size() > 1){
        output.push_back(data_local[1]);
      }
      else{
        message = sqlite3_bind_text(seen_query, 1, kmer_pair.first.c_str(), static_cast<int>(kmer_pair.first.size()),NULL);
        problemEncountered(message, "bind_text for seen_query");
        // int sqlite3_bind_int(sqlite3_stmt*, int, int);
        message = sqlite3_bind_int(seen_query, 2, data.settings.iteration);
        problemEncountered(message, "bind_int for seen_query");

        //UNSURE OF BELOW line
        message = sqlite3_step(data_query);
        if(message != SQLITE_DONE){
          cerr << "Statement is not done!\n\tEXITING" << endl;
          exit(1);
        }
        ofstream OUTF(data.cachefile);
        OUTF << kmer_pair.first << '\n';
        OUTF.close();
      }

    }
  }



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

static bool isRowReady(int message){
  return message == SQLITE_ROW;
}
