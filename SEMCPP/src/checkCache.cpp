// # This program will use cached computations to save us some compute time
// #prints previously cached items to $cachefile and non-cached to $outfile
// # also filters out kmers based on previous DNase filtering

// DBI library

#include "iterativeSEM.hpp"
extern "C"{
  #include "../lib/sqlite3/sqlite3.h"
}
#include <sys/stat.h>
#include <iostream>
using namespace std;

  // takes in_file, out_file, out_cache, cache(location of cache), current iteration,
  // all of them are file names, and are specific to nucleotide and position

  // cachefile is the output cache
  // cache is the input cache
                               // default is "none", defined in iterativeSEM.hpp

// EFFECTS: returns true if file already exists in current directory,
//          false otherwise
static bool fileExists(const string &filename);
static void problemEncountered(const int &message, const string &what);
static bool isRowReady(const int &message);


void checkCache(const Dataset &data, string outfile){
  vector<pair<string, int> > output;

  bool newcache = fileExists(data.cachefile);
    // data_local.cachefile is Cache in original algorithm!!
    // data_local.cachefile was defined in generateSNPEffectMatrix!!

  if(data.settings.verbose){
    cout << "Querying cache for processed kmers.\n";
  }

  sqlite3 *cacheDB;
  int message;
  message = sqlite3_open(data.cachefile.c_str(), &cacheDB);
  problemEncountered(message, "open");

  string msg = "SELECT count(*) FROM seen_cache WHERE kmer=? AND iter!=?";
  //const char* tail_msg;

  if(newcache){
    ofstream OUTF(outfile);
    if(!OUTF){
      cerr << "There is a problem with " << data.cachefile << "\n\tEXITING\n";
      exit(1);
    }

    sqlite3_stmt* seen_query;
    message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &seen_query, NULL);
    //if(*tail_msg) delete[] tail_msg;
    problemEncountered(message, "seen query");

    msg = "SELECT kmer, alignment FROM kmer_cache WHERE kmer=?";
    sqlite3_stmt* data_query;
    message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &data_query, NULL);
    //if(*tail_msg) delete[] tail_msg;
    problemEncountered(message, "data query");

    msg = "INSERT OR IGNORE INTO seen_cache VALUES(?, ?)";
    sqlite3_stmt* staged_query;
    message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &staged_query, NULL);
    //if(*tail_msg) delete[] tail_msg;
    problemEncountered(message, "staged query");
    // range for loop, ranges over every pair within the "map"
    for(auto kmer_pair : data.kmerHash){
      message = sqlite3_bind_text(data_query, 1, kmer_pair.first.c_str(), static_cast<int>(kmer_pair.first.size()), NULL);
      problemEncountered(message, "bind_text for data_query");

      message = sqlite3_step(data_query);
      problemEncountered(message, "step data_query");
      if(!isRowReady(message)){
        cerr << "Row isn't ready!!\n\tEXITING\n";
        exit(1);
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
        exit(1);
      }
      vector<pair<string, int> > data_local;

      message = sqlite3_column_type(data_query, 0);
      if(message != SQLITE3_TEXT){
        cerr << "first column type from data_query was not expected!!\n\tEXITING";
        exit(1);
      }
      if(message != SQLITE_INTEGER){
        cerr << "second column type from data_query was not expected!!\n\tEXITING";
        exit(1);
      }
      //const unsigned char* store = sqlite3_column_text(data_query, i);
      // IMPLEMENTATION DEPENEDENT JUST BELOW
      // sqlite3_column_text returns const unsigned char*, but C++ string library words with const char*
      const char* text = reinterpret_cast<const char*>(sqlite3_column_text(data_query, 0));
      int iter = sqlite3_column_int(data_query, 1);
      data_local.push_back({text, iter});


      if(data_local.size() > 1){
        output.push_back(data_local[1]);
      }
      else{
        message = sqlite3_bind_text(seen_query, 1, kmer_pair.first.c_str(), static_cast<int>(kmer_pair.first.size()), NULL);
        problemEncountered(message, "bind_text for seen_query");
        // int sqlite3_bind_int(sqlite3_stmt*, int, int);
        message = sqlite3_bind_int(seen_query, 2, data.settings.iteration);
        problemEncountered(message, "bind_int for seen_query");


        num_col = sqlite3_column_count(seen_query);
        if(num_col < 0){
          cerr << "num_col for seen_query is less than 1!!\n\tEXITING" << endl;
          exit(1);
        }
        vector<int> count;
        for(int i = 0; i < num_col; i++){
          // IMPLEMENTATION DEPENEDENT JUST BELOW
          //count.push_back(reinterpret_cast<const char*>(sqlite3_column_text(seen_query, i)));

          count.push_back(sqlite3_column_int(seen_query, i));
        }
        if(count[0] > 0){
          // have already seen this before, was filtered
        }
        else{
          message = sqlite3_bind_text(staged_query, 1, kmer_pair.first.c_str(), static_cast<int>(kmer_pair.first.size()), NULL);
          problemEncountered(message, "bind_text for staged_query");
          // int sqlite3_bind_int(sqlite3_stmt*, int, int);
          message = sqlite3_bind_int(staged_query, 2, data.settings.iteration);
          problemEncountered(message, "bind_int for staged_query");

          message = sqlite3_step(staged_query);
          if(message != SQLITE_DONE){
            cerr << "Statement is not done!\n\tEXITING" << endl;
            exit(1);
          }
        }

        message = sqlite3_step(data_query);
        if(message != SQLITE_DONE){
          cerr << "Statement is not done!\n\tEXITING" << endl;
          exit(1);
        }
        OUTF << kmer_pair.first << '\n';
      }

    }
    OUTF.close();
    OUTF.open(data.cachefile);
    for(auto el : output){
      OUTF << el.first << '\n';
    }
    OUTF.close();
  }
  else{
    if(data.settings.verbose){
      cout << "No existing cache!\n";
    }
    // BLOB? why is it a BLOB? I used "text" not "blob", will need to ask about this

      sqlite3_stmt *build_statement;

      msg = "CREATE TABLE kmer_cache (kmer TEXT PRIMARY KEY NOT NULL, alignment BLOB)";
      message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &build_statement, NULL);
      problemEncountered(message, "preparing build_statement");
      message = sqlite3_step(build_statement);
      if(message != SQLITE_DONE){
        cerr << "Build statement is not done!\n\tEXITING\n";
        exit(1);
      }

      msg = "CREATE UNIQUE INDEX kmerIDX ON kmer_cache(kmer)";
      message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &build_statement, NULL);
      problemEncountered(message, "preparing build_statement");
      message = sqlite3_step(build_statement);
      if(message != SQLITE_DONE){
        cerr << "Build statement is not done!\n\tEXITING\n";
        exit(1);
      }

      msg = "CREATE TABLE seen_cache (kmer TEXT PRIMARY KEY NOT NULL, iter INT NOT NULL)";
      message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &build_statement, NULL);
      problemEncountered(message, "preparing build_statement");
      message = sqlite3_step(build_statement);
      if(message != SQLITE_DONE){
        cerr << "Build statement is not done!\n\tEXITING\n";
        exit(1);
      }

      msg = "CREATE UNIQUE INDEX seenIDX ON seen_cache(kmer)";
      message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &build_statement, NULL);
      problemEncountered(message, "preparing build_statement");
      message = sqlite3_step(build_statement);
      if(message != SQLITE_DONE){
        cerr << "Build statement is not done!\n\tEXITING\n";
        exit(1);
      }

      sqlite3_stmt *staged_query;
      msg = "INSERT OR IGNORE INTO kmer_cache VALUES(?,?)";
      message = sqlite3_prepare(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &staged_query, NULL);
      problemEncountered(message, "staged_query for INSERT OR IGNORE");

      for(auto kmer_pair : data.kmerHash){
        message = sqlite3_bind_text(staged_query, 1, kmer_pair.first.c_str(), static_cast<int>(kmer_pair.first.size()), NULL);
        problemEncountered(message, "bind text for inserting into seen_cache");
        message = sqlite3_bind_int(staged_query, 2, data.settings.iteration);
        problemEncountered(message, "bind int for inserting into seen_cache");
        message = sqlite3_step(staged_query);
        if(message != SQLITE_DONE){
          cerr << "Build statement is not done!\n\tEXITING\n";
          exit(1);
        }
      }
      // do I copy over something? investigate
      // ANSWER: copying is not necessary, will approrpiately use kmerHash when it is necessary
  }
  message = sqlite3_close(cacheDB);
  problemEncountered(message, "closing the connection");
}


// taken from internet
// returns true if file exists
static bool fileExists(const string &filename){
  struct stat buffer;
  return (stat (filename.c_str(), &buffer) == 0);
}

static void problemEncountered(const int &message, const string &what){
  if(message != SQLITE_OK){
    cerr << "Problem encountered opening " << what << "!\n\tEXITING\n";
    exit(1);
  }
}

static bool isRowReady(const int &message){
  return message == SQLITE_ROW;
}
