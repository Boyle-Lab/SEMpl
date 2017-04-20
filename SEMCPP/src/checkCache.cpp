// # This program will use cached computations to save us some compute time
// #prints previously cached items to $cachefile and non-cached to $outfile
// # also filters out kmers based on previous DNase filtering

// DBI library

#include "./src/iterativeSEM.hpp"
extern "C"{
    #include "./lib/sqlite3/sqlite3.h"
}
#include "src/common.hpp"
#include <iostream>
using namespace std;

  // takes in_file, out_file, out_cache, cache(location of cache), current iteration,
  // all of them are file names, and are specific to nucleotide and position

  // signal_cache in Dataset is the output cache
  // cache is the input cache
                               // default is "none", defined in iterativeSEM.hpp

// EFFECTS: returns true if file already exists in current directory,
//          false otherwise
bool fileExists(const string &filename);
static void problemEncountered(const int message, const string &what);
static void isRowReady(const int message);
// static void prepareStmt(sqlite3 *db, string stmt, sqlite3_stmt *query);
static void checkDone(const int message, const string &s);
static const char* convert_to_const_char(const unsigned char* store);


// MODIFIES: adds appropriate kmers to the specific output_cache
//           see line 121, switch statement
// EFFECTS: checks cache located at cachefile
// NOTE: THE NAMING SCHEME IS CONFUSING
//       out_cache is the argument given to -out_cache in the original algorithm
//       cachefile is the argument given to -cache in the original algorithm
void checkCache(Dataset &data, const vector<string> &in_file, vector<string> &out_cache,
                const string &cachefile, Dataset::accumSummaryData::accumSummary_dest dest){

    bool newcache = fileExists(cachefile);
#ifdef DEBUG
    // cout << "\tDEBUG: " << cachefile;
    // if(newcache) cout << " does not exist";
    // else cout << " exists";
    // cout << '\n';
#endif

    if(data.settings.verbose){
        cout << "Querying cache for processed kmers\n";
    }

    sqlite3 *cacheDB = NULL;
    int message = 0;
    message = sqlite3_open_v2(cachefile.c_str(), &cacheDB, 
                              SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
                              NULL);
    problemEncountered(message, "open");

    string msg = "";
    //const char* tail_msg;

    if(newcache){
        if(data.settings.verbose){
            cout << "Cache does exist\n";
        }

//         int message = sqlite3_prepare_v2(db, stmt.c_str(), static_cast<int>(stmt.size()), &query, NULL);
//     problemEncountered(message, stmt);

        msg = "SELECT count(*) FROM seen_cache WHERE kmer=? AND iter!=?";
        sqlite3_stmt* seen_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &seen_query, NULL);
        problemEncountered(message, msg);
        // prepareStmt(cacheDB, msg, seen_query);

        msg = "SELECT kmer, alignment FROM kmer_cache WHERE kmer=?";
        sqlite3_stmt* data_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &data_query, NULL);
        problemEncountered(message, msg);
        // prepareStmt(cacheDB, msg, data_query);

        msg = "INSERT OR IGNORE INTO seen_cache VALUES(?, ?)";
        sqlite3_stmt* staged_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &staged_query, NULL);
        problemEncountered(message, msg);
        // prepareStmt(cacheDB, msg, staged_query);

        // range for loop, ranges over every pair within the "map"


        for(auto kmer : in_file){
            
                // cout << kmer << endl;


            message = sqlite3_bind_text(data_query, 1, kmer.c_str(),
                      static_cast<int>(kmer.size()), NULL);
            problemEncountered(message, "bind_text for data_query");

            message = sqlite3_step(data_query);
            //problemEncountered(message, "step data_query");
            // isRowReady(message);

            int num_col = sqlite3_column_count(data_query);
#ifdef DEBUG
            cout << "There are " << num_col << " columns in data_query\n";


            if(num_col < 0) {
                cerr << "Number of columns from data_query is less than 0!!\n\tEXITING" << endl;
                exit(1);
            }
            if(num_col != 2) {
                cerr << "Number of columns from data_query is 2!!\n\tEXITING" << endl;
                exit(1);
            }

#endif

            if(num_col > 0){
                // output.push_back(data_local);
                // val is NULL
                const unsigned char* val = sqlite3_column_text(data_query, 1);
                // if(!val) cerr << "val is NULL!!!\n";
                // #ifdef DEBUG
                //     cout << val << endl;
                // #endif
                // const char* text = convert_to_const_char(val);
                // #ifdef DEBUG
                //     cout << text << endl;
                // #endif
                sqlite3_free((char*)val);

#ifdef DEBUG
                cout << "\tDEBUG: text pushed to checkCache output vector\n\t"
                     << "\t   " << text << endl;
#endif

                switch (dest) {
                    case Dataset::accumSummaryData::accumSummary_dest::alignment:
                        data.signal_cache.emplace_back(text);
                    break;
                    case Dataset::accumSummaryData::accumSummary_dest::scrambled:
                        data.signal_cache_scramble.emplace_back(text);
                    break;
                    case Dataset::accumSummaryData::accumSummary_dest::enumerated:
                        data.signal_cache_enumerate.emplace_back(text);
                    break;
                    case Dataset::accumSummaryData::accumSummary_dest::none:
                        cerr << "none shouldn't happen!!" << endl;
                        exit(1);
                    break;
                    default:
                        cerr << "default shouldn't happen!!" << endl;
                        exit(1);
                    break;
                }
                free((char*)text);

            }
            else{
                message = sqlite3_bind_text(seen_query, 1, kmer.c_str(),
                          static_cast<int>(kmer.size()), NULL);
                problemEncountered(message, "bind_text for seen_query");
                // int sqlite3_bind_int(sqlite3_stmt*, int, int);
                message = sqlite3_bind_int(seen_query, 2, data.settings.iteration);
                problemEncountered(message, "bind_int for seen_query");

                num_col = sqlite3_column_count(seen_query);
                if(num_col < 0){
                    cerr << "num_col for seen_query is less than 1!!\n\tEXITING" << '\n';
                    exit(1);
                }
                message = sqlite3_step(seen_query);
                isRowReady(message);
                int count = sqlite3_column_int(seen_query, 0);
                if(count > 0){
                    // have already seen this before, was filtered
                }
                else{
                    message = sqlite3_bind_text(staged_query, 1, kmer.c_str(),
                                static_cast<int>(kmer.size()), NULL);
                    problemEncountered(message, "bind_text for staged_query");
                    // int sqlite3_bind_int(sqlite3_stmt*, int, int);
                    message = sqlite3_bind_int(staged_query, 2, data.settings.iteration);
                    problemEncountered(message, "bind_int for staged_query");

                    // return by reference
                    out_cache.push_back(kmer);

                    message = sqlite3_step(staged_query);
                    checkDone(message, "staged query execution line 138");
                    sqlite3_reset(staged_query);
                    sqlite3_clear_bindings(staged_query);
                }

                message = sqlite3_step(data_query);
                if(message != SQLITE_DONE){
                    cerr << "Statement is not done!\n\tEXITING" << '\n';
                    exit(1);
                }
                //OUTF << kmer << '\n';
            }


            sqlite3_reset(seen_query);
            sqlite3_clear_bindings(seen_query);
            sqlite3_reset(data_query);
            sqlite3_clear_bindings(data_query);
        }

        // signal_cache is filled above on line 117
    }
    else{
        if(data.settings.verbose){
            cout << "No existing cache!\n";
        }
    // BLOB? why is it a BLOB? I used "text" not "blob", will need to ask about this

        // sqlite3_stmt *build_statement = NULL;
        char *z_err_msg = NULL;
        //cout << "Creating Table kmer_cache" << endl;
        msg = "CREATE TABLE kmer_cache (kmer TEXT PRIMARY KEY NOT NULL, alignment BLOB)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        if(message != SQLITE_OK){
            cerr << "\tproblem with exec on create TABLE kmer_cache\n";
            exit(1);
        }
        sqlite3_free(z_err_msg);
       
	//cout << "Creating Unique Index kmerIDX" << endl;
        msg = "CREATE UNIQUE INDEX kmerIDX ON kmer_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        if(message != SQLITE_OK){
                cerr << "\tproblem with exec on create unique index on kmer_cache\n";
                exit(1);
        }
        sqlite3_free(z_err_msg);
        // problemEncountered(message, z_err_msg);
        // prepareStmt(cacheDB, msg, build_statement);
        // message = sqlite3_step(build_statement);
        // checkDone(message, "build statement create unique index kmer_cache");
        // sqlite3_finalize(build_statement);
	//cout << "Creating Table seen_cache" << endl;
        msg = "CREATE TABLE seen_cache (kmer TEXT PRIMARY KEY NOT NULL, iter INT NOT NULL)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        if(message != SQLITE_OK){
            cerr << "\tproblem with create TABLE seen_cache\n";
            exit(1);
        }
        sqlite3_free(z_err_msg);
        // prepareStmt(cacheDB, msg, build_statement);
        // message = sqlite3_step(build_statement);
        // checkDone(message, "build statement create table seen_cache");
        // sqlite3_finalize(build_statement);
	//cout << "Creating Unique Index seenIDX" << endl;
        msg = "CREATE UNIQUE INDEX seenIDX ON seen_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        if(message != SQLITE_OK){
                cerr << "\tproblem with create unique index on seen_cache\n";
                exit(1);
        }
        sqlite3_free(z_err_msg);

        sqlite3_stmt *staged_query = NULL;
        msg = "INSERT OR IGNORE INTO seen_cache VALUES(?,?)";
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), 
                                     static_cast<int>(msg.size()), 
                                     &staged_query, NULL);
        problemEncountered(message, msg);

        if(staged_query == NULL){
            cerr << "staged_query shouldn't be NULL\n";
            exit(1);
        }
	//cout << "Binding to cache" << endl;
        for(auto kmer : in_file){
            message = sqlite3_bind_text(staged_query, 1, kmer.c_str(),
                      static_cast<int>(kmer.size()), NULL);
            problemEncountered(message, "bind text for inserting into seen_cache");
            message = sqlite3_bind_int(staged_query, 2, data.settings.iteration);
            problemEncountered(message, "bind int for inserting into seen_cache");
	    out_cache.push_back(kmer);
            message = sqlite3_step(staged_query);
            if(message != SQLITE_DONE){
                cerr << "Build statement is not done!\n\tEXITING\n";
                exit(1);
            }
            sqlite3_reset(staged_query);
            sqlite3_clear_bindings(staged_query);
        }
	//cout << "Finalizing staged_query" << endl;
        message = sqlite3_finalize(staged_query);
        if(message != SQLITE_OK){
            cerr << "problem with finalize of staged_query" << endl;
            cerr << "Error code: " << message << endl;
        }
      // do I copy over something? investigate
      // ANSWER: copying is not necessary, will appropriately use kmerHash when it is necessary
    }


    //cout << "Closing cacheDB" << endl;
    message = sqlite3_close_v2(cacheDB);
    problemEncountered(message, "closing the connection");
}


// static void prepareStmt(sqlite3 *db, string stmt, sqlite3_stmt *query){
//     int message = sqlite3_prepare_v2(db, stmt.c_str(), static_cast<int>(stmt.size()), &query, NULL);
//     problemEncountered(message, stmt);
// }

static void problemEncountered(const int message, const string &what){
    if(message != SQLITE_OK){
        cerr << "Problem encountered with " << what << " !\n\tEXITING\n";
        cout << "\tError code: " << message << endl;
        exit(1);
    }
}

static void checkDone(const int message, const string &s){
    if(message != SQLITE_DONE){
        cerr << s << " is not done!\n\tEXITING\n";
        exit(1);
    }
}

static void isRowReady(const int message){
    if(message != SQLITE_ROW){
        cerr << message << endl;
        cerr << "Row isn't ready!!\n\tEXITING" << endl;
        exit(1);
    }
}

static const char* convert_to_const_char(const unsigned char* store){
    return reinterpret_cast<const char*>(store);
}
