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

// EFFECTS: returns true if file already exists in current directory,
//          false otherwise
bool fileExists(const string &filename);
static void problemEncountered(const int message, const string &what);
static void isRowReady(const int message);
// static void prepareStmt(sqlite3 *db, string stmt, sqlite3_stmt *query);
static void checkDone(const int message, const string &s);
// static const char* convert_to_const_char(const unsigned char* store);


// MODIFIES: adds appropriate kmers to the specific output_cache
//           see line 121, switch statement
// EFFECTS: checks cache located at cachefile
// NOTE: THE NAMING SCHEME IS CONFUSING
//       out_cache is the argument given to -out_cache in the original algorithm
//       cachefile is the argument given to -cache in the original algorithm
void checkCache(Dataset &data, const vector<string> &in_file, vector<string> &out_cache,
                const string &cachefile, Dataset::accumSummary_type::accumSummary_dest dest){

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
    // if(message != SQLITE_OK){
    //     cerr << '\t' << '\t' << sqlite3_errmsg(cacheDB) << endl;
    //     ifstream test(cachefile);
    //     cerr << '\t';
    //     if(test){
    //         cerr << "success";
    //     }
    //     else{
    //         cerr << "failure";
    //     }
    //     cerr << endl << cachefile << endl;
    //     cerr << endl;
    // }
    problemEncountered(message, "open");

    string msg = "";

    if(newcache){
        if(data.settings.verbose){
            // cout << "Cache does exist\n";
        }

        msg = "SELECT count(*) FROM seen_cache WHERE kmer=? AND iter!=?";
        sqlite3_stmt* seen_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &seen_query, NULL);
        problemEncountered(message, msg);

        msg = "SELECT kmer, alignment FROM kmer_cache WHERE kmer=?";
        sqlite3_stmt* data_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &data_query, NULL);
        problemEncountered(message, msg);

        msg = "INSERT OR IGNORE INTO seen_cache VALUES(?, ?)";
        sqlite3_stmt* staged_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), static_cast<int>(msg.size()), &staged_query, NULL);
        problemEncountered(message, msg);

        for(string kmer : in_file){
            // cout << kmer.c_str() << endl;
            message = sqlite3_bind_text(data_query, 1, kmer.c_str(),
                      -1, NULL);
            problemEncountered(message, "bind_text for data_query");

            message = sqlite3_step(data_query);

            int num_col = sqlite3_column_count(data_query);
#ifdef DEBUG
            // cout << "There are " << num_col << " columns in data_query\n";
            if(num_col < 0) {
                cerr << "Number of columns from data_query is less than 0!!\n\tEXITING" << endl;
                exit(1);
            }
            if(num_col != 2) {
                cerr << "Number of columns from data_query is not 2!!\n\tEXITING" << endl;
                exit(1);
            }
#endif
            // if(num_col > 0) {
            // CHECKING THE NUMBER OF COLS WORKS ONLY FOR THE PERL API
            // I BELIEVE THIS API RETURNS A NULL IF SOMETHING DOESN'T EXIST
            // if(message == SQLITE_ROW){

            const char* text = (char*)sqlite3_column_text(data_query, 1);


            if(text){
                // found (?)
                // result is not NULL

                switch (dest) {
                    case Dataset::accumSummary_type::accumSummary_dest::alignment:
                        data.signal_cache.emplace_back(text);
                    break;
                    case Dataset::accumSummary_type::accumSummary_dest::scrambled:
                        data.signal_cache_scramble.emplace_back(text);
                    break;
                    case Dataset::accumSummary_type::accumSummary_dest::enumerated:
                        data.signal_cache_enumerate.emplace_back(text);
                    break;
                    case Dataset::accumSummary_type::accumSummary_dest::none:
                        cerr << "none shouldn't happen!!" << endl;
                        exit(1);
                    break;
                    default:
                        cerr << "default shouldn't happen!!" << endl;
                        exit(1);
                    break;
                }

                // sqlite3_free((char*)text);
                text = NULL;
                
            }
            else{
                // not found (?)
                message = sqlite3_bind_text(seen_query, 1, kmer.c_str(),
                          -1, NULL);
                problemEncountered(message, "bind_text for seen_query");
                // int sqlite3_bind_int(sqlite3_stmt*, int, int);
                message = sqlite3_bind_int(seen_query, 2, data.settings.iteration);
                problemEncountered(message, "bind_int for seen_query");

                num_col = sqlite3_column_count(seen_query);
                if(num_col < 1){
                    cerr << "num_col for seen_query is less than 1!!\n\tEXITING" << '\n';
                    exit(1);
                }
                message = sqlite3_step(seen_query);
                isRowReady(message);

                // message is holding the int value now, 
                // count[0], not a code
                message = sqlite3_column_int(seen_query, 0);
                // cout << message << endl;
                if(message > 0){
                    // don't print for processing
                    // cout << "no print" << endl;
                }
                else{
                    // cout << "print" << endl;
                    message = sqlite3_bind_text(staged_query, 1, kmer.c_str(),
                                static_cast<int>(kmer.size()), NULL);
                    problemEncountered(message, "bind_text for staged_query");
                    // int sqlite3_bind_int(sqlite3_stmt*, int, int);
                    message = sqlite3_bind_int(staged_query, 2, data.settings.iteration);
                    problemEncountered(message, "bind_int for staged_query");

                    // return by reference
                    out_cache.push_back(kmer);

                    message = sqlite3_step(staged_query);
                    checkDone(message, "staged query execution line 171");
                    sqlite3_reset(staged_query);
                    sqlite3_clear_bindings(staged_query);
                }

                sqlite3_reset(seen_query);
                sqlite3_clear_bindings(seen_query);
            }


            sqlite3_reset(seen_query);
            sqlite3_clear_bindings(seen_query);
            // sqlite3_finalize(seen_query);
            // sqlite3_finalize(data_query);
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
        problemEncountered(message, "create unique index on kmer_cache");
        sqlite3_free(z_err_msg);
        
        msg = "CREATE TABLE seen_cache (kmer TEXT PRIMARY KEY NOT NULL, iter INT NOT NULL)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create table seen_cache");
        sqlite3_free(z_err_msg);
        
        msg = "CREATE UNIQUE INDEX seenIDX ON seen_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
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

            checkDone(message, "build_statement");

            sqlite3_reset(staged_query);
            sqlite3_clear_bindings(staged_query);
        }
	//cout << "Finalizing staged_query" << endl;
        message = sqlite3_finalize(staged_query);
        
        problemEncountered(message, "finalize staged_query");
    }

    message = sqlite3_close_v2(cacheDB);
    problemEncountered(message, "closing the connection");
}

static void problemEncountered(const int message, const string &what){
    if(message != SQLITE_OK){
        cerr << "Problem encountered with " << what << " !\n\tEXITING\n";
        cerr << "\tError code: " << message << endl;
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

// static const char* convert_to_const_char(const unsigned char* store){
//     return reinterpret_cast<const char*>(store);
// }
