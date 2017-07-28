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

  // takes in_file, out_file, to_align, cache(location of cache), current iteration,
  // all of them are file names, and are specific to nucleotide and position

  // signal_cache in Dataset is the output cache
  // cache is the input cache

// EFFECTS: returns true if file already exists in current directory,
//          false otherwise
static void problemEncountered(const int message, const string &what);
static void isRowReady(const int message);
static void checkDone(const int message, const string &s);


// MODIFIES: adds appropriate kmers to the specific output_cache
//           see line 121, switch statement
// EFFECTS: checks cache located at cachefile
//       to_align is the argument given to -out_cache in the original algorithm
//       cachefile is the argument given to -cache in the original algorithm
//       -out_file is built into the function within the switch statements
// IMPORTANT: to_align is the kmers that need to be aligned to genome
//            signal_cache_whatever are the alignments!!!
void checkCache(Dataset &data, const vector<string> &in_file, vector<string> &to_align,
                const string &cachefile, Dataset::accumSummary_type::accumSummary_dest dest,
                int position, char bp){

    vector<string> signal_cache_data;


    bool newcache = fileExists(cachefile);

    to_align.clear();

    if(data.settings.verbose){
        cout << "Querying cache for processed kmers..." << flush;
    }
    // should this be here?
    // ANS: I don't think so, this is already processed data
    switch (dest) {
        case Dataset::accumSummary_type::accumSummary_dest::alignment:
        break;
        case Dataset::accumSummary_type::accumSummary_dest::scrambled:
            data.signal_cache_scramble.clear();
        break;
        case Dataset::accumSummary_type::accumSummary_dest::enumerated:
            data.signal_cache_enumerate.clear();
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

    sqlite3 *cacheDB = NULL;
    int message = 0;
    message = sqlite3_open_v2(cachefile.c_str(), &cacheDB, 
                              SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
                              NULL);
    problemEncountered(message, "open");

    string msg = "";

    if(newcache){
        if(data.settings.verbose){
            // cout << "Cache does exist\n";
        }

        msg = "SELECT count(*) FROM seen_cache WHERE kmer=? AND iter!=?";
        sqlite3_stmt* amount_seen_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), 
                                     static_cast<int>(msg.size()), 
                                     &amount_seen_query, NULL);
        problemEncountered(message, msg);

        msg = "SELECT kmer, alignment FROM kmer_cache WHERE kmer=?";
        sqlite3_stmt* cache_signal_data_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), 
                                     static_cast<int>(msg.size()), 
                                     &cache_signal_data_query, NULL);
        problemEncountered(message, msg);

        msg = "INSERT OR IGNORE INTO seen_cache VALUES(?, ?)";
        sqlite3_stmt* insert_into_seen_cache_query = NULL;
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), 
                                     static_cast<int>(msg.size()), 
                                     &insert_into_seen_cache_query, NULL);
        problemEncountered(message, msg);

        // int num_col = 0;

        for(const string &kmer : in_file){
            // cout << kmer.c_str() << endl;
            message = sqlite3_bind_text(cache_signal_data_query, 1, kmer.c_str(),
                      -1, SQLITE_TRANSIENT);
            problemEncountered(message, "bind_text for cache_signal_data_query");

            message = sqlite3_step(cache_signal_data_query);

            // num_col = sqlite3_column_count(cache_signal_data_query);
#ifdef DEBUG
            // cout << "There are " << num_col << " columns in cache_signal_data_query\n";
            // if(num_col < 0) {
            //     cerr << "Number of columns from cache_signal_data_query is less than 0!!\n\tEXITING" << endl;
            //     exit(1);
            // }
            // if(num_col != 2) {
            //     cerr << "Number of columns from cache_signal_data_query is not 2!!\n\tEXITING" << endl;
            //     exit(1);
            // }
#endif
            // CHECKING THE NUMBER OF COLS WORKS ONLY FOR THE PERL API
            // I BELIEVE THIS API RETURNS A NULL IF SOMETHING DOESN'T EXIST
            // UPDATE: GOOGLE THIS

//          grabs the alignment of current kmers
            const char* text = (char*)sqlite3_column_text(cache_signal_data_query, 1);
            
            if(text){
                #ifdef DEBUG
                // cerr << "\tfound: #" << kmer << '#' << endl
                     // << "\tcorresponding align: #" << text << '#' << endl;
                #endif
                signal_cache_data.emplace_back(text);
                
                text = NULL;
                
            }
            else{
                // not found
                #ifdef DEBUG
                // cerr << "not found: " << kmer << endl;
                #endif
                message = sqlite3_bind_text(amount_seen_query, 1, kmer.c_str(),
                          -1, SQLITE_TRANSIENT);
                problemEncountered(message, "bind_text for amount_seen_query");
                // int sqlite3_bind_int(sqlite3_stmt*, int, int);
                message = sqlite3_bind_int(amount_seen_query, 2, data.settings.iteration);
                problemEncountered(message, "bind_int for amount_seen_query");

                // num_col = sqlite3_column_count(amount_seen_query);
                // if(num_col < 1){
                //     cerr << "num_col for amount_seen_query is less than 1!!\n\tEXITING" << '\n';
                //     exit(1);
                // }
                message = sqlite3_step(amount_seen_query);
                isRowReady(message);

                // message is holding the int value now, 
                // count[0], not a code
                #ifdef DEBUG
                    if(sqlite3_column_type(amount_seen_query, 0) != SQLITE_INTEGER){
                        cerr << "incorrect column_type on amount_seen_query\n\tEXITING"
                             << endl;
                        exit(1);
                    }
                #endif
                message = sqlite3_column_int(amount_seen_query, 0);
                // cout << message << endl;
                if(message > 0){
                    // don't print for processing
                    #ifdef DEBUG
                    // cerr << "no print" << endl;
                    #endif
                }
                else{
                    #ifdef DEBUG
                    // cerr << "print" << endl;
                    #endif
                    message = sqlite3_bind_text(insert_into_seen_cache_query, 1, kmer.c_str(),
                                static_cast<int>(kmer.size()), SQLITE_TRANSIENT);
                    problemEncountered(message, "bind_text for insert_into_seen_cache_query");
                    // int sqlite3_bind_int(sqlite3_stmt*, int, int);
                    message = sqlite3_bind_int(insert_into_seen_cache_query, 2, data.settings.iteration);
                    problemEncountered(message, "bind_int for insert_into_seen_cache_query");

                    // return by reference
                    to_align.push_back(kmer);

                    message = sqlite3_step(insert_into_seen_cache_query);
                    checkDone(message, "staged query execution line 197");
                    sqlite3_reset(insert_into_seen_cache_query);
                    sqlite3_clear_bindings(insert_into_seen_cache_query);
                }

                sqlite3_reset(amount_seen_query);
                sqlite3_clear_bindings(amount_seen_query);
            }

            sqlite3_reset(amount_seen_query);
            sqlite3_clear_bindings(amount_seen_query);
            sqlite3_reset(cache_signal_data_query);
            sqlite3_clear_bindings(cache_signal_data_query);
        }

        switch (dest) {
            case Dataset::accumSummary_type::accumSummary_dest::alignment:
                data.signal_cache[ {position, bp} ] = signal_cache_data;
            break;
            case Dataset::accumSummary_type::accumSummary_dest::scrambled:
                data.signal_cache_scramble = signal_cache_data;
            break;
            case Dataset::accumSummary_type::accumSummary_dest::enumerated:
                data.signal_cache_enumerate = signal_cache_data;
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

        sqlite3_stmt *insert_into_seen_cache_query = NULL;
        msg = "INSERT OR IGNORE INTO seen_cache VALUES(?,?)";
        message = sqlite3_prepare_v2(cacheDB, msg.c_str(), 
                                     static_cast<int>(msg.size()), 
                                     &insert_into_seen_cache_query, NULL);
        problemEncountered(message, msg);

        if(insert_into_seen_cache_query == NULL){
            cerr << "insert_into_seen_cache_query shouldn't be NULL\n";
            exit(1);
        }
	//cout << "Binding to cache" << endl;

        to_align = in_file;
        for(const auto &kmer : in_file){
            message = sqlite3_bind_text(insert_into_seen_cache_query, 1, 
                      kmer.c_str(), -1, SQLITE_TRANSIENT);
            problemEncountered(message, "bind text for inserting into seen_cache");
            message = sqlite3_bind_int(insert_into_seen_cache_query, 2, data.settings.iteration);
            problemEncountered(message, "bind int for inserting into seen_cache");

            message = sqlite3_step(insert_into_seen_cache_query);

            checkDone(message, "insert_into_seen_cache_query");

            sqlite3_reset(insert_into_seen_cache_query);
            sqlite3_clear_bindings(insert_into_seen_cache_query);
        }
	//cout << "Finalizing insert_into_seen_cache_query" << endl;
        message = sqlite3_finalize(insert_into_seen_cache_query);
        
        problemEncountered(message, "finalize insert_into_seen_cache_query");
    }

    message = sqlite3_close_v2(cacheDB);
    problemEncountered(message, "closing the connection");

    if(data.settings.verbose){
        cout << "FINISH" << endl;
    }
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
