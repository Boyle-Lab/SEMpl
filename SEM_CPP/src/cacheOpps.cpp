// These functions will use cached computations to save us some compute time
// Refactor of checkCache and writeCache into a more cohesive set of functions

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


// This builds a new cache if it does not exist and connects to the cache
void connectCache(Dataset &data, const string &cachefile, sqlite3 &cacheDB) {

    bool newcache = fileExists(cachefile);
    int message = 0;
    message = sqlite3_open_v2(cachefile.c_str(), &cacheDB,
                              SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
                              NULL);
    problemEncountered(message, "open");

    // Keep journal in memory (may want to load entire DB into memory)
    sqlite3_exec(cacheDB, "PRAGMA journal_mode = MEMORY", NULL, NULL, NULL);

    if(newcache) {
        if(data.settings.verbose){
            cout << "Existing cache used.\n" << flush;
        }
    } else {
        if(data.settings.verbose){
            cout << "No existing cache, creating new cache for use.\n" << flush;
        }

        char *z_err_msg = NULL;
        //cout << "Creating Table kmer_cache" << flush;
        msg = "CREATE TABLE kmer_cache (kmer TEXT NOT NULL, alignment BLOB)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        if(message != SQLITE_OK){
            cerr << "\tproblem with exec on create TABLE kmer_cache\n";
            exit(1);
        }
        sqlite3_free(z_err_msg);

        //cout << "Creating Unique Index kmerIDX" << flush;
        msg = "CREATE INDEX kmerIDX ON kmer_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create index on kmer_cache");
        sqlite3_free(z_err_msg);

        msg = "CREATE TABLE seen_cache (kmer TEXT PRIMARY KEY NOT NULL, iter INT NOT NULL)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create table seen_cache");
        sqlite3_free(z_err_msg);

        msg = "CREATE UNIQUE INDEX seenIDX ON seen_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);
    }

}

// Closes the connection to the cacheDB
void closeCache(Dataset &data, const string &cachefile, sqlite3 &cacheDB) {
    message = sqlite3_close_v2(cacheDB);
    problemEncountered(message, "closing the connection");
}

// MODIFIES: adds appropriate kmers to the specific output_cache
//           see line 121, switch statement
// EFFECTS: checks cache located at cachefile
//       to_align is the argument given to -out_cache in the original algorithm
//       cachefile is the argument given to -cache in the original algorithm
//       -out_file is built into the function within the switch statements
// IMPORTANT: to_align is the kmers that need to be aligned to genome
//            signal_cache_whatever are the alignments!!!
void checkCache(Dataset &data, vector<string> &in_file, vector<string> &to_align,
                sqlite3 &cacheDB, Dataset::accumSummary_type::accumSummary_dest dest,
                int position, char bp){

    vector<string> signal_cache_data;
    to_align.clear();

    if(data.settings.verbose){
        cout << "Querying cache for processed kmers..." << flush;
    }

    string msg = "";

    // Prepare SQL queries
//   msg = "SELECT count(*) FROM seen_cache WHERE kmer=? AND iter!=?";
    msg = "SELECT count(*) FROM seen_cache WHERE kmer=?";
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

  	//utilize sqlite transactions to speed this all up
	  sqlite3_exec(cacheDB, "BEGIN TRANSACTION", NULL, NULL, NULL);

    for(const string &kmer : in_file){
        message = sqlite3_bind_text(cache_signal_data_query, 1, kmer.c_str(),
                  -1, SQLITE_TRANSIENT);
        problemEncountered(message, "bind_text for cache_signal_data_query");
        message = sqlite3_step(cache_signal_data_query);

//      grabs the alignment of current kmers
        char* text = (char*)sqlite3_column_text(cache_signal_data_query, 1);

        if(text){  //The kmer was found so keep all alignments for this kmer
            do {
                signal_cache_data.emplace_back(text);
                text = NULL;

                message = sqlite3_step(cache_signal_data_query);
                text = (char*)sqlite3_column_text(cache_signal_data_query, 1);
            } while (message == SQLITE_ROW);
        }
        else{
            // not found

            //NOTE: I think that we don't need to do this query since we are INSERT OR IGNORE for the insert
            message = sqlite3_bind_text(amount_seen_query, 1, kmer.c_str(),
                      -1, SQLITE_TRANSIENT);
            problemEncountered(message, "bind_text for amount_seen_query");
            message = sqlite3_step(amount_seen_query);
            isRowReady(message);

            // message is holding the int value now,
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

        sqlite3_reset(cache_signal_data_query);
        sqlite3_clear_bindings(cache_signal_data_query);
    }

    sqlite3_exec(cacheDB, "COMMIT TRANSACTION", NULL, NULL, NULL);

	  //probably should finalize all of these - this is a memory leak otherwise
    message = sqlite3_finalize(insert_into_seen_cache_query);
    problemEncountered(message, "finalize insert_into_seen_cache_query");
    message = sqlite3_finalize(amount_seen_query);
    problemEncountered(message, "finalize amount_seen_query");
    message = sqlite3_finalize(cache_signal_data_query);
    problemEncountered(message, "finalize cache_signal_data_query");

    // store computations found in cache
    switch (dest) {
        case Dataset::accumSummary_type::accumSummary_dest::alignment:
            swap(data.signal_cache[ {position, bp} ] , signal_cache_data);
        break;
        case Dataset::accumSummary_type::accumSummary_dest::scrambled:
            swap(data.signal_cache_scramble, signal_cache_data);
        break;
        case Dataset::accumSummary_type::accumSummary_dest::enumerated:
            swap(data.signal_cache_enumerate, signal_cache_data);
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

    if(data.settings.verbose){
        cout << "FINISH" << endl;
    }
}

// REQUIRES: accumSummary_scale is filled with the correct data
// EFFECTS: writes output of accumSummary_scale to cache, based upon dest
void writeCache(Dataset &data, sqlite3 &cacheDB,
                Dataset::accumSummary_type::accumSummary_dest dest){


    // points to an output vector from running accumSummary_scale(args)
    const vector<string> *ptr = nullptr;
    // const vector<double> *max_ptr = nullptr;
    switch (dest) {
        case Dataset::accumSummary_type::accumSummary_dest::none:
            cerr << "dest shouldn't be none!!!!" << endl;
            exit(1);
        break;
        case Dataset::accumSummary_type::accumSummary_dest::enumerated:
            ptr = &data.accumSummary_data.enum_accum_lines;
            // max_ptr = &data.accumSummary_data.enum_accum_max;
        break;
        case Dataset::accumSummary_type::accumSummary_dest::scrambled:
            ptr = &data.accumSummary_data.scramble_accum_lines;
            // max_ptr = &data.accumSummary_data.scramble_accum_max;
        break;
        case Dataset::accumSummary_type::accumSummary_dest::alignment:
            ptr = &data.accumSummary_data.align_accum_lines;
            // max_ptr = &data.accumSummary_data.align_accum_max;
        break;
        default:
            cerr << "there is no default for dest's switch statement!!!" << endl;
            exit(1);
        break;
    }


    //utilize sqlite transactions to speed this all up
    sqlite3_exec(cacheDB, "BEGIN TRANSACTION", NULL, NULL, NULL);

    // build queries
    sqlite3_stmt *staged_query = nullptr;
    msg = "INSERT INTO kmer_cache VALUES(?,?)";
    sqlite3_prepare_v2(cacheDB, msg.c_str(), -1,
                       &staged_query, NULL);
    problemEncountered(message, msg);

    // Remove later per note below
    msg = "SELECT count(*) FROM kmer_cache WHERE kmer=? AND alignment=?";
    sqlite3_stmt *amount_seen_query = NULL;
    message = sqlite3_prepare_v2(cacheDB, msg.c_str(),
                                 static_cast<int>(msg.size()),
                                 &amount_seen_query, NULL);
    problemEncountered(message, msg);

    string temp = "", val = "";
    #ifdef DEBUG
    ofstream debug("written.txt");
    #endif

    for(size_t idx = 0; idx < ptr->size(); ++idx){

        val = ptr->at(idx);

        grab_string_4_index(val, temp);

        #ifdef DEBUG
        debug << "\tkmer: " << "first: #" << temp << "#\tsecond:#" << val << '#' << endl;
        #endif

        //currently we have kmers that have been seen in the input to this function
        // so these need to be filtered out here.
        // This should be edited later to not ever re-query these
        message = sqlite3_bind_text(amount_seen_query, 1, temp.c_str(),
                  -1, SQLITE_TRANSIENT);
        problemEncountered(message, "bind_text for amout_seen_query");

        message = sqlite3_bind_text(amount_seen_query, 2, val.c_str(),
                  -1, SQLITE_TRANSIENT);
        problemEncountered(message, "bind_text for amout_seen_query");
        message = sqlite3_step(amount_seen_query);

        #ifdef DEBUG
        if(sqlite3_column_type(amount_seen_query, 0) != SQLITE_INTEGER){
            cerr << "incorrect column_type on amount_seen_query\n\tEXITING"
                 << endl;
            exit(1);
         }
         #endif

        message = sqlite3_column_int(amount_seen_query, 0);

        if(message > 0) {
            //something found
        #ifdef DEBUG
            cerr << "DUP: " << message << endl;
        #endif
        } else {
            message = sqlite3_bind_text(staged_query, 1, temp.c_str(), -1, SQLITE_TRANSIENT);
            problemEncountered(message, "bind text 1 for staged_query, writeCache");
            message = sqlite3_bind_text(staged_query, 2, val.c_str(), -1, SQLITE_TRANSIENT);
            problemEncountered(message, "bind text 2 for staged_query, writeCache");
            message = sqlite3_step(staged_query);
            if(message != SQLITE_DONE){
                cerr << "Build statement is not done!\n\tEXITING" << endl;
                exit(1);
            }
            sqlite3_reset(staged_query);
            sqlite3_clear_bindings(staged_query);
        }

        sqlite3_reset(amount_seen_query);
        sqlite3_clear_bindings(amount_seen_query);
    }

    sqlite3_exec(cacheDB, "COMMIT TRANSACTION", NULL, NULL, NULL);

    message = sqlite3_finalize(staged_query);
    problemEncountered(message, "finalize staged_query");

    message = sqlite3_finalize(amount_seen_query);
    problemEncountered(message, "finalize amount_seen_query");

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
