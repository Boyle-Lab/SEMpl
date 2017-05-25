#include "iterativeSEM.hpp"
#include "../lib/sqlite3/sqlite3.h"
#include "common.hpp"
//#include <sstream>
using namespace std;

bool fileExists(const string &filename);
static void problemEncountered(const int &message, const string &what);
//static void isRowReady(const int &message);
// static void prepareStmt(sqlite3 *db, string stmt, sqlite3_stmt *query);
// static void checkDone(const int &message, const string &s);

// REQUIRES: accumSummary_scale is filled with the correct data
// EFFECTS: writes output of accumSummary_scale to cache, based upon dest
void writeCache(Dataset &data, const string &cache,
                Dataset::accumSummary_type::accumSummary_dest dest){

    if(data.settings.verbose){
        cout << "Building cache for processed kmers.\n";
    }
    // wants the third space, indexed from 0
    map<string, string> kmers;

    const vector<string> *ptr = nullptr;
    switch (dest) {
        case Dataset::accumSummary_type::accumSummary_dest::none:
            cerr << "dest shouldn't be none!!!!\n";
            exit(1);
        break;
        case Dataset::accumSummary_type::accumSummary_dest::enumerated:
            ptr = &data.accumSummary_data.enum_accum_lines;
        break;
        case Dataset::accumSummary_type::accumSummary_dest::scrambled:
            ptr = &data.accumSummary_data.scramble_accum_lines;
        break;
        case Dataset::accumSummary_type::accumSummary_dest::alignment:
            ptr = &data.accumSummary_data.align_accum_lines;
        break;
        default:
            cerr << "there is no default for dest's switch statement!!!\n";
            exit(1);
        break;
    }

    for(string line : *ptr){
#ifdef DEBUG
        // cout << "string at index 3 of string: " << line
            // << "\nis: " << grab_string_at_index(line, 3, string("\t")) << endl;
#endif
        kmers[grab_string_at_index(line, 3, string("\t"))] = line;
    }

    bool newcache = fileExists(cache);

    sqlite3 *cacheDB;
    int message = 0;
    string msg;
    message = sqlite3_open(cache.c_str(), &cacheDB);
    problemEncountered(message, "open");

    if(newcache){

    }
    else{
        if(data.settings.verbose){
            cout << "No existing cache!\n";
        }
    // BLOB? why is it a BLOB? I used "text" not "blob", will need to ask about this

        // sqlite3_stmt *build_statement = nullptr;

        char *z_err_msg = NULL;

        msg = "CREATE TABLE kmer_cache (kmer TEXT PRIMARY KEY NOT NULL, alignment BLOB)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);
        // message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        // // prepareStmt(cacheDB, msg, build_statement);
        // // message = sqlite3_step(build_statement);
        // checkDone(message, "build statement create table kmer_cache");
        // sqlite3_finalize(build_statement);

        msg = "CREATE UNIQUE INDEX kmerIDX ON kmer_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);
        // message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        // // prepareStmt(cacheDB, msg, build_statement);
        // // message = sqlite3_step(build_statement);
        // checkDone(message, "build statement create unique index kmer_cache");
        // sqlite3_finalize(build_statement);

        msg = "CREATE TABLE seen_cache (kmer TEXT PRIMARY KEY NOT NULL, iter INT NOT NULL)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);
        // message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        // // prepareStmt(cacheDB, msg, build_statement);
        // // message = sqlite3_step(build_statement);
        // checkDone(message, "build statement create table seen_cache");
        // sqlite3_finalize(build_statement);

        msg = "CREATE UNIQUE INDEX seenIDX ON seen_cache(kmer)";
        message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        problemEncountered(message, "create unique index on seen_cache");
        sqlite3_free(z_err_msg);
        // message = sqlite3_exec(cacheDB, msg.c_str(), NULL, NULL, &z_err_msg);
        // // prepareStmt(cacheDB, msg, build_statement);
        // // message = sqlite3_step(build_statement);
        // checkDone(message, "build statement create unique index on seen_cache");
    }
    sqlite3_stmt *staged_query = nullptr;
    msg = "INSERT OR IGNORE INTO kmer_cache VALUES(?,?)";
    sqlite3_prepare_v2(cacheDB, msg.c_str(), static_cast<int>(msg.size()),
                       &staged_query, NULL);
    problemEncountered(message, msg);
    // prepareStmt(cacheDB, msg, staged_query);


    //stringstream ss;
    #ifdef DEBUG
    cout << "\tkmers" << endl;
    #endif

    for(auto val1 : kmers){
        #ifdef DEBUG
        cout << "\tkmer: " << val1.first << ' ' << val1.second << endl;
        #endif

        message = sqlite3_bind_text(staged_query, 1, val1.first.c_str(), -1, NULL);
        problemEncountered(message, "bind text 1 for staged_query, writeCache");
        message = sqlite3_bind_text(staged_query, 2, val1.second.c_str(), -1, NULL);
        problemEncountered(message, "bind text 2 for staged_query, writeCache");
        message = sqlite3_step(staged_query);
        if(message != SQLITE_DONE){
            cerr << "Build statement is not done!\n\tEXITING\n";
            exit(1);
        }
        sqlite3_reset(staged_query);
        sqlite3_clear_bindings(staged_query);
    }

    message = sqlite3_finalize(staged_query);
    problemEncountered(message, "finalize staged_query");

    message = sqlite3_close(cacheDB);
    problemEncountered(message, "closing the connection");
}



// static void prepareStmt(sqlite3 *db, string stmt, sqlite3_stmt *query){
//     //sqlite3_prepare(sqlite3 *db, const char *zSql, int nByte, sqlite3_stmt **ppStmt, const char **pzTail)
//     int message = sqlite3_prepare_v2(db, stmt.c_str(), static_cast<int>(stmt.size()), &query, nullptr);
//     problemEncountered(message, stmt);
// }

static void problemEncountered(const int &message, const string &what){
    if(message != SQLITE_OK){
        cerr << "Problem encountered with " << what << "!\n\tEXITING\n";
        cerr << "code: " << message << endl;
        exit(1);
    }
}

// static void checkDone(const int &message, const string &s){
//     if(message != SQLITE_DONE){
//         cerr << s << " is not done!\n\tEXITING" << endl;
//         exit(1);
//     }
// }

//  static void isRowReady(const int &message){
//      if(message != SQLITE_ROW){
//          cerr << "Row isn't ready!!\n\tEXITING\n";
//          exit(1);
//      }
//  }

