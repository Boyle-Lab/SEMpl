// Cody Morterud and Colten Williams


// Required Options:
//  --PWM PWM input file
//  --merge_file DNAse-seq file
//  --big_wig ChIP-seq signal file
//  --TF_name TF name

#include "./src/iterativeSEM.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include "./src/common.hpp"
using namespace std;

void find_signal(Dataset &data, int length);
void create_baselines(Dataset &data, int length);
void align_to_genome(Dataset &data);
void generate_output(Dataset &data);
int generate_kmers(Dataset &data);
void Enumerate_kmer(Dataset &data);

void generateSNPEffectMatrix(Dataset &data) {
    // default options are built into settings within data
    // in iterativeSEM.hpp

    // Set a default cache for scores if none has been provided
    if(data.cachefile.empty()){
        data.cachefile = data.base_dir + "/CACHE.db";
    }

    if(data.settings.verbose){
        cout << "\nGenerating SEM for " << data.TF_name << endl;
    }

    if(data.settings.verbose){
        cout << "\tBuilding directories" << endl;
    }
    string cmd = "mkdir -p " + data.output_dir;
    system(cmd.c_str());


    // clears the result caches, as in the caches that get filled
    // with pre-computed kmers
    data.signal_cache.clear();
    data.signal_cache_scramble.clear();
    data.signal_cache_enumerate.clear();

    connectCache(data, data.cachefile, data.cacheDB);

    //-------------------------------------
    // Primary compute steps:
    //-------------------------------------

    //Step 1: Generate Enumerated k-mers
    cout << "Creating enumerated kmers from PWM file" << endl;
    cout << "\tstep one" << endl;
    int length = generate_kmers(data); // data.kmerHash is now filled in!!!!

    //Step 2: Change one base at each location in k-mers and align to genome
    // ALSO: print output to file
//NOTE: We may want to not write this to a file. Also it might be a speedup if we go straight to signal instead of doing all of
//  the kmers here first then doing signal because of duplicate lookups.
//  Check how many times this would be duplicated in our test.

    cout << "\tstep two" << endl;
    align_to_genome(data);

    //Step 3: Find the signal using ChIP-seq data
    cout << "\tstep three" << endl << flush;
    find_signal(data, length);

    //Step 4: Generate baselines
    cout << "\tstep four" << endl;
    create_baselines(data, length);

    //Step 5: Create R plot(s) and a SEM output
    cout << "\tstep five" << endl;
    generate_output(data);

    closeCache(data, data.cacheDB);

    if(data.settings.verbose){
        cout << "The SNP Effect Matrix has been completed for " << data.TF_name << endl;
    }

}


// Generates an initial list of kmers based on the provided PWM and threshold
int generate_kmers(Dataset &data){
    if(data.settings.verbose){
        cout << "Creating enumerated kmers" << endl;
    }

    Enumerate_kmer(data); // data.kmerHash is now filled in!!!

    // length of kmers is constant for a single run of the program
    return data.PWM_data.matrix_arr[0].size();

   // threshold is stored within data.settings.threshold
}

// Calls wrapper program to align kmers and stores this
void align_to_genome(Dataset &data){
    if(data.settings.verbose){
        cout << "Aligning SNPs in kmers to the genome\n";

    }

    // align all to genome
    alignToGenomeWrapper(data, data.settings.iteration, "./data/hg19");
}


// assumes alignment is filled from previous function
// and that the output is sorted, and contains only unique string values

// NOTE: iteratively constructs the output from
//       findMaximumAverageSignalWrapper(args),
//       as opposed to aggregately, as done in the original algorithm
void find_signal(Dataset &data, int length){
    if(data.settings.verbose){
        cout << "Finding the average signal" << endl;
    }

    // signal is big_wig
    data.accumSummary_data.align_accum_lines.clear();

    vector<string> files;
    string cachefile = "";

//NOTE: We can keep all of this in memory instead of writing the alignments to file
    GetFilesInDirectory(files, data.output_dir + "/ALIGNMENT/");

#ifdef DEBUG
    sort(files.begin(), files.end());
#endif

    char bp = '\0';
    const char *pos = nullptr;
    const char *end = nullptr;

    // clear data.sig_deets*
    data.sig_deets_maximum.clear();
    data.sig_deets_counter.clear();
    data.sig_deets_stdev.clear();
    data.sig_deets_sterr.clear();

    char* arr = nullptr;
    int position = 0;

    for(const string &file : files){
        if(file.find("bed") == string::npos){
            continue;
        }
        if(data.settings.verbose){
            cout << "\tfile: " << file << endl;
        }


        // grabbing bp and position from file name
        size_t val = file.find("_pos");
        if(val == string::npos){
            cerr << "\"_pos\" not found within file\n\tEXITING";
            exit(1);
        }
        try{
            bp = file[val - 1];
            pos = file.c_str() + val - 1 + 5;
            end = pos;
            while(*end != '.'){
                ++end;
            }

            arr = new char[end - pos + 1];
            arr[end - pos] = '\0';
            strncpy(arr, file.c_str() + val - 1 + 5, end - pos);
            position = atoi(arr);
        }
        catch(...){
            cerr << "pointer errors likely\n\tEXITING" << endl;
            exit(1);
        }

        // Calculate summary score for all alignments
        try{
            if(data.settings.verbose){
                cout << "\taccumSummary_scale(args) is running..." << flush;
            }
            accumSummary_scale(data, data.bigwig_file, file, length,
                               Dataset::accumSummary_type::accumSummary_dest::alignment);
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }
        }
        catch(...){
            cerr << "Problem with accumSummary_scale\n\tEXITING" << endl;
            exit(1);
        }
        // do not use N/A files

        // Write NEW computed scores to our cache so that we don't need to do previous steps again
        try{
            if(data.settings.verbose){
                cout << "\twriteCache(args) is running..." << flush;
            }
            writeCache(data, data.cacheDB,
                       Dataset::accumSummary_type::accumSummary_dest::alignment);
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }
        }
        catch(...){
            cerr << "problem with writeCache\n\tEXITING" << endl;
            exit(1);
        }

        // SHOULD THERE BE AN ERROR CHECK IF signal_cache_enumerate IS EMPTY????
        try{
            if(data.settings.verbose){
                cout << "\tsorting..." << flush;
            }
            sort(data.signal_cache[ {position, bp} ].begin(),
                 data.signal_cache[ {position, bp} ].end() );
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }
            data.signal_output.resize(data.signal_cache[ {position, bp} ].size()
                                    + data.accumSummary_data.align_accum_lines.size());
            // returns iterator to one past the location of the last copy
            if(data.settings.verbose){
                cout << "\tcopying..." << flush;
            }
            auto iter = data.signal_output.begin();
            if(data.accumSummary_data.align_accum_lines.size() != 0 ){
                iter = copy(data.accumSummary_data.align_accum_lines.begin(),
                                 data.accumSummary_data.align_accum_lines.end(),
                                 data.signal_output.begin() );
            }
            if(data.settings.verbose){
                cout << "FINISH (calculated)" << data.accumSummary_data.align_accum_lines.size() << endl;
            }

            //  FILLS data.signal_enumerate_output
            if(data.settings.verbose){
                cout << "\tunique copying..." << flush;
            }
            auto end_iter3 = unique_copy(data.signal_cache[ {position, bp} ].begin(),
                        data.signal_cache[ {position, bp} ].end(),
                        iter);
            if(data.settings.verbose){
                cout << "FINISH (cache)" << data.signal_cache[ {position, bp} ].size() << endl;
            }
            data.signal_output.resize(end_iter3 - data.signal_output.begin() );
        }
        catch(...){
            cerr << "problem with algorithm usage" << endl;
            exit(1);
        }

        #ifdef DEBUG
            //cerr << "\tsignal_output" << endl;
            // for(auto val : data.signal_output){
            //     cerr << "\tval: #" << val << '#' << endl;
            // }
            //cerr << "\tone" << endl;
            //ofstream one("one.txt");
            //for(auto val : data.signal_cache[ {position, bp} ] ){
            //    one << val << endl;
            //}
            //ofstream two("two.txt");
            //cerr << "\ttwo" << endl;
            //for(auto val : data.accumSummary_data.align_accum_lines){
            //    two << val << endl;
            //}

        #endif

        // Now summarize the scores from all alignments (including those from cache)
        // into a single maximum, counter, stdev, and sterr
        try{
            if(data.settings.verbose){
                cout << "\tfinding findMaximumAverageSignal..." << flush;
            }
            findMaximumAverageSignalWrapper(data.signal_output,
                                            data.Signal_data.alignment_maximum,
                                            data.Signal_data.alignment_counter,
                                            data.Signal_data.alignment_stdev,
                                            data.Signal_data.alignment_sterr);
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }
        }
        catch(...){
            cerr << "problem with findMaximumAverageSignalWrapper(args)" << endl;
            exit(1);
        }
            // iteratively fills sig_deets


            // filled Signal_data data for appropriate type,
            // (alignment, scrambled, or enumerated)

            // first pair type

#ifdef DEBUG
            if(data.settings.verbose){
                cout << "\tpos: " << position << "  bp: " << bp << endl;
            }
#endif
            try{
                auto iter = data.sig_deets_maximum.insert( { {position, bp},
                                                    data.Signal_data.alignment_maximum} );

                    // cout << '\t' << data.Signal_data.alignment_maximum << "  max inserted" << endl;

                if(!iter.second){
                    cerr << "duplicate key inserted into sig_deets_maximum" << endl;
                    exit(1);
                }
                // second pair type
                auto iter1 = data.sig_deets_counter.insert( { {position, bp},
                                                    data.Signal_data.alignment_counter} );
                    // cout << '\t' << data.Signal_data.alignment_counter << "  counter inserted" << endl;
                if(!iter1.second){
                    cerr << "duplicate key inserted into sig_deets_counter" << endl;
                    exit(1);
                }
                iter = data.sig_deets_stdev.insert( { {position, bp},
                                                    data.Signal_data.alignment_stdev} );
                    // cout << '\t' << data.Signal_data.alignment_stdev << "  stdev inserted" << endl;
                if(!iter.second){
                    cerr << "duplicate key inserted into sig_deets_stdev" << endl;
                    exit(1);
                }
                iter = data.sig_deets_sterr.insert( { {position, bp},
                                                    data.Signal_data.alignment_sterr} );
                if(!iter.second){
                    cerr << "duplicate key inserted into sig_deets_sterr" << endl;
                    exit(1);
                }
            }
            catch(...){
                cerr << "exception thrown when inserting alignment data" << endl;
                exit(1);
            }



        free(arr);
    }

}

void create_baselines(Dataset &data, int length){
    if (data.settings.verbose){
        cout << "Creating directory BASELINE" << endl;
    }
    string cmd = "rm -rf " + data.output_dir + "/BASELINE";

    system(cmd.c_str());
    cmd = "mkdir " + data.output_dir + "/BASELINE";

    system(cmd.c_str());

    //Create baseline from scrambled k-mers

    // grab keys of kmerHash



    // scrambled kmers
    vector<string> scramble_kmers;
    scramble_kmers.reserve(data.kmerHash.size());

    // enumerated kmers
    vector<string> enumerate_kmers;

    // final output from processing scrambled kmer data
    vector<string> scramble_cache_to_align;
    // final output from processing enumerated kmer data
    vector<string> enumerate_cache_to_align;


    for(const auto &pair : data.kmerHash){
        enumerate_kmers.push_back(pair.first);
        // cout << "first: " << pair.first << endl << "second: " << pair.second << endl;
    }

    if(!data.settings.fastrun){
        scramble_kmers = enumerate_kmers;
        for(auto it = scramble_kmers.begin();
            it != scramble_kmers.end();
            ++it){
            random_shuffle(it->begin(), it->end());
        }
        // scramble_kmers IS NOW SCRAMBLED!!!!

        checkCache(data, scramble_kmers, scramble_cache_to_align, data.cacheDB,
                    Dataset::accumSummary_type::accumSummary_dest::scrambled);
        seq_col_to_fa(scramble_cache_to_align,
                      data.output_dir + "/BASELINE/Scrambled_kmer.fa");
        bowtie_genome_map(length, "../data/hg19",
                          data.output_dir + "/BASELINE/Scrambled_kmer.fa",
                          data.output_dir + "/BASELINE/Scrambled_kmer.bed",
                          data.DNase_file, data.settings.verbose);

        // NEED TO CHECK THAT THIS IS THE RIGHT RELATIVE DIRECTORY
        //cmd = "./bin/bedtools intersect -a " + data.output_dir
        //        + "BASELINE/Scrambled_kmer.bed -b "
        //        + data.DNase_file + " -wa -u > " + data.output_dir
        //        + "BASELINE/Scrambled_kmer_filtered.bed";
        //system(cmd.c_str());

        if(data.settings.verbose){
            cout << "\tRunning accumSummary_scale(args)..." << flush;
        }
        try{
            accumSummary_scale(data, data.bigwig_file,
                               data.output_dir + "BASELINE/Scrambled_kmer.bed",
                               length,
                               Dataset::accumSummary_type::accumSummary_dest::scrambled);
        }
        catch(...){
            cerr << "problem with scrambled accumSummary_scale(args)\n\tEXITING" << endl;
            exit(1);
        }
        if(data.settings.verbose){
            cout << "FINISH" << endl;
        }

        if(data.settings.writecache){
            if(data.settings.verbose){
                cout << "\twriting to cache..." << flush;
            }
            try{
                writeCache(data, data.cacheDB,
                           Dataset::accumSummary_type::accumSummary_dest::scrambled);
            }
            catch(...){
                cerr << "problem with writeCache(args)\n\tEXITING" << endl;
                exit(1);
            }
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }
        }
        // SHOULD THERE BE AN ERROR CHECK IF signal_cache_enumerate IS EMPTY????
        // SIGNAL CACHE SCRAMBLE DOES NOT EXIST OMG

        data.signal_scramble_output.clear();

        sort(data.signal_cache_scramble.begin(),
             data.signal_cache_scramble.end());
        data.signal_scramble_output.resize(data.signal_cache_scramble.size()
                                          + data.accumSummary_data.scramble_accum_lines.size());
        // returns iterator to one past the location of the last copy
        auto iter = data.signal_scramble_output.begin();
        if(data.accumSummary_data.scramble_accum_lines.size() != 0){
            iter = copy(data.accumSummary_data.scramble_accum_lines.begin(),
                             data.accumSummary_data.scramble_accum_lines.end(),
                             data.signal_scramble_output.begin());
        }
        //  FILLS data.signal_scramble_output !!!!!!!!!!!!
        // iter is the next position to have a value inserted at of
        // data.signal_scramble_output
        auto end_iter2 = unique_copy(data.signal_cache_scramble.begin(),
                    data.signal_cache_scramble.end(),
                    iter);
        data.signal_scramble_output.resize(end_iter2 - data.signal_scramble_output.begin());
        #ifdef DEBUG
            ofstream debug1(data.output_dir + "/BASELINE/Scrambled_kmer.signal");
            for(auto val : data.signal_scramble_output){
                debug1 << val << endl;
            }
            debug1.close();
        #endif

    } // !data.settings.fastrun


    // checkCache(Dataset &data, const std::vector<std::string> &in_file,
    //            std::vector<std::string> &out_cache, const std::string &cachefile,
    //            Dataset::accumSummary_type::accumSummary_dest dest)

    // PWM kmer hits from enumerate kmers
    checkCache(data, enumerate_kmers, enumerate_cache_to_align, data.cacheDB,
               Dataset::accumSummary_type::accumSummary_dest::enumerated);

//    if(data.settings.verbose){
//        cout << "Aligned (cache): " << data.signal_cache[ {position, bp} ].size() << endl;
//    }

    if(!enumerate_cache_to_align.empty()){
        seq_col_to_fa(enumerate_cache_to_align,
                      data.output_dir + "/BASELINE/Enumerated_kmer.fa");
        bowtie_genome_map(length, "../data/hg19",
                          data.output_dir + "/BASELINE/Enumerated_kmer.fa",
                          data.output_dir + "/BASELINE/Enumerated_kmer.bed",
                          data.DNase_file, data.settings.verbose);
        try{

            accumSummary_scale(data, data.bigwig_file,
                               data.output_dir + "/BASELINE/Enumerated_kmer.bed",
                              length, Dataset::accumSummary_type::accumSummary_dest::enumerated);
        }
        catch(...){
            cerr << "problem with accumSummary_scale on enumerated!!\n\tEXITING" << endl;
            exit(1);
        }
        if(data.settings.writecache){
            writeCache(data, data.cacheDB,
                       Dataset::accumSummary_type::accumSummary_dest::enumerated);
        }
    }

    sort(data.signal_cache_enumerate.begin(),
         data.signal_cache_enumerate.end());
    data.signal_enumerate_output.resize(data.signal_cache_enumerate.size()
                                      + data.accumSummary_data.enum_accum_lines.size());
    // returns iterator to one past the location of the last copy
    auto iter = data.signal_enumerate_output.begin();
    if( data.accumSummary_data.enum_accum_lines.size() != 0 ){
        iter = copy(data.accumSummary_data.enum_accum_lines.begin(),
                         data.accumSummary_data.enum_accum_lines.end(),
                         data.signal_enumerate_output.begin());
    }

    //  FILLS data.signal_enumerate_output !!!!!!!!!!!!
    auto end_iter1 = unique_copy(data.signal_cache_enumerate.begin(),
                data.signal_cache_enumerate.end(),
                iter);
    data.signal_enumerate_output.resize(end_iter1 - data.signal_enumerate_output.begin());

    #ifdef DEBUG
        ofstream debug(data.output_dir + "/BASELINE/Enumerated_kmer.signal");
        for(auto val : data.signal_enumerate_output){
            debug << val << endl;
        }
        debug.close();
    #endif



    if(!data.settings.fastrun){
        if(data.settings.verbose) cout << "scram" << flush;
        findMaximumAverageSignalWrapper(data.signal_scramble_output,
                                        data.Signal_data.scramble_maximum,
                                        data.Signal_data.scramble_counter,
                                        data.Signal_data.scramble_stdev,
                                        data.Signal_data.scramble_sterr);
    }
    #ifdef DEBUG
    ofstream enum_out("Enumerated_kmer.signal");
    for(auto val : data.signal_enumerate_output){
        enum_out << val << endl;
    }
    enum_out.close();
    #endif


    if(data.settings.verbose) cout << "enum" << flush;
    findMaximumAverageSignalWrapper(data.signal_enumerate_output,
                                    data.Signal_data.enumerate_maximum,
                                    data.Signal_data.enumerate_counter,
                                    data.Signal_data.enumerate_stdev,
                                    data.Signal_data.enumerate_sterr);

}

void generate_output(Dataset &data){
    if (data.settings.verbose){
        cout << "Generating Output" << endl;
    }

    generateSEM(data);

    if(!data.settings.fastrun){
        generateRplot(data);
        quality_control(data);
    }
}
