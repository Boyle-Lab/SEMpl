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

void create_baselines(Dataset &data, int length);
void align_to_genome(Dataset &data);
void generate_output(Dataset &data);
int generate_kmers(Dataset &data);
void Enumerate_kmer(Dataset &data);

void generateSNPEffectMatrix(Dataset &data) {
    // default options are built into settings within data
    // in iterativeSEM.hpp

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

    //Step 4: Generate baselines
    cout << "\tstep four" << endl;
    create_baselines(data, length);

    //Step 5: Create R plot(s) and a SEM output
    cout << "\tstep five" << endl;
    generate_output(data);

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
                          data.DNase_file, data.settings.threads, data.settings.verbose);

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
                          data.DNase_file, data.settings.threads, data.settings.verbose);
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
        //generateRplot(data);
        quality_control(data);
    }
}
