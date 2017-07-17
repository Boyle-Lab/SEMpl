
// generateSNPEffectMatrix.pl --PWM <PWM_file>
//    --merge_file <merge file> --big_wig <bigwig file>
//    --TF_name <TF_name> [options] ...
// Required Options:
//  --PWM PWM input file
//  --merge_file DNAse-seq file
//  --big_wig ChIP-seq signal file
//  --TF_name TF name
/*

my $savedCmd = join(" ", @ARGV);

#other command line options
my $PWM;
my $DNase;
my $Signal; #in big wig format
my $TFname;
my $Threshold;
my $OutputFolder; #defaults to results/TFname/
my $Cache;

# Default Options
my $delSNPList = 1;
my $delAlignmentBed = 1;
my $delFilteredBed = 1;
my $delSignalFile = 0;
my $writecache = 0;
my $fastrun = 0;
my $iteration = -1;
*/
// ./generateSNPEffectMatrix.pl -PWM $PWM -merge_file $dnase -big_wig $chip
//     -TF_name $tf -output $output -threshold $threshold -iteration $iterID
//     -writecache -readcache $CACHE -verbose";

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
  //
  // f(data.output_dir.empty()){
  // data.output_dir = "results/" + data.TF_name + "/";
  //
  //
  // ifstream test_file(data.output_dir);
  // if(test_file){
  //   // directory exists
  // }
  // else{
  //   string s = "mkdir -p ";
  //   s += data.output_dir;
  //   system(s.c_str());
  // }

	if(data.cachefile.empty()){
		data.cachefile = data.output_dir + "/CACHE.db";
	}

// main script content
//---------------------------------------------------------

    if(data.settings.verbose){
		cout << "\nGenerating SEM for " << data.TF_name << endl;
	}

	//Step 1: Generate Enumerated k-mers
	cout << "Creating enumerated kmers from PWM file" << endl;
    cout << "\tstep one" << endl;
    int length = generate_kmers(data);
    // generate_kmers(data);
    // data.kmerHash is now filled in!!!!

	//Step 2: Change one base at each location in k-mers and align to genome
    // ALSO: print output to file


    // temporarily commented for testing
    // temporarily commented for testing
    // temporarily commented for testing
    cout << "\tstep two" << endl;
    align_to_genome(data);

    //Step 3: Filter using DNase data and finding the signal at each location
    // ALSO: read in output of filterDNaseWrapper back to memory
    // writes the *_filtered files
    if(data.settings.verbose){
        cout << "Filtering using DNase data and finding the signal" << endl;
    }
    // temporarily commented for testing
    // temporarily commented for testing
    // temporarily commented for testing
    cout << "\tstep three" << endl;
    filterDNaseWrapper(data);

    //Step 4: Find the signal using chIP-seq data
    cout << "\tstep four" << endl;
    find_signal(data, length);

    //Step 5: Generate baselines
    cout << "\tstep five" << endl;
    create_baselines(data, length);

    //Step 6: Create R plot(s) and a SEM output
    cout << "\tstep six" << endl;
    generate_output(data);

    if(data.settings.verbose){
        cout << "The SNP Effect Matrix has been completed for " << data.TF_name << endl;
    }

}

int generate_kmers(Dataset &data){
    if(data.settings.verbose){
        cout << "Creating enumerated kmers" << endl;
    }

    Enumerate_kmer(data);
  // data.kmerHash is now filled in!!!

    return Dataset::PWM::NUM_ROWS;

  // convert_PWM_format.pl is effectively performed within Enumerate_kmer(args)

  // threshold is stored within data.settings.threshold

  // length = number of lines in example transcription factor(?) file - 2;
}

void align_to_genome(Dataset &data){
    if(data.settings.verbose){
        cout << "Aligning SNPs in kmers to the genome\n";
    }
        // align all to genome
    alignToGenomeWrapper(data, data.settings.iteration, "./data/hg19");
}

// assumes filterDNaseWrapper_output is filled from previous function
// and that the output is sorted, and contains only unique string values

// NOTE: iteratively constructs the output from 
//       findMaximumAverageSignalWrapper(args),
//       as opposed to aggregately, as done in the original algorithm
void find_signal(Dataset &data, int length){
    if(data.settings.verbose){
        cout << "Finding the average signal" << endl;
    }
    // this function, in the original algorithm, iterates through files, where those
    // files correspond to the length of something, then each nucleotide letter
    // A, C, T, G

    // the original implementation seems to gather input from filterDNaseWrapper,
    // and print whatever output there is to bedfile. The input from filterDNaseWrapper
    // is sorted, and unique

    // signal is big_wig
    data.accumSummary_data.align_accum_lines.clear();
    data.accumSummary_data.align_accum_max.clear();

    vector<string> files;
    string cachefile = "";

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
    int value = 0;

    for(const string &file : files){
        if(file.find("filtered") == string::npos){
            continue;
        }
        if(data.settings.verbose){
            cout << "\tfile: " << file << endl;
        }

        size_t val = file.find("_pos");
        if(val == string::npos){
            cerr << "\"_pos\" not found within file\n\tEXITING";
            exit(1);
        }

        try{
            // grabbing bp and position from file name
            bp = file[val - 1];
            pos = file.c_str() + val - 1 + 5;
            end = pos;
            while(*end != '.'){
                ++end;
            }

            arr = new char[end - pos + 1];
            arr[end - pos] = '\0';
            strncpy(arr, file.c_str() + val - 1 + 5, end - pos);
            value = atoi(arr);
        }
        catch(...){
            cerr << "pointer errors likely\n\tEXITING" << endl;
            exit(1);
        }

        // cout << "\t" << file << "\n";

        // cin.get();

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
        try{
            if(data.settings.verbose){
                cout << "\twriteCache(args) is running..." << flush;
            }
            writeCache(data, data.cachefile,
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
            sort(data.signal_cache.begin(), data.signal_cache.end());
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }
            data.signal_output.resize(data.signal_cache.size()
                                    + data.accumSummary_data.align_accum_lines.size());
            // returns iterator to one past the location of the last copy
            if(data.settings.verbose){
                cout << "\tcopying..." << flush; 
            }
            auto iter = copy(data.accumSummary_data.align_accum_lines.begin(),
                             data.accumSummary_data.align_accum_lines.end(),
                             data.signal_output.begin());
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }

            //  FILLS data.signal_enumerate_output !!!!!!!!!!!!
            if(data.settings.verbose){
                cout << "\tunique copying..." << flush;
            }
            unique_copy(data.signal_cache_enumerate.begin(),
                        data.signal_cache_enumerate.end(),
                        iter);
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }
        }
        catch(...){
            cerr << "problem with algorithm usage" << endl;
            exit(1);
        }

        try{
            if(data.settings.verbose){
                cout << "\tfinding findMaximumAverageSignal..." << flush;
            }
            findMaximumAverageSignalWrapper(data,
                                            Dataset::accumSummary_type::accumSummary_dest::alignment);
            if(data.settings.verbose){
                cout << "FINISH" << endl;
            }

            // iteratively fills sig_deets


            // filled Signal_data data for appropriate type,
            // (alignment, scrambled, or enumerated)

            // first pair type

#ifdef DEBUG
            if(data.settings.verbose){
                cout << "\tpos: " << value << "  bp: " << bp << endl;
            }
#endif

            auto iter = data.sig_deets_maximum.insert( { {value, bp}, 
                                                data.Signal_data.alignment_maximum} );
            
                // cout << '\t' << data.Signal_data.alignment_maximum << "  max inserted" << endl;
            
            if(!iter.second){
                cerr << "duplicate key inserted into sig_deets_maximum" << endl;
                exit(1);
            }
            // second pair type
            auto iter1 = data.sig_deets_counter.insert( { {value, bp}, 
                                                data.Signal_data.alignment_counter} );
                // cout << '\t' << data.Signal_data.alignment_counter << "  counter inserted" << endl;
            if(!iter1.second){
                cerr << "duplicate key inserted into sig_deets_counter" << endl;
                exit(1);
            }
            iter = data.sig_deets_stdev.insert( { {value, bp}, 
                                                data.Signal_data.alignment_stdev} );
                // cout << '\t' << data.Signal_data.alignment_stdev << "  stdev inserted" << endl;
            if(!iter.second){
                cerr << "duplicate key inserted into sig_deets_stdev" << endl;
                exit(1);
            }
            iter = data.sig_deets_sterr.insert( { {value, bp}, 
                                                data.Signal_data.alignment_sterr} );
            if(!iter.second){
                cerr << "duplicate key inserted into sig_deets_sterr" << endl;
                exit(1);
            }

        }
        catch(...){
            cerr << "problem with findMaximumAverageSignalWrapper(args)" << endl;
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
        scramble_kmers.push_back(pair.first);
        // cout << "first: " << pair.first << endl << "second: " << pair.second << endl;
    }
    enumerate_kmers = scramble_kmers;
    // scramble_kmers and enumerate_kmers have same elements

    if(!data.settings.fastrun){
        for(auto it = scramble_kmers.begin(); 
            it != scramble_kmers.end(); 
            ++it){
            random_shuffle(it->begin(), it->end());
        }
        // scramble_kmers IS NOW SCRAMBLED!!!!

        checkCache(data, scramble_kmers, scramble_cache_to_align, data.cachefile,
                    Dataset::accumSummary_type::accumSummary_dest::scrambled);
        seq_col_to_fa(scramble_cache_to_align,
                      data.output_dir + "/BASELINE/Scrambled_kmer.fa");
        bowtie_genome_map(length, "../data/hg19",
                          data.output_dir + "/BASELINE/Scrambled_kmer.fa",
                          data.output_dir + "/BASELINE/Scrambled_kmer.bed",
                          data.settings.verbose);

        // NEED TO CHECK THAT THIS IS THE RIGHT RELATIVE DIRECTORY
        cmd = "./bin/bedtools intersect -a " + data.output_dir
                + "/BASELINE/Scrambled_kmer.bed -b "
                + data.DNase_file + " -wa -u > " + data.output_dir
                + "/BASELINE/Scrambled_kmer_filtered.bed";
        system(cmd.c_str());

        if(data.settings.verbose){
            cout << "\tRunning accumSummary_scale(args)..." << flush;
        }
        try{
            accumSummary_scale(data, data.bigwig_file,
                               data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed",
                               length,
                               Dataset::accumSummary_type::accumSummary_dest::scrambled);
        }
        catch(...){
            cerr << "problem with scrambled accumSummary_scale(args)\n\tEXITING" << endl;
            exit(1);
        }
        if(data.settings.verbose){
            cout << "FINISHed running accumSummary_scale(args) on scrambled" << endl;
        }

        if(data.settings.writecache){
            if(data.settings.verbose){
                cout << "\twriting to cache..." << flush;
            }
            try{
                writeCache(data, data.cachefile,
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
        auto iter = copy(data.accumSummary_data.scramble_accum_lines.begin(),
                         data.accumSummary_data.scramble_accum_lines.end(),
                         data.signal_scramble_output.begin());
        //  FILLS data.signal_scramble_output !!!!!!!!!!!!
        // iter is the next position to have a value inserted at of
        // data.signal_scramble_output
        unique_copy(data.signal_cache_scramble.begin(),
                    data.signal_cache_scramble.end(),
                    iter);

    } // !data.settings.fastrun


    // checkCache(Dataset &data, const std::vector<std::string> &in_file,
    //            std::vector<std::string> &out_cache, const std::string &cachefile,
    //            Dataset::accumSummary_type::accumSummary_dest dest)

    checkCache(data, enumerate_kmers, enumerate_cache_to_align, data.cachefile,
               Dataset::accumSummary_type::accumSummary_dest::enumerated);

    if(!enumerate_cache_to_align.empty()){
        seq_col_to_fa(enumerate_cache_to_align,
                      data.output_dir + "/BASELINE/Enumerated_kmer.fa");
        bowtie_genome_map(length, "../data/hg19", 
                          data.output_dir + "/BASELINE/Enumerated_kmer.fa",
                          data.output_dir + "/BASELINE/Enumerated_kmer_filtered.bed",
                          data.settings.verbose);
        accumSummary_scale(data, data.bigwig_file,
                           data.output_dir + "/BASELINE/Enumerated_kmer_filtered.bed",
                           length, Dataset::accumSummary_type::accumSummary_dest::enumerated);
        if(data.settings.writecache){
            writeCache(data, data.cachefile,
                       Dataset::accumSummary_type::accumSummary_dest::enumerated);
        }
    }

    // SHOULD THERE BE AN ERROR CHECK IF signal_cache_enumerate IS EMPTY????
    if(data.signal_cache_enumerate.empty()){
        cerr << "\tdata.signal_cache_enumerate is empty!!!! IDK" << endl;
    }
    sort(data.signal_cache_enumerate.begin(),
         data.signal_cache_enumerate.end());
    data.signal_enumerate_output.resize(data.signal_cache_enumerate.size()
                                      + data.accumSummary_data.enum_accum_lines.size());
    // returns iterator to one past the location of the last copy
    auto iter = copy(data.accumSummary_data.enum_accum_lines.begin(),
                     data.accumSummary_data.enum_accum_lines.end(),
                     data.signal_enumerate_output.begin());

    //  FILLS data.signal_enumerate_output !!!!!!!!!!!!
    unique_copy(data.signal_cache_enumerate.begin(),
                data.signal_cache_enumerate.end(),
                iter);

    findMaximumAverageSignalWrapper(data,
                    Dataset::accumSummary_type::accumSummary_dest::enumerated);
    findMaximumAverageSignalWrapper(data,
                    Dataset::accumSummary_type::accumSummary_dest::scrambled);


}

void generate_output(Dataset &data){
    if (data.settings.verbose){
        cout << "Generating Output" << endl;
    }

    // generateSEM(data);

    if(!data.settings.fastrun){
        generateRplot(data);
        quality_control(data);
    }
}
