
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

#include "src/iterativeSEM.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include "src/common.hpp"
using namespace std;

void find_signal(Dataset &data, int length);
void create_baselines(Dataset &data, int length);
void align_to_genome(Dataset &data);
void generate_output(Dataset &data);
int generate_kmers(Dataset &data);
void Enumerate_kmer(Dataset &data);

void generateSNPEffectMatrix(Dataset &data) {
	// default options are built into settings within data
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
		cout << "\nGenerating SEM for " + data.TF_name + "\n"
			<< "Command: \n" << data.command << "\n";
	}

	//Step 1: Generate Enumerated k-mers
	cout << "Creating enumerated kmers from PWM file\n";

    int length = generate_kmers(data);
    // data.kmerHash is now filled in!!!!

	//Step 2: Change one base at each location in k-mers and align to genome
    // ALSO: print output to file
    if ( data.settings.verbose ) {
        cout << "Aligning SNPs in kmers to the genome\n";
    }
    align_to_genome(data);

    //Step 3: Filter using DNase data and finding the signal at each location
    // ALSO: read in output of filterDNaseWrapper back to memory
    if(data.settings.verbose){
        cout << "Filtering using DNase data and finding the signal\n";
    }
    filterDNaseWrapper(data);
    // data.filterDNaseWrapper_output is filled!!!!!!

    //Step 4: Find the signal using chIP-seq data
    find_signal(data, length);

    //Step 5: Generate baselines
    create_baselines(data, length);

    //Step 6: Create R plot(s) and a SEM output
    generate_output(data);

    if(data.settings.verbose){
        cout << "The SNP Effect Matrix has been completed  for " << data.TF_name << '\n';
    }

}

int generate_kmers(Dataset &data){
    if(data.settings.verbose){
        cout << "Creating enumerated kmers\n";
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
    alignToGenomeWrapper(data, data.settings.iteration, "../data/hg19");
}

// assumes filterDNaseWrapper_output is filled from previous function
// and that the output is sorted, and contains only unique string values
void find_signal(Dataset &data, int length){
    if(data.settings.verbose){
        cout << "Finding the average signal\n";
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

    for(const auto &file : files){
        accumSummary_scale(data, data.bigwig_file, file, length,
                           Dataset::accumSummaryData::accumSummary_dest::alignment);
#ifdef DEBUG
        cout << "\tDeleting " << file << '\n';
        int val = system(("rm -rf " + file).c_str());
        assert(val == 0);
#else
        system(("rm -rf " + file).c_str());
#endif
        // write to cache
        // -in_file and -cache are built into data
        writeCache(data, data.cachefile,
                   Dataset::accumSummaryData::accumSummary_dest::alignment);

        // SHOULD THERE BE AN ERROR CHECK IF signal_cache_enumerate IS EMPTY????
        sort(data.signal_cache.begin(), 
             data.signal_cache.end());
        data.signal_output.resize(data.signal_cache.size() 
                                + data.accumSummary_data.align_accum_lines.size());
        // returns iterator to one past the location of the last copy
        auto iter = copy(data.accumSummary_data.align_accum_lines.begin(),
                         data.accumSummary_data.align_accum_lines.end(),
                         data.signal_output.begin());

        // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
        // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
        // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
        // and corresponding other data

        //  FILLS data.signal_enumerate_output !!!!!!!!!!!!
        unique_copy(data.signal_cache_enumerate.begin(), 
                    data.signal_cache_enumerate.end(),
                    iter);


        findMaximumAverageSignalWrapper(data,
                                        Dataset::accumSummaryData::accumSummary_dest::alignment);

    }

}

void create_baselines(Dataset &data, int length){
    if (data.settings.verbose){
        cout << "Creating directory BASELINE\n";
    }
    string cmd = "rm -rf " + data.output_dir + "/BASELINE";

    system(cmd.c_str());
    cmd = "mkdir " + data.output_dir + "/BASELINE";

    system(cmd.c_str());

    //Create baseline from scrambled k-mers

    // grab keys of kmerHash


    std::vector<string> scramble_cache_output;

    for(const auto &pair : data.kmerHash){
        data.scramble_kmers.push_back(pair.first);
    }

    if(!data.settings.fastrun){
        scramble_kmer(data);
        // scramble_kmers IS NOW SCRAMBLED!!!!                          $Cache

        checkCache(data, data.scramble_kmers, scramble_cache_output, data.cachefile,
                    Dataset::accumSummaryData::accumSummary_dest::scrambled);
        seq_col_to_fa(data.signal_cache_scramble,
                      data.output_dir + "/BASELINE/Scrambled_kmer.fa");
        bowtie_genome_map(length, "../data/hg19",
                          data.output_dir + "/BASELINE/Scrambled_kmer.fa",
                          data.output_dir + "/BASELINE/Scrambled_kmer.bed");

        // NEED TO CHECK THAT THIS IS THE RIGHT RELATIVE DIRECTORY
        cmd = "./bin/bedtools intersect -a " + data.output_dir
                + "/BASELINE/Scrambled_kmer.bed -b "
                + data.DNase_file + " -wa -u > "+ data.output_dir
                + "/BASELINE/Scrambled_kmer_filtered.bed";
        system(cmd.c_str());

        accumSummary_scale(data, data.bigwig_file,
                           data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed",
                           length,
                           Dataset::accumSummaryData::accumSummary_dest::scrambled);

        if(data.settings.writecache){
            writeCache(data, data.cachefile, 
                       Dataset::accumSummaryData::accumSummary_dest::scrambled);
        }
        // SHOULD THERE BE AN ERROR CHECK IF signal_cache_enumerate IS EMPTY????
        sort(data.signal_cache_scramble.begin(), 
             data.signal_cache_scramble.end());
        data.signal_enumerate_output.resize(data.signal_cache_scramble.size() 
                                          + data.accumSummary_data.scramble_accum_lines.size());
        // returns iterator to one past the location of the last copy
        auto iter = copy(data.accumSummary_data.scramble_accum_lines.begin(),
                         data.accumSummary_data.scramble_accum_lines.end(),
                         data.signal_scramble_output.begin());

        // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
        // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
        // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
        // and corresponding other data

        //  FILLS data.signal_enumerate_output !!!!!!!!!!!!
        unique_copy(data.signal_cache_scramble.begin(), 
                    data.signal_cache_scramble.end(),
                    iter);

    } // !data.settings.fastrun


    // checkCache(Dataset &data, const std::vector<std::string> &in_file,
    //            std::vector<std::string> &out_cache, const std::string &cachefile,
    //            Dataset::accumSummaryData::accumSummary_dest dest)
    checkCache(data, data.scramble_kmers, data.signal_cache_enumerate, data.cachefile,
               Dataset::accumSummaryData::accumSummary_dest::enumerated);

    if(!data.signal_cache_enumerate.empty()){
        seq_col_to_fa(data.signal_cache_enumerate,
                      data.output_dir + "/BASELINE/Enumerated_kmer.fa");
        bowtie_genome_map(length,
                          "../data/hg19", data.output_dir + "/BASELINE/Enumerated_kmer.fa",
                          data.output_dir + "/BASELINE/Enumerated_kmer_filtered.bed");
        accumSummary_scale(data, data.bigwig_file,
                           data.output_dir + "/BASELINE/Enumerated_kmer_filtered.bed",
                           length, Dataset::accumSummaryData::accumSummary_dest::enumerated);
        if(data.settings.writecache){
            writeCache(data, data.cachefile,
                       Dataset::accumSummaryData::accumSummary_dest::enumerated);
        }
    }

    // SHOULD THERE BE AN ERROR CHECK IF signal_cache_enumerate IS EMPTY????
    sort(data.signal_cache_enumerate.begin(), 
         data.signal_cache_enumerate.end());
    data.signal_enumerate_output.resize(data.signal_cache_enumerate.size() 
                                      + data.accumSummary_data.enum_accum_lines.size());
    // returns iterator to one past the location of the last copy
    auto iter = copy(data.accumSummary_data.enum_accum_lines.begin(),
                     data.accumSummary_data.enum_accum_lines.end(),
                     data.signal_enumerate_output.begin());

    // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
    // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
    // CHANGE REQUIRED REGARDING WHAT IS IN DATA.SIGNAL_ENUMERATE..._OUTPUT
    // and corresponding other data

    //  FILLS data.signal_enumerate_output !!!!!!!!!!!!
    unique_copy(data.signal_cache_enumerate.begin(), 
                data.signal_cache_enumerate.end(),
                iter);

    findMaximumAverageSignalWrapper(data,
                                    Dataset::accumSummaryData::accumSummary_dest::enumerated);

    if(data.settings.delAlignmentBed){
        system(("rm -f " + data.output_dir + "/BASELINE/Scrambled_kmer.bed").c_str());
        system(("rm -f " + data.output_dir + "/BASELINE/Enumerated_kmer.bed").c_str());
    }
    if(data.settings.delFilteredBed){
        system(("rm -f " + data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed").c_str());
        system(("rm -f " + data.output_dir + "/BASELINE/Enumerated_kmer_filtered.bed").c_str());
    }
    if(data.settings.delSNPList){
        system(("rm -f  " + data.output_dir + "/BASELINE/*.scrambled").c_str());
		system(("rm -f  " + data.output_dir + "/BASELINE/*.fa").c_str());
		system(("rm -f  " + data.output_dir + "/BASELINE/*.sm.txt").c_str());
		system(("rm -f  " + data.output_dir + "/BASELINE/*.cache").c_str());
		system(("rm -f  " + data.output_dir + "/BASELINE/Enumerated_kmer.txt").c_str());
    }

}

void generate_output(Dataset &data){
    if (data.settings.verbose){
        cout << "Generating Output\n";
    }

    generateSEM(data);

    if(!data.settings.fastrun){
        generateRplot(data);
        quality_control(data);
    }
}
