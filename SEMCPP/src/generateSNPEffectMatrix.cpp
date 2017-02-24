
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

#include "iterativeSEM.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
#include "common.h"
using namespace std;

void find_signal(Dataset &data);
void create_baselines(Dataset &data, int length);
void align_to_genome(Dataset &data);
void generate_output(Dataset &data);
int generate_kmers(Dataset &data);
void Enumerate_kmer(Dataset &data);


void generateSNPEffectMatrix(Dataset &data){
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
	cout << "Creating enumerated kmers from PWM file" << '\n';

    int length = generate_kmers(data);
    // data.kmerHash is now filled in!!!!

	//Step 2: Change one base at each location in k-mers and align to genome
    if(data.settings.verbose){
        cout << "Aligning SNPs in kmers to the genome" << '\n';
    }
    align_to_genome(data);

    //Step 3: Filter using DNase data and finding the signal at each location
    if(data.settings.verbose){
        cout << "Filtering using DNase data and finding the signal" << '\n';
    }
    filterDNaseWrapper(data);

    //Step 4: Find the signal using chIP-seq data
    find_signal(data);

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
    cout << "Creating enumerated kmers" << '\n';
  }

  Enumerate_kmer(data);
  // data.kmerHash is now filled in!!!

  return Dataset::PWM::NUM_ROWS;

  // convert_PWM_format.pl is effectively performed within Enumerate_kmer(args)

  // threshold is stored within data.settings.threshold

  // length = number of lines in example transcription factor(?) file - 2;
}

void align_to_genome(Dataset &data){
        if(data.settings.verbose) cout << "Aligning SNPs in kmers to the genome\n";
        // align all to genome
        alignToGenomeWrapper(data, data.settings.iteration);
}

// assumes filterDNaseWrapper_output is filled from previous function
// assumes filterDNaseWrapper_output is filled from previous function
// assumes filterDNaseWrapper_output is filled from previous function
// and that the output is sorted, and contains only unique string values
// and that the output is sorted, and contains only unique string values
// and that the output is sorted, and contains only unique string values
void find_signal(Dataset &data){ //This function call needs to be checked or altered
    if(data.settings.verbose){
        cout << "Finding the average signal" << '\n';
    }
    // this function, in the original algorithm, iterates through files, where those
    // files correspond to the length of something, then each nucleotide letter
    // A, C, T, G

    // the original implementation seems to gather input from filterDNaseWrapper,
    // and print whatever output there is to bedfile. The input from filterDNaseWrapper
    // is sorted, and unique

    // I will need to understand the correct output of filterDNaseWrapper, then
    // I can work on this. The exact output of filterDNaseWrapper is important

    // http://bedtools.readthedocs.io/en/latest/content/example-usage.html

    // perl passes variables to subroutines by alias
    // chip is the file corresponding to the -big_wig argument
    // chip is a bigwig file
    // Signal is chip
    // first argument to accumSummary_scale is signal, then filteredBedfile, then length
    //                                  bigwig file,    output of filterDNaseWrapper,  length
    //                                      hfile,          cfile,              scale

    for(auto val : data.filterDNaseWrapper_output){
        //accumSummary_scale(data, const std::string &hfile, const std::string &cfile, int scale);
        // will have to see how Colten implements bedtools, then make a decision about this
        // and accumSummary_scale
    }



    //Processes bigwig
    // string targetDir = data.output_dir + "/ALIGNMENT/";
    // string filteredBedfile = "";
    // // I don't think this can be done, I think you need to list the files, then
    // // open each individually
    // ifstream directory(targetDir);
    // string file = "";
    // while (directory >> file ){
    //     if(file == "/pos/"){
    //         targetDir += file + "/";
    //         filteredBedfile = targetDir + file + "_filtered.bed";
    //         //accumSummary_scale(data);
    //         if(data.settings.writecache){
    //             writeCache(data);
    //         }
    //         // Not sure on the implementation of the following line currently
    //     }
    //     string cachefile = targetDir + file + "signal.cache";
    //
    // }
    // directory.close();
    //
    // if (data.settings.verbose){
    //     cout << "Creating directory SIGNAL " << '\n';
    // }
    // string cmd = "rm -rf " + data.output_dir + "/SIGNAL";
    //
    //
    // system(cmd.c_str());
    // cmd = "mkdir " + data.output_dir + "/SIGNAL";
    //
    // system(cmd.c_str());
    // cmd = "cp " + data.output_dir
    //      + "/ALIGNMENT/*/*/signal " + data.output_dir + "/SIGNAL/";
    //
    // system(cmd.c_str());
    //
    //
    // //build signal summary
    // findMaximumAverageSignalWrapper(data);
    //
    // // cmd = "rm " + data.output_dir + "/SIGNAL/*.signal";
    // //
    // // system(cmd.c_str());
}

void create_baselines(Dataset &data, int length){
    if (data.settings.verbose){
        cout << "Creating directory BASELINE" << '\n';
    }
    string cmd = "rm -rf " + data.output_dir + "/BASELINE";

    system(cmd.c_str());
    cmd = "mkdir " + data.output_dir + "/BASELINE";

    system(cmd.c_str());

    //Create baseline from scrambled k-mers
    cmd = "cut -f2 " + data.output_dir + "/Enumerated_kmer.txt > "
            + data.output_dir + "/BASELINE/Enumerated_kmer.txt";

    system(cmd.c_str());

    if(!data.settings.fastrun){
        //scramble_kmer(data);
        //checkCache(data);
        //seq_col_to_fa(data, 0);
        //bowtie_genome_map(data);
        cmd = "./bin/bedtools intersect -a " + data.output_dir + "/BASELINE/Scrambled_kmer.bed -b "
                + data.DNase_file + " -wa -u > "+ data.output_dir
                + "/BASELINE/Scrambled_kmer_filtered.bed";
        system(cmd.c_str());
        //accumSummary_scale(data, data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed", )
        if (data.settings.writecache){
            writeCache(data);
        }
    }

    cmd = "cat " + data.output_dir + "/BASELINE/Enumerated_kmer.cache | sort | uniq >> "
            + data.output_dir + "/BASELINE/Enumerated_kmer_filtered.signal";

    system(cmd.c_str());

    findMaximumAverageSignalWrapper(data);

    if(data.settings.delAlignmentBed){
        cmd = "rm -f " + data.output_dir + "/BASELINE/Scrambled_kmer.bed";

        system(cmd.c_str());
        cmd = "rm -f " + data.output_dir + "/BASELINE/Enumerated_kmer.bed";

        system(cmd.c_str());
    }

    if(data.settings.delFilteredBed){
        cmd = "rm -f " + data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed";

        system(cmd.c_str());
        cmd = "rm -f " + data.output_dir + "/BASELINE/Enumerated_kmer_filtered.bed";

        system(cmd.c_str());
    }

    // if(data.settings.delSNPList){
    //     cmd = "rm -f " + data.output_dir + "/BASELINE/*.scrambled";
    //
    //     system(cmd.c_str());
    //     cmd = "rm -f " + data.output_dir + "/BASELINE/*.fa";
    //
    //     system(cmd.c_str());
    //     cmd = "rm -f " + data.output_dir + "/BASELINE/*.sm.txt";
    //
    //     system(cmd.c_str());
    //     cmd = "rm -f " + data.output_dir + "/BASELINE/*.cache";
    //
    //     system(cmd.c_str());
    //     cmd = "rm -f " + data.output_dir + "/BASELINE/*.Enumerated_kmer.txt";
    //
    //     system(cmd.c_str());
    // }
}

void generate_output(Dataset &data){
    if (data.settings.verbose){
        cout << "Generating Output" << '\n';
    }

    generateSEM(data);

    if(!data.settings.fastrun){
        generateRplot(data);
        quality_control(data);
    }
}
