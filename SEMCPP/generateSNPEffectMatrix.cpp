
// generateSNPEffectMatrix.pl --PWM <PWM_file> --merge_file <merge file> --big_wig <bigwig file> --TF_name <TF_name> [options] ...
// Required Options:
//  --PWM PWM input file
//  --merge_file DNAse-seq file
//  --big_wig ChIP-seq signal file
//  --TF_name TF name
/*

=head1 DESCRIPTION

This is the current version of the technique to generate a SNP Effect Matrix.

=head1 OPTIONS

=over 8

=item B<-h, --help>

Print this brief help message from the command line.

=item B<-d, --debug>

Print debug output.

=item B<-v, --verbose>

Verbose output.

=item B<--output <target>>

Define the output location [Default: results/TFname/]

=item B<--Threshold <value>>

Set PWM cutoff threshold (if not a pre-computed PWM)

=item B<--delSNP>

delete SNP file [default]

=item B<--delAlign>

delete Alignment Bed file

=item B<--delFiltered>

delete DNase filtered Bed file [default]

=item B<--delSignal>

delete Signal files

=back

=head1 OUTPUT

 (1) <TF_name>.sem -- has the numerical values of the SNP Effect Matrix

 (2) <TF_name>_semplot.pdf -- graphical representation of the SNP Effect Matrix

 Bed and signal files from intermediate steps can also be kept

/////////////////////////////////////////

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
// ./generateSNPEffectMatrix.pl -PWM $PWM -merge_file $dnase -big_wig $chip -TF_name $tf -output $output -threshold $threshold -iteration $iterID -writecache -readcache $CACHE -verbose";

#include "iterativeSEM.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>
using namespace std;

void find_signal(Dataset &data);
void create_baselines(Dataset &data);
void generate_output(Dataset &data);
int generate_kmers(Dataset &data);
void Enumerate_kmer(Dataset &data);

string revCompDNA(string);


void generateSNPEffectMatrix(Dataset &data){
	// default options are built into settings within data

	if(data.output_dir.empty()){
		data.output_dir = "results/" + data.TF_name + "/";
	}

  ifstream test_file(data.output_dir);
  if(test_file){
    // directory exists
  }
  else{
    string s = "mkdir -p ";
    s += data.output_dir.c_str();
    system(s.c_str());
  }

	if(data.cache_file.empty()){
		data.cache_file = data.output_dir + "/CACHE.db";
	}

// main script content
//---------------------------------------------------------

	if(data.settings.verbose){
		cout << "\nGenerating SEM for " + data.TF_name + "\n"
			<< "Command: \n" << data.command << "\n";
	}

	//Step 1: Generate Enumerated k-mers
	cout << "Creating enumerated kmers from PWM file" << endl;

    Enumerate_kmer(data);
    // data.kmerHash is now filled in!!!!

	//Step 2: Change one base at each location in k-mers and align to genome
    if(data.settings.verbose){
        cout << "Aligning SNPs in kmers to the genome" << endl;
    }
    alignToGenomeWrapper(data);

    //Step 3: Filter using DNase data and finding the signal at each location
    if(data.settings.verbose){
        cout << "Filtering using DNase data and finding the signal" << endl;
    }
    filterDNaseWrapper(data);

    //Step 4: Find the signal using chIP-seq data
    find_signal(data);

    //Step 5: Generate baselines
    create_baselines(data);

    //Step 6: Create R plot(s) and a SEM output
    generate_output(data);

    if(data.settings.verbose){
        cout << "The SNP Effect Matrix has been completed  for " << data.TF_name << endl;

    }

}

int generate_kmers(Dataset &data){
  if(data.settings.verbose) cout << "Creating enumerated kmers" << '\n';

  Enumerate_kmer(data);
  // data.kmerHash is now filled in!!!

  return Dataset::PWM::NUM_ROWS;

  // convert_PWM_format.pl is effectively performed within Enumerate_kmer(args)

  // threshold is stored within data.settings.threshold

  // length = numer of lines in example transcription factor(?) file - 2;
}

void align_to_genome(Dataset &data, int iteration){
  if(data.settings.verbose) cout << "Aligning SNPs in kmers to the genome\n";
  // align all to genome
  alignToGenomeWrapper(data);
}

void find_signal(Dataset &data){ //This function call needs to be checked or altered
    if(data.settings.verbose){
        cout << "Finding the average signal" << endl;
    }

    //Processes bigwig
    string targetDir = data.output_dir + "/ALIGNMENT/";
    string filteredBedfile = "";
    ifstream directory(targetDir);
    string file;
    while (directory >> file ){
        if(file == "/pos/"){
            targetDir += file + "/";
            filteredBedfile = targetDir + file + "_filtered.bed";
            accumSummary_scale(data);
            if(data.settings.writecache){
                writeCache(data);
            }
            // Not sure on the implementation of the following line currently
        }
        string cachefile = targetDir + file + "signal.cache";

    }
    directory.close();

    if (data.settings.verbose){
        cout << "Creating directory SIGNAL " << endl;
    }
    string cmd = "rm -rf " + data.output_dir +"/SIGNAL";
    char *convert = new char[cmd.length() +1];
    strcpy(convert, cmd.c_str());
    system(convert);
    cmd = "mkdir " + data.output_dir +"/SIGNAL";
    strcpy(convert, cmd.c_str());
    system(convert);
    cmd = "cp " + data.output_dir +"/ALIGNMENT/*/*/signal " + data.output_dir + "/SIGNAL/";
    strcpy(convert, cmd.c_str());
    system(convert);


    //build signal summary
    findMaximumAverageSignalWrapper(data);

    cmd = "rm " + data.output_dir + "/SIGNAL/*.signal";
    strcpy(convert, cmd.c_str());
    system(convert);
}

void create_baselines(Dataset &data){
    if (data.settings.verbose){
        cout << "Creating directory Baseline" << endl;
    }
    string cmd = "rm -rf " + data.output_dir +"/BASELINE";
    char *convert = new char[100];
    strcpy(convert, cmd.c_str());
    system(convert);
    cmd = "mkdir " + data.output_dir +"/BASELINE";
    strcpy(convert, cmd.c_str());
    system(convert);

    //Create baseline from scrambled k-mers
    cmd = "cut -f2 " + data.output_dir + "/Enumerated_kmer.txt > " + data.output_dir + "/BASELINE/Enumerated_kmer.txt";
    strcpy(convert, cmd.c_str());
    system(convert);

    if(!data.settings.fastrun){
        scramble_kmer(data);
        checkCache(data);
        seq_col_to_fa(data);
        bowtie_genome_map(data);
        cmd = "./bin/bedtools intersect -a " + data.output_dir + "/BASELINE/Scrambled_kmer.bed -b " + data.DNase_file + " -wa -u > "+ data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed";
        strcpy(convert, cmd.c_str());
        system(convert);
        //accumSummary_scale(data, data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed", )
        if (data.settings.writecache){
            writeCache(data);
        }
    }

    cmd = "cat " + data.output_dir + "/BASELINE/Enumerated_kmer.cache | sort | uniq >> " + data.output_dir + "/BASELINE/Enumerated_kmer_filtered.signal";
    strcpy(convert, cmd.c_str());
    system(convert);
    findMaximumAverageSignalWrapper(data);

    if(data.settings.delAlignmentBed){
        cmd = "rm -f " + data.output_dir + "/BASELINE/Scrambled_kmer.bed";
        strcpy(convert, cmd.c_str());
        system(convert);
        cmd = "rm -f " + data.output_dir + "/BASELINE/Enumerated_kmer.bed";
        strcpy(convert, cmd.c_str());
        system(convert);
    }

    if(data.settings.delFilteredBed){
        cmd = "rm -f " + data.output_dir + "/BASELINE/Scrambled_kmer_filtered.bed";
        strcpy(convert, cmd.c_str());
        system(convert);
        cmd = "rm -f " + data.output_dir + "/BASELINE/Enumerated_kmer_filtered.bed";
        strcpy(convert, cmd.c_str());
        system(convert);
    }

    if(data.settings.delSNPList){
        cmd = "rm -f " + data.output_dir + "/BASELINE/*.scrambled";
        strcpy(convert, cmd.c_str());
        system(convert);
        cmd = "rm -f " + data.output_dir + "/BASELINE/*.fa";
        strcpy(convert, cmd.c_str());
        system(convert);
        cmd = "rm -f " + data.output_dir + "/BASELINE/*.sm.txt";
        strcpy(convert, cmd.c_str());
        system(convert);
        cmd = "rm -f " + data.output_dir + "/BASELINE/*.cache";
        strcpy(convert, cmd.c_str());
        system(convert);
        cmd = "rm -f " + data.output_dir + "/BASELINE/*.Enumerated_kmer.txt";
        strcpy(convert, cmd.c_str());
        system(convert);
    }
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
