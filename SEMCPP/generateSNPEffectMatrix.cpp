
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
using namespace std;

void generateSNPEffectMatrix(Dataset &data){
	// default options are built into settings within data
	
	if(data.output_dir.empty()){
		data.output_dir = "results/" + data.TF_name + "/";
	}
	/*
	DIR* dir = opendir(data.output_dir);
	if(dir){
		// directory exists
	}
	else{
		system("mkdir -p " + data.output_dir.c_str());
	}
	closedir(dir);
	
	if(data.cache_file.empty()){
		data.cache_file = data.output_dir + "/CACHE.db";
	}
	*/

// main script content
//---------------------------------------------------------

	if(data.settings.verbose){
		cout << "\nGenerating SEM for " + data.TF_name + "\n"
			<< "Command: \n" << data.command << "\n";
	}

	// STEPS UNFINISHED

}

