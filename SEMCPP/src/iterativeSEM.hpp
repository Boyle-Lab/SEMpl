//
//  iterativeSEM.hpp
//  SEMCPP
//
//  Created by Cody Morterud and Colten Williams on 11/9/16.
//  Copyright © 2016 Boyle Lab. All rights reserved.
//

#ifndef iterativeSEM_hpp
#define iterativeSEM_hpp

#include <string>
#include <array>
#include <fstream>
#include <vector>
#include <cassert>
#include <map>

/*
 example execution from command line
 "./iterativeSEM.pl -PWM examples/MA0114.1.pwm
    -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak
    -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig
    -TF_name HNF4A -output examples/HNF4A/"
*/

// overview
// a struct to contain an instance of PWM data and DNase data
// data members made public for ease of access, otherwise wouldn't, given more time

struct Dataset {

    Dataset(const Dataset &other) = delete;
    Dataset& operator=(const Dataset &other) = delete;


    Dataset(){ }
	// overview
	// a struct to contain and manage the PWM data as given in the example file
	struct PWM{

	    static const int NUM_ROWS = 13, NUM_COLUMNS = 5;


	    // holds first three inputs as given in example
	    std::string first_input, second_input, third_input;
	    // holds the characters at the end of each matrix row,
        // holds the two characters at the end
	    std::string end_of_line_char, end_of_matrix_char;
	    // holds the integer values of the matrix
	    std::array<std::array<int, NUM_COLUMNS>, NUM_ROWS> matrix_arr;
	    // holds the modified SEM version of the PWM
	    std::array<std::array<int, NUM_COLUMNS>, NUM_ROWS> sem_arr;

	};
    /*
	struct DNase{

	    static const int LINES_IN_FILE = 116018;

	    // "chr" and the chromosome number
	    std::array<std::string, LINES_IN_FILE> chromosome;
	    // first two numbers given
	    std::array<int, LINES_IN_FILE> first_num, second_num;
	    // chromosome from above and the specific section
	    std::array<std::string, LINES_IN_FILE> chr_section;
	    // third number
	    std::array<int, LINES_IN_FILE> third_num;
	    // fourth and fifth numbers, double
	    std::array<double, LINES_IN_FILE> fourth_num, fifth_num;
	    // sixth and seventh number given
	    std::array<int, LINES_IN_FILE> sixth_num, seventh_num;


	};
    */
	struct TFMdata{
		// a c g t
		static const int LETTER_NUM = 4;
		// first letter is a, then c, then g, then t, at least in example
		std::array<std::vector<int>, LETTER_NUM> letter_array;
    // should be int or char? int for now
	};
	struct accumSummaryData{
	    //Alignment summary data
        enum class accumSummary_dest{alignment, scrambled, enumerated, none};
		// lines of output from accumSummary_scale.pl
		std::vector<std::string> align_accum_lines;
		// max of output from accumSummary_scale.pl
		std::vector<double> align_accum_max;

		//baseline scramble kmer summary data
		std::vector<std::string> scramble_accum_lines;
		std::vector<double> scramble_accum_max;

		//baseline enumerated kmer summary data
		std::vector<std::string> enum_accum_lines;
		std::vector<double> enum_accum_max;
	};
  // contains default settings
	struct SettingsForSNPEffectMatrix{
		bool delSNPList = true, delAlignmentBed = true, delFilteredBed = true;
		bool delSignalFile = false, writecache = false, fastrun = false,
             verbose = false;
        int iteration = -1;
		double threshold = 0.0;
	};
	// contains data from findMaximumAverageSignalWrapper
	struct MaximumAverageSignalData{
        //baseline signal data
        double scramble_maximum = 0;
        int scramble_counter = 0;
        double scramble_stdev = 0;
        double scramble_sterr = 0;

        //alignment signal data, name may need to be changed
        double alignment_maximum = 0;
        int alignment_counter = 0;
        double alignment_stdev = 0;
        double alignment_sterr = 0;

        //enumerate signal data, name may need to be changed
        double enumerate_maximum = 0;
        int enumerate_counter = 0;
        double enumerate_stdev = 0;
        double enumerate_sterr = 0;
	};

	// DNase DNase_data;
	PWM PWM_data;
	TFMdata TFM_data;
	accumSummaryData accumSummary_data;
	SettingsForSNPEffectMatrix settings;
	MaximumAverageSignalData signal_Data;

  // name of original command passed in
	std::string command = "";

  // name of transcription factor
	std::string TF_name = "";

  // name of various files and directories
	std::string PWM_file = "";
	std::string bigwig_file = "";
	std::string DNase_file = "";
	std::string output_dir = "";
	std::string cachefile = "";

    std::map<std::string, int> kmerHash;
    std::vector<std::string> scramble_kmers;

    std::vector<std::string> filterDNaseWrapper_output;

    std::vector<std::string> signal_cache;
    std::vector<std::string> signal_cache_scramble;
    std::vector<std::string> signal_cache_enumerate;

    std::vector<std::string> signal_output;
    std::vector<std::string> signal_scramble_output;
    std::vector<std::string> signal_enumerate_output;

#ifdef DEBUG
    size_t size_of_kmerHash = 0;
#endif
};

//Declare functions in header to be used by other functions

//main files
void generateSNPEffectMatrix(Dataset &data);

//src files
void accumSummary_scale(Dataset &data, const std::string &hfile,
                        const std::string &cfile, int scale,
                        Dataset::accumSummaryData::accumSummary_dest dest);
// check accumSummary_scale calls in steps before find_signal
void alignToGenomeWrapper(Dataset &data, int iteration,
                            std::string genome);
void bowtie_genome_map(int length, const std::string& genome,
                        const std::string& file, const std::string& final_output);
void changeBase(const Dataset &data, int position,
                            std::string nucleotide,
                            std::vector<std::string> &new_kmer_vec,
                            std::string const &genome);
void checkCache(Dataset &data, const std::vector<std::string> &in_file,
                std::vector<std::string> &out_cache, const std::string &cachefile,
                Dataset::accumSummaryData::accumSummary_dest dest);
void combineBedFiles(Dataset &data);
void Enumerate_kmer(Dataset &data);
void filterDNaseWrapper(Dataset &data);
void findMaximumAverageSignalWrapper(Dataset &data,
                                     Dataset::accumSummaryData::accumSummary_dest dest);
void generatePWMfromSEM(Dataset &data);
void generateRmeplot(Dataset &data);
void generateRplot(Dataset &data);
void generateSelfInfo(Dataset &data);
void generateSEM(const Dataset &data);
void generateSignalMethylTable(Dataset &data);
//                                      as noted in original implementation
double get_threshold(Dataset &data, double pval);
void pwm_to_tfm(Dataset &data);
void quality_control(Dataset &data);
void scramble_kmer(Dataset &data);
bool seq_col_to_fa(const std::vector<std::string> &column,
                    const std::string &file);
void writeCache(Dataset &data, const std::string &cache,
                Dataset::accumSummaryData::accumSummary_dest dest);
//Library converted function
void bedtools(Dataset &data);

#endif /* iterativeSEM_hpp */
