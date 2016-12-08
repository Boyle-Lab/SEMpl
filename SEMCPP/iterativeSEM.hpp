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

/*
 example execution from command line
 "./iterativeSEM.pl -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output examples/HNF4A/"
*/

// overview
// a struct to contain an instance of PWM data and DNase data
// data members made public for ease of access, otherwise wouldn't given more time

struct Dataset {
	// overview
	// a struct to contain and manage the PWM data as given in the example file
	struct PWM{

	    static const int NUM_ROWS = 13, NUM_COLUMNS = 5;


	    // holds first three inputs as given in example
	    std::string first_input, second_input, third_input;
	    // holds the characters at the end of each matrix row, holds the two characters at the end
	    std::string end_of_line_char, end_of_matrix_char;
	    // holds the integer values of the matrix
	    std::array<std::array<int, NUM_COLUMNS>, NUM_ROWS> matrix_arr;
	    // holds the modified SEM version of the PWM
	    std::array<std::array<int, NUM_COLUMNS>, NUM_ROWS> sem_arr;

	};
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
	struct TFMdata{
		// a c g t
		static const int LETTER_NUM = 4;
		// first letter is a, then c, then g, then t, at least in example
		std::array<std::vector<int>, LETTER_NUM> letter_array;
    // should be int or char? int for now 
	};
	struct accumSummaryData{
		// lines of output from accumSummary_scale.pl
		std::vector<std::string> accum_lines;
		// max of output from accumSummary_scale.pl
		std::vector<double> accum_max;
	};
	struct SettingsForSNPEffectMatrix{
		bool delSNPList = true, delAlignmentBed = true, delFilteredBed = true;
		bool delSignalFile = false, writecache = false, fastrun = false, verbose = false;
		int iteration = -1;
		double threshold;
    bool debug = false;
	};


	DNase DNase_data;
	PWM PWM_data;
	TFMdata TFM_data;
	accumSummaryData accumSummary_data;
	SettingsForSNPEffectMatrix settings;

  // name of original command passed in
	std::string command = "";

  // name of transcription factor
	std::string TF_name = "";

  // name of various files and directories
	std::string PWM_file = "";
	std::string bigwig_file = "";
	std::string DNase_file = "";
	std::string output_dir = "";
	std::string cache_file = "";
};

//Declare functions in header to be used by other functions

//main files
void generateSNPEffectMatrix(Dataset &data);

//src files
void accumSummary_scale(Dataset &data);
void alignToGenomeWrapper(Dataset &data);
void bowtie_genome_map(Dataset &data);
void changeBase(Dataset &data);
void checkCache(Dataset &data);
void combineBedFiles(Dataset &data);
void common(Dataset &data);
void Enumerate_kmer(Dataset &data);
void filterDNaseWrapper(Dataset &data);
void findMaximumAverageSignalWrapper(Dataset &data);
void generatePWMfromSEM(Dataset &data);
void generateRmeplot(Dataset &data);
void generateRplot(Dataset &data);
void generateSelfInfo(Dataset &data);
void generateSEM(Dataset &data);
void generateSignalMethylTable(Dataset &data);
double get_threshold(Dataset &data);
void pwm_to_tfm(Dataset &data);
void quality_control(Dataset &data);
void scramble_kmer(Dataset &data);
void seq_col_to_fa(Dataset &data);
void writeCache(Dataset &data);


#endif /* iterativeSEM_hpp */
