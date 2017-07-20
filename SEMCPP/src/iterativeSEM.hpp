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
#include <exception>

#include <iostream>



#define NAN_VALUE -28.0

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

	    static const int NUM_ROWS = 13, NUM_COLUMNS = 4;

	    // holds the integer values of the matrix
	    std::array<std::array<int, NUM_COLUMNS>, NUM_ROWS> matrix_arr;
	    // holds the modified SEM version of the PWM
	    // std::array<std::array<int, NUM_COLUMNS>, NUM_ROWS> sem_arr;

	};
	struct TFMdata{
		// a c g t
		static const int LETTER_NUM = 4;
		// first letter is a, then c, then g, then t, at least in example
		std::array<std::vector<int>, LETTER_NUM> letter_array;
    // should be int or char? int for now
	};
	struct accumSummary_type{
	    //Alignment summary data
        enum class accumSummary_dest{alignment, scrambled, enumerated, none};
		// lines of output from accumSummary_scale.pl
        // making change to support storing all accum_summary_data
		// std::vector<std::string> align_accum_lines;
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
		double threshold = -1.0;
        // negative threshold value indicates not defined
	};
	// contains data from findMaximumAverageSignalWrapper
	struct MaximumAverageSignalData{
        
        double scramble_maximum = 0.0;
        int scramble_counter = 0;
        double scramble_stdev = 0.0;
        double scramble_sterr = 0.0;

        // Maxiumum average signal data will be built iteratively
        // as in all the calculations of average... will not occur
        // at once, but once the data is available to be calculated
        double alignment_maximum = 0.0;
        int alignment_counter = 0;
        double alignment_stdev = 0.0;
        double alignment_sterr = 0.0;

        
        double enumerate_maximum = 0.0;
        int enumerate_counter = 0;
        double enumerate_stdev = 0.0;
        double enumerate_sterr = 0.0;
	};

    std::map< std::pair<int, char>, double> sig_deets_maximum;
    std::map< std::pair<int, char>, int> sig_deets_counter;
    std::map< std::pair<int, char>, double> sig_deets_stdev;
    std::map< std::pair<int, char>, double> sig_deets_sterr;

	// DNase DNase_data;
	PWM PWM_data;
	TFMdata TFM_data;
	accumSummary_type accumSummary_data;
	SettingsForSNPEffectMatrix settings;
	MaximumAverageSignalData Signal_data;

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

    std::map<std::string, double> kmerHash;

    // std::vector<std::pair<std::pair<int, char>, 
    //             std::vector<std::string> > > signal_cache;
    std::map< std::pair<int, char>, std::vector<std::string> > signal_cache;

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
                        Dataset::accumSummary_type::accumSummary_dest dest);
// check accumSummary_scale calls in steps before find_signal
void alignToGenomeWrapper(Dataset &data, int iteration,
                            std::string genome);
void bowtie_genome_map(int length, const std::string& genome,
                        const std::string& file, const std::string& final_output,
                        bool verbose);
void changeBase(const Dataset &data, int position, const char nucleotide, 
                std::vector<std::string> &new_kmer_vec, 
                const std::string &genome);
void checkCache(Dataset &data, const std::vector<std::string> &in_file,
                std::vector<std::string> &out_cache, const std::string &cachefile,
                Dataset::accumSummary_type::accumSummary_dest dest,
                int position = -1, char bp = 'Q');
void combineBedFiles(Dataset &data);
void Enumerate_kmer(Dataset &data);
void filterDNaseWrapper(const Dataset &data);
void findMaximumAverageSignalWrapper(Dataset &data,
                                     Dataset::accumSummary_type::accumSummary_dest dest);
void generatePWMfromSEM(const Dataset &data);
void generateRmeplot(Dataset &data);
void generateRplot(const Dataset &data);
void generateSelfInfo(Dataset &data);
void generateSEM(const Dataset &data);
void generateSignalMethylTable(Dataset &data);
//                                      as noted in original implementation
double get_threshold(Dataset &data, double pval);
void pwm_to_tfm(Dataset &data);
void quality_control(const Dataset &data);
void scramble_kmer(Dataset &data);
bool seq_col_to_fa(const std::vector<std::string> &column,
                    const std::string &file);
void writeCache(Dataset &data, const std::string &cache,
                Dataset::accumSummary_type::accumSummary_dest dest);

std::string read_pwm(Dataset &data, std::string file);

#endif /* iterativeSEM_hpp */
