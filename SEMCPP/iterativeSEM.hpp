//
//  iterativeSEM.hpp
//  SEMCPP
//
//  Created by Cody Morterud on 11/9/16.
//  Copyright © 2016 Boyle Lab. All rights reserved.
//

#ifndef iterativeSEM_hpp
#define iterativeSEM_hpp

#include <string>
#include <array>
#include <fstream>
#include <vector>

/* 
 example execution from command line 
 "./iterativeSEM.pl -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output examples/HNF4A/"
*/

// overview
// a struct to contain an instance of PWM data and DNase data
// data members made public for ease of access, otherwise wouldn't with more time

struct Dataset {


	// overview
	// a struct to contain and manage the PWM data as given in the example file
	struct PWM{
	
	    static const int MATRIX_SIZE = 13, ROW_SIZE = 5;
	    
	
	    // holds first three inputs as given in example
	    std::string first_input, second_input, third_input;
	    // holds the characters at the end of each matrix row, holds the two characters at the end
	    std::string end_of_line_char, end_of_matrix_char;
	    // holds the integer values of the matrix
	    std::array<std::array<int, ROW_SIZE>, MATRIX_SIZE> matrix_arr;
	    
	
	    

	    // Modifies: all member variables
	    // Effects: initializes all member variables to 0 or nullptr
	    PWM();
	    
	    // Requires: input is a valid input stream with correctly formatted data
	    // Modifies: all member variables, input
	    // Effects: initializes member variables as given in the example .pwm file
	    //          in the input ifstream object
	    PWM(std::ifstream &input);
	    
	    // Requires: other is completely defined
	    // Modifies: all member variables
	    // Effects: copy constructor to copy over all data from other
	    //          uses deep copy
	    PWM(const PWM &other);
	    
	    // Requires: other is completely defined
	    // Modifies: all member variables
	    // Effects: assignment to copy over all data from other
	    //          uses deep copy
	    PWM operator=(const PWM &other);
	    
	    // Effects: does nothing, as member variables are objects with destructors
	    ~PWM();
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
	    
		
	    
	    // Modifies: all member variables
	    // Effects: initializes all member variables to 0 or nullptr
	    DNase();
	    
	    // Requires: input is a valid input stream with correctly formatted data in file
	    // Modifies: all member variables, input
	    // Effects: initializes member variables as given in the example .narrowpeak.gz file
	    //          in the input ifstream object
	    DNase(std::ifstream &input);
	    
	    // Requires: other is completely defined
	    // Modifies: all member variables
	    // Effects: copy constructor to copy over all data from other
	    //          uses deep copy
	    DNase(const DNase &other);
	    
	    // Requires: other is completely defined
	    // Modifies: all member variables
	    // Effects: assignment to copy over all data from other
	    //          uses deep copy
	    DNase operator=(const DNase &other);
	    
	    // Effects: does nothing, as member variables are objects without dynamic memory
	    ~DNase();
	};
	struct TFMdata{
		// a c g t
		static const int LETTER_NUM = 4;
		// first letter is a, then c, then g, then t
		std::array<std::vector<char>, LETTER_NUM> letter_array;
	};

	DNase DNase_data;
	PWM PWM_data;
	TFMdata TFM_data;
	
	// Requires: input is a valid ifstream object
	// Modifies: DNase_data, PWM_data
	// Effects: initializes PWM and DNase with their corresponding constructor taking an ifstream object
	Dataset(std::ifstream& input);

	

};

#endif /* iterativeSEM_hpp */
