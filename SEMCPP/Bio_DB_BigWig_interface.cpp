

// $hfile is BigWig file containing sequence reads
// my $wig = Bio::DB::BigWig->new(-bigwig, $hfile)

// Should only translated as far as to make the function work as used in original file

#include <string>
#include <ifstream>
#include <cassert>
using namespace std;

//
unknown_type new_function(string flag, string hfile){
	assert(flag.at(0) == '-');
	string bw_path = hfile;

	ifstream input(hfile);
	if(!input){
		cout << "BigWig file unable to open" << endl;
		exit(EXIT_FAILURE);
	}
	else{
		input.close();
	}

	
}
		
