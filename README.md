# SEM_CPP
C++ implementation of the SEM algorithm

#Compilation
Copy entire repository to computer, then cd into lib/libBigWig-master/ and run "make".

cd back into directory where makefile is present, or from the previous location "cd .." twice.

Run "make" again.

#Tasks

Continue C++ implementation of algorithm.
Use single header file or multiple, but plan to equate
each .cpp file with each .pl file
Use descriptive variable names
	
ALWAYS INITIALIZE VARIABLES
	can be quite difficult to debug uninitialized variables
Use compiler flags as provided in Makefile
Will use struct as opposed to class for simplicity
Will use public data members for simplicity of struct

Use assert(bool) to check REQUIRES clauses and for error checking
	
More comments
	
#Notes

12-1-2016
Appropriate interface found, will need to eventually test.

Also found appropriate pv2sc library, sucessfully compiles, need to test

11-16-2016
In process of finding a C++/C interface for BigWig files/software, as opposed to attempting to rewrite perl module,

Below notes are now irrelevant

#Resources

http://search.cpan.org/~lds/Bio-BigFile/lib/Bio/DB/BigWig.pm

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922891/

For tfmpvalue binary source code, use below

tfmpvalue.cpp Bonsai

Possible BigWig file/software interface below

https://github.com/deltadev/bbi

github.com/jayhesselberth/libBigWig
