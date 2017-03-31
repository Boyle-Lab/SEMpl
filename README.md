# SEM_CPP
C++ implementation of the SEM algorithm

# Compilation
Clone to computer, then cd into lib/libBigWig-master and run "make".

cd into lib/TFM-Pvalue and run "make".

cd into lib/bowtie-1.0.0 and run "make".

mv all .so files into lib/\n
There should be three .so files.

mkdir obj for executable files.

Run "make" again.

# Tasks

Continue C++ implementation of algorithm.
Use descriptive variable names
	
ALWAYS INITIALIZE VARIABLES
	can be quite difficult to debug uninitialized variables

Use assert(bool) to check REQUIRES clauses and for error checking

# Resources

http://search.cpan.org/~lds/Bio-BigFile/lib/Bio/DB/BigWig.pm

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922891/

For tfmpvalue binary source code, use below

tfmpvalue.cpp Bonsai

https://github.com/deltadev/bbi

github.com/jayhesselberth/libBigWig

https://github.com/BenLangmead/bowtie/tree/85fad1d77775e773611c3d16f3d1fb527c32f4b7

sqlite.org
