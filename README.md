# SEM_CPP
C++ implementation of the SEM algorithm

# Expectations
1. bowtie binary in ./bin
2. bedtools binary in ./bin
3. argument to "-readcache" is not optional.

# Compilation
1. Clone to computer, then cd into lib/libBigWig-master and run "make".
2. cd into lib/TFM-Pvalue and run "make".
3. cd into lib/bowtie-1.0.0 and run "make".
4. mv all .so files into lib/
5. There should be three .so files.
6. mkdir obj for object files.
7. mkdir bin and place expected binaries.
8. run "export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:/home/cmorteru/SEM_CPP/SEMCPP/lib"" without the outermost quotes so that dynamic libraries can be found.
9. Run "make" without quotes to compile everything.

# Running
	./iterativeSEM -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output results/HNF4A

# Testing
Run "make test" to compile and run on example input, assuming .so libraries are in the correct directory.

# Tasks

1. Continue C++ implementation of algorithm.
2. Use descriptive variable names
3. ALWAYS INITIALIZE VARIABLES
   * can be quite difficult to debug uninitialized variables
4. Use assert(bool) to check REQUIRES clauses and for error checking

# Resources

1. http://search.cpan.org/~lds/Bio-BigFile/lib/Bio/DB/BigWig.pm
2. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922891/
3. For tfmpvalue binary source code, use below
   * tfmpvalue.cpp Bonsai
4. https://github.com/deltadev/bbi
5. github.com/jayhesselberth/libBigWig
5. https://github.com/BenLangmead/bowtie/tree/85fad1d77775e773611c3d16f3d1fb527c32f4b7
6. sqlite.org
