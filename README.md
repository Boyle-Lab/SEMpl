# SEM_CPP
C++ implementation of the SEM algorithm

# System Requirements

## Hardware Requirements
Generation of a SEM requires variable RAM and disk storage based on the size of the initial PWM being considered. For minimal performance, we recommend a computer with the following specs:

RAM: 64+ GB  
CPU: 8+ cores, 3.4+ GHz/core

The runtime on this minimal system is approximately XX CPU hours. Compile time is approximately 35 seconds.

## Software Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Ubuntu 18.04  
Packages: libcurl4-dev

## Demo

We include a small of generation of the SEM for HNF4A in HepG2 cells. Execution time of this demo is approximately XX seconds on 20 threads. The expected output is:
```
Running Iterative SEM building..
        PWM: examples/MA0114.1.pwm
        merge_file: examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak.gz
        bigwig: examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig
        TF_name: HNF4A
         output: results/HNF4A/
        cachefile flag: results/HNF4A/HNF4A.cache.db
        verbose
....
```

# Installation
Clone a copy of the SEMpl repository and submodules:

```
git clone --recurse-submodules git@github.com:Boyle-Lab/SEM_CPP.git
```

Build external libraries (this should be shored up)
```
cd SEM_CPP/SEM_CPP/lib/libBigWig-master
make
cd ../TFM-Pvalue
make SEMCPPobj
cd ../bowtie-1.0.0
make
cd ..
mv */*.so ../
cd ..
```

Symlink to bowtie index location
```
ln -s /data/genomes/hg19/bowtie_index/ data
```

Build SEMpl
```
make
```
 
# Usage information
SEMpl runs as an iterative process and requires specific input files (need more details). The following example will build the SEM for HNF4a in HepG2 cells given the example data
```
./iterativeSEM -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output results/HNF4A
```

# Testing
Run "make test" to compile and run this input example.

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
