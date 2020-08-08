# SEMpl
C++ implementation of the SEM algorithm

Sierra S Nishizaki, Natalie Ng, Shengcheng Dong, Robert S Porter, Cody Morterud, Colten Williams, Courtney Asman, Jessica A Switzenberg, Alan P Boyle, Predicting the effects of SNPs on transcription factor binding affinity, Bioinformatics, Volume 36, Issue 2, 15 January 2020, Pages 364–372, https://doi.org/10.1093/bioinformatics/btz612

We have made all of the SEMs generated as part of this work available [here](SEMs/).

# System Requirements

## Hardware Requirements
Generation of a SEM requires variable RAM and disk storage based on the size of the initial PWM being considered. For minimal performance, we recommend a computer with the following specs:

RAM: 64+ GB  
CPU: 8+ cores, 3.4+ GHz/core

The runtime on this minimal system is approximately 38 CPU hours. Compile time is approximately 35 seconds.

## Software Requirements

The package development version is tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Ubuntu 18.04  
Packages: libcurl4-dev

## Demo

We include a small of generation of the SEM for HNF4A in HepG2 cells. Execution time of this demo is approximately 6791 seconds on 20 threads. The expected output is:
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
git clone --recurse-submodules https://github.com/Boyle-Lab/SEMpl.git
```

Build external libraries:
```
cd SEMpl/lib/libBigWig
make
cd ..
make
mv */*.so .
cd ..
```

Symlink to bowtie index location (use your own index location):
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
./iterativeSEM -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -genome data/hg19 -output results/HNF4A
```

# Testing
Run "make test" to compile and run this input example.

