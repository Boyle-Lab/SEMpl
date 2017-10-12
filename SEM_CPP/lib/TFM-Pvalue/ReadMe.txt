TFMpvalue
Copyright 2007 CRIStAL(LIFL)-USTL-INRIA. All rights reserved.
Contact: tfmpvalue@lifl.fr


DESCRIPTION

TFMpvalue is a software suite to perform exact P-value and score threshold 
computation onto Position Weight Matrices (PWMs for short).

Four programs are available:

- TFMpvalue-pv2sc takes as input a matrix and a pvalue and outputs the
  score threshold associated to the Pvalue
  
- TFMpvalue-sc2pv takes as input a matrix and a score and outputs the
  P-value corresponding to the score
  
- TFMpvalue-distrib takes as input a matrix and two scores min and max 
  and computes the distribution of scores between the bounds
  
- TFMpvalue-fastpvalue takes as input a matrix and a score and outputs the
  P-value corresponding to the score using the FastPvalue algorithm for a given
  granularity (rounding)
  
All programs take as input the background probabilities for each
letter (the sum must equal 1) and a matrix. The matrix contains real
values and can be either a count matrix, a frequency matrix or a
position weight matrix (ie log ratio matrix). By default, the program
assumes that the matrix is a simple count matrix and transforms it into
a log-ratio matrix including a pseudocount.

The programs TFMpavlue-sc2pv and TFMpvalue-pv2sc output a line with the matrix 
length, the score threshold, the P-value and the granularity separated by a 
space. TFMpvalue-distrib outputs for each accessible score between min and max
the score, the P-value and the probability of such a score separated by a space.


BUILDING

  Simply run 'make' to compile the programs. You can activate the VERBOSE
  option by adding -DVERBOSE=1 to the CFLAGS line in the Makefile.
  
  You can test the software by running 'make testpv2sc', 'make testsc2pv' and 'make
  testdistrib'. These commmands respectively compute the score associated to a P-value of 
  1e-5, the P-value associated to a score of 8.77 and the score distribution 
  between 8 and 9 for the MA0045.pfm matrix file (from JASPAR database : 
  http://jaspar.genereg.net/)

SYNOPSIS

  TFMpvalue-pv2sc -a X -t X -g X -c X -m matrix_filename  -p pvalue [-w]
  TFMpvalue-sc2pv sc2pv -a X -t X -g X -c X -m matrix_filename -s threhold [-w]
  TFMpvalue-distrib -a X -t X -g X -c X -m matrix_filename -s min_score -S max_score -G granularity

  -a -t -c -g : background probabilies, a number between 0 and 1. The sum should be 1.

  -m : matrix file which is assumed to be in the following format :
  
sc1A sc2A sc3A ...
sc1C sc2C sc3C ...
sc1G sc2G sc3G ...
sc1T sc2T sc3T ...

  -w : the matrix is already a weight matrix, otherwise it is assumed to be a count matrix

  -p : requested pvalue

  -s : score threshold (minimum score in case of TFMpvalue-distrib)
  
  -S : maximum score (only for TFMpvalue-distrib)
  
  -G : granularity. Is used to round the matrix, a number less than 1.
  
  
REFERENCE

Efficient and accurate P-value computation for Position Weight Matrices.
Algorithms Mol Biol. 2007 Dec 11;2:15.
Touzet H, Varré JS.


