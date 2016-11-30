/*
 *  Matrix.h
 *  pvalue
 *
 *  Created by Jean-Stéphane Varré on 02/07/07.
 *  Copyright 2007 LIFL-USTL-INRIA. All rights reserved.
 *
 */

#ifndef __MATRIX__
#define __MATRIX__

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <cstdlib>

#include "FileException.h"
#include "ParseException.h"

using namespace std;

#define ROUND_TO_INT(n) ((long long)floor(n))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

extern map<char,int> OPTIONS;

//#define PRINTLOGRATIO

class Matrix {
  
private:
  
  
  /**
  * Split a string following delimiters
   */
  void tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
    
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  }
  
  
public:
  

  // used for efficiency tests
  long long totalMapSize;
  long long totalOp;
  
  double ** mat; // the matrix as it is stored in the matrix file
  int length;
  double granularity; // the real granularity used, greater than 1
  long long ** matInt; // the discrete matrix with offset
  double errorMax;
  long long *offsets; // offset of each column
  long long offset; // sum of offsets
  long long *minScoreColumn; // min discrete score at each column
  long long *maxScoreColumn; // max discrete score at each column
  long long *sum;
  long long minScore;  // min total discrete score (normally 0)
  long long maxScore;  // max total discrete score
  long long scoreRange;  // score range = max - min + 1
  long long *bestScore;
  long long *worstScore;
  double background[4];
  
  Matrix() {
    granularity = 1.0;
    offset = 0;
    background[0] = background[1] = background[2] = background[3] = 0.25;
  }
  
  Matrix(double pA, double pC, double pG, double pT) {
    granularity = 1.0;
    offset = 0;
    background[0] = pA;
    background[1] = pC;
    background[2] = pG;
    background[3] = pT;  
  }
    
  void toLogOddRatio () {
    for (int p = 0; p < length; p++) {
      double sum = mat[0][p] + mat[1][p] + mat[2][p] + mat[3][p];
      for (int k = 0; k < 4; k++) {
        mat[k][p] = log((mat[k][p] + 0.25) /(sum + 1)) - log (background[k]); 
      }
    }
#ifdef PRINTLOGRATIO
    for (int k = 0; k < 4; k++ ) {
      for (int i = 0 ; i < length; i++) {
        cout << mat[k][i] << "\t";
      }
      cout << endl;
    }
#endif
  }
  
  /**
    * Transforms the initial matrix into an integer and offseted matrix.
   */
  void computesIntegerMatrix (double granularity, bool sortColumns = true);
  
  // computes the complete score distribution between score min and max
  void showDistrib (long long min, long long max) {
    map<long long, double> *nbocc = calcDistribWithMapMinMax(min,max); 
    map<long long, double>::iterator iter;
    
    if (OPTIONS['h']) {
      cout << "Scores and p-values between " << min << " and " << max << endl;
    } 
    
    // computes p values and stores them in nbocc[length] 
    double sum = 0;
    map<long long, double>::reverse_iterator riter = nbocc[length-1].rbegin();
    while (riter != nbocc[length-1].rend()) {
      sum += riter->second;
      nbocc[length][riter->first] = sum;
      riter++;      
    }
    
    iter = nbocc[length].begin();
    while (iter != nbocc[length].end() && iter->first <= max) {
      cout << (((iter->first)-offset)/granularity) << " " << (iter->second) << " " << nbocc[length-1][iter->first] << endl;
      iter ++;
    }
  }
  
  /**
    * Computes the pvalue associated with the threshold score requestedScore.
    */
  void lookForPvalue (long long requestedScore, long long min, long long max, double *pmin, double *pmax);
    
  /**
    * Computes the score associated with the pvalue requestedPvalue.
    */
  long long lookForScore (long long min, long long max, double requestedPvalue, double *rpv, double *rppv);
    
  /** 
    * Computes the distribution of scores between score min and max as the DP algrithm proceeds 
    * but instead of using a table we use a map to avoid computations for scores that cannot be reached
    */
  map<long long, double> *calcDistribWithMapMinMax (long long min, long long max); 
    
  
  /**
    * Computes the pvalue for a given score and at a fixed granularity
   */
  long long fastPvalue (Matrix *m, long long alpha);
  
  void readJasparMatrix (string filename) {
    
    ifstream f(filename.data());
    if (!f) {
      throw new FileException();
    } 
    
    string str;
    this->length = 0;
    vector<string> v;
    mat = new double*[4];
    
    for (int j = 0; j < 4; j++) {
      getline(f,str);
      tokenize(str,v," \t|");
      this->length = v.size();
      this->mat[j] = new double[this->length];
      for (unsigned int i = 0; i < v.size(); i++) {
        mat[j][i] = atof(v.at(i).data());
      }
      v.clear();
    }
    
    f.close();
    
    
#ifdef PRINTVERBOSE
    cout << "INITIAL MATRIX" << endl;    
    for (int k = 0; k < 4; k++ ) {
      for (int i = 0 ; i < length; i++) {
        cout << mat[k][i] << "\t";
      }
      cout << endl;
    }
#endif
    
    
  }
  
  
  void readHorizontalMatrix (string filename) {
    
    ifstream f(filename.data());
    if (!f) {
      throw new FileException();
    } 
    
    string str;
    this->length = 0;
    // comment out for JASPAR matrices
    getline(f,str); // line with matrix name and family
    vector<string> v;
    mat = new double*[4];
    
    for (int j = 0; j < 4; j++) {
      getline(f,str);
      tokenize(str,v," \t|");
      this->length = v.size() -1; // not -1 for JASPAR
      this->mat[j] = new double[this->length];
      for (unsigned int i = 1; i < v.size(); i++) { // 1 if not JASPAR
        mat[j][i-1] = atof(v.at(i).data());
      }
      v.clear();
    }
    
    f.close();
    
  }
  
}; /* Matrix */




#endif
