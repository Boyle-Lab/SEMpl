#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
//#include <GetOpt.h>
#include <cstdlib>
#include <stdio.h>

#include "Matrix.h"
#include "ArgumentException.h"
#include "TFMpvalue.h"

#include <Rcpp.h>
using namespace std;

/********************************************************************
 * Free the Matric class
 * *****************************************************************/
void freeMatrix(Matrix m, int nrow){
  // free the memory allocated, not typical Rcpp way
  for(int i=0; i<nrow; i++){
    delete[] m.mat[i];
    delete[] m.matInt[i];
  }
  delete[] m.matInt;
  delete[] m.mat;
  delete[] m.offsets;
  delete[] m.minScoreColumn;
  delete[] m.maxScoreColumn;
  delete[] m.sum;
  delete[] m.bestScore;
  delete[] m.worstScore;
}

/********************************************************************
 * .Call() Entry points  sc2pv
 * *****************************************************************/
RcppExport SEXP sc2pv (SEXP mat, SEXP Rscore, SEXP bg, SEXP type){
  
  Rcpp::NumericVector background(bg);
  Rcpp::NumericMatrix matrix(mat);
  Rcpp::NumericVector RScoreVec(Rscore);
  Rcpp::CharacterVector Type(type);
  // Fill with background
  Matrix m(background[0], background[1], background[2], background[3]);
  // Fill with matrix
  int i=0, j=0;
  m.mat = new double*[4];
  int ncol = matrix.ncol();
  int nrow = matrix.nrow();
  m.length = ncol;
  for(i=0; i<nrow; i++){
    m.mat[i] = new double[ncol];
    for(j=0; j<ncol; j++){
      m.mat[i][j] = matrix[j*nrow+i];
    }
  }
  /*cout << "INITIAL MATRIX" << endl;
  for(j=0; j<ncol; j++){
    for(i=0; i<nrow; i++){
      cout << m.mat[i][j] << "\t";
    }
    cout << endl;
  }*/
  //cout << "Matrix length  : " << m.length << endl;
  
  // toPWM when it is PFM
  if(!strcmp(Type[0], "PFM")){
    m.toLogOddRatio();
  }
  
  // testScoreToPvalue
  double initialGranularity = 0.1;
  bool forcedGranularity = false;
  double maxGranularity = 1e-9;
  double requestedScore = RScoreVec[0];
  qlonglong totalSize = 0;
  qlonglong totalOp = 0;
  qlonglong max;
  qlonglong min;
  double ppv;
  double pv;
  qlonglong score;
  for (double granularity = initialGranularity; granularity >= maxGranularity; granularity /= 10){
    //cout << "Computing rounded matrix with granularity " << granularity << endl;
    m.computesIntegerMatrix(granularity);
    max = requestedScore*m.granularity + m.offset + m.errorMax+1;
    min = requestedScore*m.granularity + m.offset - m.errorMax-1;
    score = requestedScore*m.granularity + m.offset;
    /*cout << "Score range : " << m.scoreRange << endl;
    cout << "Min         : " << min << endl;
    cout << "Max         : " << max << endl;
    cout << "Precision   : " << m.granularity << endl;
    cout << "Error max   : " << m.errorMax << endl;
    cout << "Computing pvalue for requested score " << requestedScore << " " << score << endl;*/
#ifdef MEMORYCOUNT
    m.totalMapSize = 0;
    m.totalOp = 0;
#endif
    m.lookForPvalue(score,min,max,&ppv,&pv);
 /*   cout << "Prev. Pvalue  : " << ppv << endl;
    cout << "Pvaluex       : " << pv << endl;
    cout << "Comp. score   : " << score << endl;*/
#ifdef MEMORYCOUNT
    totalSize += m.totalMapSize;
    totalOp += m.totalOp;
#endif
    //cout << "***********************************************" << endl;
    if (ppv == pv) {
      if (!forcedGranularity) {
        break;
      }
    }
  }
  Rcpp::NumericVector ans(1);
  ans[0] = pv;
  // free the memory allocated, not typical Rcpp way
  //for(i=0; i<nrow; i++){
  //  delete[] m.mat[i];
  //}
  //delete[] m.mat;
  freeMatrix(m, nrow);

  return Rcpp::wrap(ans);
}

/********************************************************************
 * .Call() Entry points pv2sc
 * *****************************************************************/
RcppExport SEXP pv2sc (SEXP mat, SEXP Rpvalue, SEXP bg, SEXP type){
  Rcpp::NumericVector background(bg);
  Rcpp::NumericMatrix matrix(mat);
  Rcpp::NumericVector RPvalueVec(Rpvalue);
  Rcpp::CharacterVector Type(type);
  // Fill with background
  Matrix m(background[0], background[1], background[2], background[3]);
  // Fill with matrix
  int i=0, j=0;
  m.mat = new double*[4];
  int ncol = matrix.ncol();
  int nrow = matrix.nrow();
  m.length = ncol;
  for(i=0; i<nrow; i++){
    m.mat[i] = new double[ncol];
    for(j=0; j<ncol; j++){
      m.mat[i][j] = matrix[j*nrow+i];
    }
  }
  if(!strcmp(Type[0], "PFM")){
    m.toLogOddRatio();
  }
  // PvalueToScore
  double initialGranularity = 0.1;
  bool forcedGranularity = false;
  double requestedPvalue = RPvalueVec[0];
  double maxGranularity = 1e-10;
  bool sortColumns = false;
  qlonglong decrgr = 10;

  m.computesIntegerMatrix(initialGranularity);
  qlonglong max = m.maxScore+ceil(m.errorMax+0.5);
  qlonglong min = m.minScore;
  double pv;
  qlonglong score;
  for (double granularity = initialGranularity; granularity >= maxGranularity; granularity /= decrgr) {
    m.computesIntegerMatrix(granularity);
    double ppv;
    score = m.lookForScore(min,max,requestedPvalue,&pv,&ppv);
    min = (score - ceil(m.errorMax+0.5)) * decrgr;
    max = (score + ceil(m.errorMax+0.5)) * decrgr;
    if (pv == ppv) {
      if (!forcedGranularity) {
        break;
      }
    }
  }
  Rcpp::NumericVector ans(1);
  ans[0] = ((score-m.offset)/m.granularity);
  freeMatrix(m, nrow);
  return Rcpp::wrap(ans);
}

/********************************************************************
 * .Call() Entry point FastPvalue
 * *****************************************************************/
RcppExport SEXP FastPvalue(SEXP mat, SEXP Rscore, SEXP bg, SEXP type, 
    SEXP Rgranularity){
  Rcpp::NumericVector background(bg);
  Rcpp::NumericMatrix matrix(mat);
  Rcpp::NumericVector RscoreVec(Rscore);
  Rcpp::CharacterVector Type(type);
  Rcpp::NumericVector granularityVec(Rgranularity);

  // Fill with background
  Matrix m(background[0], background[1], background[2], background[3]);
  // Fill with matrix
  int i=0, j=0;
  m.mat = new double*[4];
  int ncol = matrix.ncol();
  int nrow = matrix.nrow();
  m.length = ncol;
  for(i=0; i<nrow; i++){
    m.mat[i] = new double[ncol];
    for(j=0; j<ncol; j++){
      m.mat[i][j] = matrix[j*nrow+i];
    }
  }
  if(!strcmp(Type[0], "PFM")){
    m.toLogOddRatio();
  }  

  // Fast Pvalue
  double granularity = granularityVec[0];
  double score = RscoreVec[0];
  m.computesIntegerMatrix(granularity,true);
  double pvalue = m.fastPvalue(&m,(qlonglong)(score * m.granularity + m.offset));
  Rcpp::NumericVector ans(1);
  ans[0] = pvalue;
  freeMatrix(m, nrow);
  return Rcpp::wrap(ans);
}

/********************************************************************
 * .Call() Entry point LAZY
 * *****************************************************************/
RcppExport SEXP lazyScore(SEXP mat, SEXP Rpvalue, SEXP bg, SEXP type,
    SEXP Rgranularity){
  Rcpp::NumericVector background(bg);
  Rcpp::NumericMatrix matrix(mat);
  Rcpp::NumericVector RpvalueVec(Rpvalue);
  Rcpp::CharacterVector Type(type);
  Rcpp::NumericVector granularityVec(Rgranularity);

  // Fill with background
  Matrix m(background[0], background[1], background[2], background[3]);
  // Fill with matrix
  int i=0, j=0;
  m.mat = new double*[4];
  int ncol = matrix.ncol();
  int nrow = matrix.nrow();
  m.length = ncol;
  for(i=0; i<nrow; i++){
    m.mat[i] = new double[ncol];
    for(j=0; j<ncol; j++){
      m.mat[i][j] = matrix[j*nrow+i];
    }
  }
  if(!strcmp(Type[0], "PFM")){
    m.toLogOddRatio();
  }
  
  // LAZY
  double granularity = granularityVec[0];
  double requestedPvalue = RpvalueVec[0];
  m.computesIntegerMatrix(granularity,true);
  map<qlonglong, double> *nbOcc = new map<qlonglong, double> [m.length+1];
  map<qlonglong, double> *pbuf = new map<qlonglong, double> [m.length+1];
  qlonglong score = m.maxScore+ceil(m.errorMax);
  qlonglong d = 0;
  double pv = 0;
  nbOcc[m.length][score] = pv;
  while (pv <= requestedPvalue) {
    score --;
    pv += _beckstette(m,&nbOcc,&pbuf,m.length-1,score,d);
    nbOcc[m.length][score] = pv;
    d++;
  }
  pv = nbOcc[m.length][score];

  Rcpp::NumericVector ans(1);
  ans[0] = ((score-m.offset)/m.granularity);
  freeMatrix(m, nrow);
  delete[] nbOcc;
  delete[] pbuf;

  return Rcpp::wrap(ans);
}

