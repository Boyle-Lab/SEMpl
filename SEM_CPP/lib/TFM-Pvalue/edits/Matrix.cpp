/*
 *  Matrix.cpp
 *  pvalue
 *
 *  Created by Jean-Stéphane Varré on 02/07/07.
 *  Copyright 2007 LIFL-USTL-INRIA. All rights reserved.
 *
 */

#include "Matrix.h"
#include <map>
#include <memory>

//#define PRINTVERBOSE 
//#define SHOWCERR
//#define VERBOSE

void Matrix::computesIntegerMatrix (double granularity, bool sortColumns) {
  double minS = 0, maxS = 0;
  double scoreRange;
  
  // computes precision
  for (int i = 0; i < length; i++) {
    double min = mat[0][i];
    double max = min;
    for (int k = 1; k < 4; k++ )  {
      min = ((min < mat[k][i])?min:(mat[k][i]));
      max = ((max > mat[k][i])?max:(mat[k][i]));
    }
    minS += min;
    maxS += max;
  } 
  
  // score range
  scoreRange = maxS - minS + 1;
  
  if (granularity > 1.0) {
    this->granularity = granularity / scoreRange;
  } else if (granularity < 1.0) {
    this->granularity = 1.0 / granularity;
  } else {
    this->granularity = 1.0;
  }
  
  matInt = new qlonglong *[length];
  for (int k = 0; k < 4; k++ ) {
    matInt[k] = new qlonglong[length];
    for (int p = 0 ; p < length; p++) {
      matInt[k][p] = ROUND_TO_INT((double)(mat[k][p]*this->granularity)); 
    }
  }
  
#ifdef PRINTVERBOSE
  /*cout << "SCORE RANGE : " << minS << " -> " << maxS << endl;
  
  cout << "PRECISION " << this->granularity << endl;
  
  cout << "INTEGER MATRIX WITHOUT OFFSET" << endl;
  
  for (int k = 0; k < 4; k++ ) {
    for (int i = 0 ; i < length; i++) {
      cout << matInt[k][i] << "\t";
    }
    cout << endl;
  }*/
#endif
  
  this->errorMax = 0.0;
  for (int i = 1; i < length; i++) {
    double maxE = mat[0][i] * this->granularity - (matInt[0][i]);
    for (int k = 1; k < 4; k++) {
      maxE = ((maxE < mat[k][i] * this->granularity - matInt[k][i])?(mat[k][i] * this->granularity - (matInt[k][i])):(maxE));
    }
    this->errorMax += maxE;
  }
#ifdef PRINTVERBOSE
  //cout << "  ERROR MAX : " << this->errorMax << endl;
#endif
  
  if (sortColumns) {
    // sort the columns : the first column is the one with the greatest value
    qlonglong min = 0;
    for (int i = 0; i < length; i++) {
      for (int k = 0; k < 4; k++) {
        min = MIN(min,matInt[k][i]);
      }
    }
    min --;
    qlonglong *maxs = new qlonglong [length];
    for (int i = 0; i < length; i++) {
      maxs[i] = matInt[0][i];
      for (int k = 1; k < 4; k++) {
        if (maxs[i] < matInt[k][i]) {
          maxs[i] = matInt[k][i];
        }
      }
    }
    qlonglong **mattemp = new qlonglong *[4];
    for (int k = 0; k < 4; k++) {        
      mattemp[k] = new qlonglong [length];
    }
    for (int i = 0; i < length; i++) {
      qlonglong max = maxs[0];
      int p = 0;
      for (int j = 1; j < length; j++) {
        if (max < maxs[j]) {
          max = maxs[j];
          p = j;
        }
      }
      maxs[p] = min;
      for (int k = 0; k < 4; k++) {        
        mattemp[k][i] = matInt[k][p];
      }
    }
#ifdef PRINTVERBOSE
    /*cout << "INTEGER MATRIX WITHOUT OFFSET ORDERED" << endl;
    for (int k = 0; k < 4; k++)  {
      for (int i = 0; i < length; i++) {
        cout << mattemp[k][i] << "\t";
      }
      cout << endl;    
    }*/
#endif
    
    for (int k = 0; k < 4; k++)  {
      for (int i = 0; i < length; i++) {
        matInt[k][i] = mattemp[k][i];
      }
    }
    for(int i=0; i<4; i++){
      delete[] mattemp[i];
    }
    delete[] mattemp;
    delete[] maxs;
  }
  
  // computes offsets
  this->offset = 0;
  offsets = new qlonglong [length];
  for (int i = 0; i < length; i++) {
    qlonglong min = matInt[0][i];
    for (int k = 1; k < 4; k++ )  {
      min = ((min < matInt[k][i])?min:(matInt[k][i]));
    }
    offsets[i] = -min;
    for (int k = 0; k < 4; k++ )  {
      matInt[k][i] += offsets[i];  
    }
    this->offset += offsets[i];
  }
  
#ifdef PRINTVERBOSE
  //cout << "OFFSET : " << this->offset << endl;
#endif
  
#ifdef PRINTVERBOSE
  /*cout << "INTEGER MATRIX WITH OFFSET" << endl;
  for (int k = 0; k < 4; k++ )  {
    for (int i = 0; i < length; i++) {
      cout << matInt[k][i] << "\t";
    }
    cout << endl;    
  }*/
#endif
  
  
  // look for the minimum score of the matrix for each column
  minScoreColumn = new qlonglong [length];
  maxScoreColumn = new qlonglong [length];
  sum            = new qlonglong [length];
  minScore = 0;
  maxScore = 0;
  for (int i = 0; i < length; i++) {
    minScoreColumn[i] = matInt[0][i];
    maxScoreColumn[i] = matInt[0][i];
    sum[i] = 0;
    for (int k = 1; k < 4; k++ )  {
      sum[i] = sum[i] + matInt[k][i];
      if (minScoreColumn[i] > matInt[k][i]) {
        minScoreColumn[i] = matInt[k][i];
      }
      if (maxScoreColumn[i] < matInt[k][i]) {
        maxScoreColumn[i] = matInt[k][i];
      }
    }
    minScore = minScore + minScoreColumn[i];
    maxScore = maxScore + maxScoreColumn[i];
    //cout << "minScoreColumn[" << i << "] = " << minScoreColumn[i] << endl;
    //cout << "maxScoreColumn[" << i << "] = " << maxScoreColumn[i] << endl;
  }
  this->scoreRange = maxScore - minScore + 1;
  
#ifdef PRINTVERBOSE
  //cout << "SCORE RANGE : " << minScore << " - " << maxScore << " : " << this->scoreRange << endl;
#endif
  
  bestScore = new qlonglong[length];
  worstScore = new qlonglong[length];
  bestScore[length-1] = maxScore;
  worstScore[length-1] = minScore;
  for (int i = length - 2; i >= 0; i--) {
    bestScore[i]  = bestScore[i+1]  - maxScoreColumn[i+1];
    worstScore[i] = worstScore[i+1] - minScoreColumn[i+1];
  }
  
}




/**
* Computes the pvalue associated with the threshold score requestedScore.
 */
void Matrix::lookForPvalue (qlonglong requestedScore, qlonglong min, qlonglong max, double *pmin, double *pmax) {
/*  
  map<qlonglong, double> *nbocc = calcDistribWithMapMinMax(min,max); 
  map<qlonglong, double>::iterator iter;
  
#ifdef SHOWCERR
  //cerr << "  Looks for Pvalue between " << min << " and " << max << " for score " << requestedScore << endl;
#endif
  // computes p values and stores them in nbocc[length] 
  double sum = nbocc[length][max+1];
  qlonglong s = max + 1;
  map<qlonglong, double>::reverse_iterator riter = nbocc[length-1].rbegin();
  while (riter != nbocc[length-1].rend()) {
    sum += riter->second;
    if (riter->first >= requestedScore) s = riter->first;
    nbocc[length][riter->first] = sum;
    riter++;      
  }
  //cout << "   s found : " << s << endl;
  
  iter = nbocc[length].find(s);
  while (iter != nbocc[length].begin() && iter->first >= s - errorMax) {
    iter--;      
  }
  //cout << "   s - E found : " << iter->first << endl;
  
#ifdef MEMORYCOUNT
  // for tests, store the number of memory bloc necessary
  for (int pos = 0; pos <= length; pos++) {
    totalMapSize += nbocc[pos].size();
  }
#endif
  
  *pmax = nbocc[length][s];
  *pmin = iter->second;
  
  delete[] nbocc;
  */
}



/**
* Computes the score associated with the pvalue requestedPvalue.
 */
qlonglong Matrix::lookForScore (qlonglong min, qlonglong max, double requestedPvalue, double *rpv, double *rppv) {
  
  std::map<qlonglong, std::shared_ptr<double>> *nbocc = calcDistribWithMapMinMax(min,max); 
  std::map<qlonglong, std::shared_ptr<double>>::iterator iter;
#ifdef SHOWCERR
  //cerr << "  Looks for score between " << min << " and " << max << endl;
#endif
cerr << "after calc" << endl;
cin.get();
  // computes p values and stores them in nbocc[length] 
  double sum = 0.0;
  std::map<qlonglong, std::shared_ptr<double>>::reverse_iterator riter = nbocc[length-1].rbegin();
  qlonglong alpha = riter->first+1;
  qlonglong alpha_E = alpha;
  nbocc[length][alpha] = std::make_shared<double>(0.0);
  while (riter != nbocc[length-1].rend()) {
    sum += *(riter->second);
    //cout << "Pv(S) " << riter->first << " " << sum << " " << requestedPvalue << endl;
    nbocc[length][riter->first] = std::make_shared<double>(sum);
    if (sum >= requestedPvalue) { 
      break;
    }
    riter++;      
  }
  //cout << "BREAK Pv(S) " << riter->first << " " << sum << " " << requestedPvalue << endl;
  if (sum > requestedPvalue) {
    alpha_E = riter->first;
    riter--;
    alpha = riter->first; 
  } else {
    if (riter == nbocc[length-1].rend()) { // path following the remark of the mail
      riter--;
      alpha = alpha_E = riter->first;
    } else {
      alpha = riter->first;
      riter++;
      sum += *(riter->second);
      alpha_E = riter->first;
    }
    nbocc[length][alpha_E] = std::make_shared<double>(sum);  
    //cout << "Pv(S) " << riter->first << " " << sum << endl;   
  }
#ifdef VERBOSE
  //cerr << riter->first << "      ALPHA found at score " << alpha << " and P-value " << nbocc[length][alpha] << endl;
  //cerr << riter->first << "      ALPHA-E found at score " << alpha_E << " and P-value " << nbocc[length][alpha_E] << endl;
#endif    
  
  // affichage des pvaleurs
  /*iter = nbocc[length].begin();
  while (iter != nbocc[length].end()) {
    cerr << iter->first << "[" << iter->second << "]" << endl;
    iter++;
  }*/
  
#ifdef MEMORYCOUNT
  // for tests, store the number of memory bloc necessary
  for (int pos = 0; pos <= length; pos++) {
    totalMapSize += nbocc[pos].size();
  }
#endif
  
  if (alpha - alpha_E > errorMax) alpha_E = alpha;
  cout << *nbocc[length][alpha] << " " << *nbocc[length][alpha_E] << flush;
  *rpv = *nbocc[length][alpha];
  *rppv = *nbocc[length][alpha_E];   

  
  // computes q values for scores greater or equal than min
  for (int pos = 0; pos <= length; pos++) {
	nbocc[pos].erase(nbocc[pos].begin(), nbocc[pos].end());
    cerr << "        map size for " << pos << " " << nbocc[pos].size() << endl;
  }

  delete [] nbocc;
  return alpha;
  
}


// computes the distribution of scores between score min and max as the DP algrithm proceeds 
// but instead of using a table we use a map to avoid computations for scores that cannot be reached
std::map<qlonglong, std::shared_ptr<double>> *Matrix::calcDistribWithMapMinMax (qlonglong min, qlonglong max) { 
  
  // maps for each step of the computation
  // nbocc[length] stores the pvalue
  // nbocc[pos] for pos < length stores the qvalue
  std::map<qlonglong, std::shared_ptr<double>> *nbocc = new std::map<qlonglong, std::shared_ptr<double>> [length+1];
  std::map<qlonglong, std::shared_ptr<double>>::iterator iter;
  
  qlonglong *maxs = new qlonglong[length+1]; // @ pos i maximum score reachable with the suffix matrix from i to length-1
  
#ifdef VERBOSE    
  //cerr << "  Calc distrib between " << min << " and " << max << endl;
#endif

  maxs[length] = 0;
  for (int i = length-1; i >= 0; i--) {
    maxs[i] = maxs[i+1] + maxScoreColumn[i];
  }
  
  // initializes the map at position 0
  for (int k = 0; k < 4; k++) {
    if (matInt[k][0]+maxs[1] >= min) {
      if(nbocc[0][matInt[k][0]] != nullptr) {
        nbocc[0][matInt[k][0]] = std::make_shared<double>(*nbocc[0][matInt[k][0]] + background[k]);
      } else {
        nbocc[0][matInt[k][0]] = std::make_shared<double>(background[k]);
      }
    }
  }
cerr << "!" << endl;
cin.get();
  
  // computes q values for scores greater or equal than min
  nbocc[length-1][max+1] = std::make_shared<double>(0.0);
  for (int pos = 1; pos < length; pos++) {
    iter = nbocc[pos-1].begin();
    while (iter != nbocc[pos-1].end()) {
      for (int k = 0; k < 4; k++) {
        qlonglong sc = iter->first + matInt[k][pos];
        if (sc+maxs[pos+1] >= min) {
          // the score min can be reached
          if (sc > max) {
            // the score will be greater than max for all suffixes
            if(nbocc[length-1][max+1] != nullptr) {
              nbocc[length-1][max+1] = std::make_shared<double>(*nbocc[length-1][max+1] + *nbocc[pos-1][iter->first] * background[k]); //pow(4,length-pos-1) ;
            } else {
              nbocc[length-1][max+1].reset(new double);
              *nbocc[length-1][max+1] = *nbocc[pos-1][iter->first] * background[k]; //pow(4,length-pos-1) ;
            }
            totalOp++;
          } else {              
            if(nbocc[pos][sc] != nullptr) {
              nbocc[pos][sc] = std::make_shared<double>(*nbocc[pos][sc] + *nbocc[pos-1][iter->first] * background[k]);
            } else {
              nbocc[pos][sc].reset(new double);
              *nbocc[pos][sc] = *nbocc[pos-1][iter->first] * background[k];
            }
            totalOp++;
          }
        } 
      }
      iter++;      
    }      
    cerr << "        map size for " << pos << " " << nbocc[pos].size() << endl;
  }
cerr << "!" << endl;
cin.get();
  
  
  delete[] maxs;
  
  return nbocc;
  
  
}




qlonglong Matrix::fastPvalue (Matrix *m, qlonglong alpha) {
  
  
  map<qlonglong, qlonglong> *q = new map<qlonglong, qlonglong> [m->length+1];
  map<qlonglong, qlonglong>::iterator iter;
  
  qlonglong P = 0;
  
  qlonglong *maxm = new qlonglong[m->length+1]; // @ pos i maximum score reachable with the suffix matrix from i to length-1
  
  maxm[m->length] = 0;
  for (int i = m->length-1; i >= 0; i--) {
    maxm[i] = maxm[i+1] + m->maxScoreColumn[i];
  }
  
  // initializes the map at position 0
  for (int k = 0; k < 4; k++) {
    if (m->matInt[k][0]+maxm[1] >= alpha) {
      //cout << "FP: Set " << m->matInt[k][0] <<  " to ";
      q[0][m->matInt[k][0]] += 1;
      //cout << q[0][m->matInt[k][0]] << endl;        
    }
  }
  
  // computes q values for scores strictly greater than alpha
  for (int pos = 1; pos < m->length; pos++) {
    iter = q[pos-1].begin();
    while (iter != q[pos-1].end()) {
      for (int k = 0; k < 4; k++) {
        qlonglong scm = iter->first + m->matInt[k][pos];
        if (scm > alpha) {
          //cout << "Update P from " << P;
          P += iter->second * (qlonglong)pow(4.0,m->length-pos-1); 
          //cout << " to P " << P << endl;
        } else if (scm + maxm[pos+1] > alpha) {
          q[pos][scm] += iter->second; 
        }
      } 
      iter++;      
    }      
    q[pos-1].erase(q[pos-1].begin(),q[pos-1].end());
  }
  
  
  delete[] maxm;
  
  return P;
  
}
