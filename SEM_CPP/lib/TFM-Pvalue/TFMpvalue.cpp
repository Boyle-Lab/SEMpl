/*
 *  TFMpvalue.cpp
 *  pvalue
 *
 *  Created by Jean-Stéphane Varré on 02/07/07.
 *  Copyright 2007 LIFL-USTL-INRIA. All rights reserved.
 *
 */

#include "TFMpvalue.h"

map<char,int> OPTIONS;
static string REQUIRED[6] = { "a:t:c:g:m:p:", "a:t:c:g:m:s:" , "a:c:t:g:m:", "a:c:t:g:m:s:S:G:", "a:c:t:g:m:s:G:", "a:c:t:g:m:p:G:" };
static string OPTIONAL[6] = { "whi", "wh", "whG:", "wh", "wh", "wh" };

void stop () {
  string str;
  getline(cin,str);
}

void enumScoreFloatPvalue (Matrix *m, int pos, double score, map<double,int> *t, qlonglong *nbocc, qlonglong pval) {
  
  if (*nbocc < pval) {
    if (pos == m->length) {
      (*t)[score] = 1;
      (*nbocc)++;
    } else {
      for (int k = 0; k < 4; k++) {
        enumScoreFloatPvalue(m,pos+1,score+m->mat[k][pos],t,nbocc,pval);
      }
    }
  }
}

void enumScoreFloat (Matrix *m, int pos, double score, map<double,int> *t) {
  
  if (pos == m->length) {
    (*t)[score] += 1;
  } else {
    for (int k = 0; k < 4; k++) {
      enumScoreFloat(m,pos+1,score+m->mat[k][pos],t);
    }
  }
}

void enumScore (Matrix *m, int pos, qlonglong score, map<qlonglong,int>*t) {
  
  if (pos == m->length) {
    (*t)[score] += 1;
  } else {
    for (int k = 0; k < 4; k++) {
      enumScore(m,pos+1,score+m->matInt[k][pos],t);
    }
  }
}


/**
 * LAZY DISTRIBUTION
 */

double _beckstette (Matrix m, map<qlonglong, double> **nbOcc, map<qlonglong, double> **pbuf, int pos, qlonglong score, qlonglong d);

double _beckstettePbuf (Matrix m, map<qlonglong, double> **nbOcc, map<qlonglong, double> **pbuf, int pos, qlonglong score, qlonglong d) {
  //cout << "d=" << d << " Pbuf_" << pos << " (" << score << ") = ";
  if (pos == -1) { return 0; }
  map<qlonglong,double>::iterator iterPbuf;
  iterPbuf = (*pbuf)[pos].find(score);
  double nb;
  if (iterPbuf == ((*pbuf)[pos]).end()) {
    nb = 0;
  } else {
    nb = iterPbuf->second; // set at the old value
  }
  // compute Pbuf[pos][score]
  for (int k = 0; k < 4; k++) {
    if (m.matInt[k][pos] < m.maxScoreColumn[pos] - d) {
      qlonglong s = score - m.matInt[k][pos];
      //cout << "(" << k << "," << pos << ")" << "->" << matInt[k][pos] << " ";
      if (s <= m.bestScore[pos-1] && s >= 0) {
        nb += _beckstette(m,nbOcc,pbuf,pos-1,s,d) * m.background[k];
      }
    }  
  }
  (*pbuf)[pos][score] = nb;
  return nb;
}

double _beckstette (Matrix m, map<qlonglong, double> **nbOcc, map<qlonglong, double> **pbuf, int pos, qlonglong score, qlonglong d) {
  //cout << "Q_" << pos << " (" << score << ")" << endl;
  if (score < 0 || pos == -1) {
    if (score == 0) return 1;
    else return 0;
  }
  map<qlonglong ,double>::iterator iterNbOcc;
  iterNbOcc = (*nbOcc)[pos].find(score);
  if (iterNbOcc == ((*nbOcc)[pos]).end()) {
    // first compute pbuf
    double nb = _beckstettePbuf(m,nbOcc,pbuf,pos,score,d);
    //      qlonglong nb = (*pbuf)[pos][score];
    //cout << nb << endl;
    // then compute NbOcc
    for (int k = 0; k < 4; k++) {
      if (m.matInt[k][pos] >= m.maxScoreColumn[pos] - d) {
        qlonglong s = score - m.matInt[k][pos];
        if (s <= m.bestScore[pos-1] && s >= 0) {
          nb += _beckstette(m,nbOcc,pbuf,pos-1,s,d) * m.background[k];
        }
      }  
    }
    (*nbOcc)[pos][score] = nb;
  }  
  return (*nbOcc)[pos][score];
}

void testLazyDistrib (Matrix m, double granularity, double requestedPvalue) {

#ifdef MEMORYCOUNT
  qlonglong totalSize = 0;
  qlonglong totalOp = 0;
#endif


#ifdef VERBOSE
  //cerr << "### LAZY DISTRIB (pvalue=" << requestedPvalue << ", with granularity=" << granularity << ") ############################################" << endl;
#endif

  m.computesIntegerMatrix(granularity,true);
  map<qlonglong, double> *nbOcc = new map<qlonglong, double> [m.length+1];
  map<qlonglong, double> *pbuf = new map<qlonglong, double> [m.length+1];
  qlonglong score = m.maxScore+ceil(m.errorMax);
  qlonglong d = 0;
  double pv = 0;
  nbOcc[m.length][score] = pv;
  while (pv <= requestedPvalue) {
    score --;
    //cout << requestedPvalue << "***** BECK for pv " << pv << " and score " << score << " " << (score - m.offset) / m.granularity << endl;
    pv += _beckstette(m,&nbOcc,&pbuf,m.length-1,score,d);
    //      pv += nbOcc[length-1][score];
    nbOcc[m.length][score] = pv;
    d++;
  }
  //cout << requested_pvalue << "***** BECK for pv " << pv << " and score " << score << " " << (score - offset) / granularity << endl;
  //    score++;
  pv = nbOcc[m.length][score];

  /*
  totalMapSize = 0;
  for (int pos = 0; pos <= length; pos++) {
    totalMapSize += nbOcc[pos].size() + pbuf[pos].size();
  }
  */
  
  if (OPTIONS['h']) {
    /*cout << "Score          : " << ((score-m.offset)/m.granularity) << endl;
    cout << "Pvalue         : " << pv << endl;
    cout << "Granularity    : " << m.granularity << endl;
#ifdef MEMORYCOUNT
    cout << "Total map size : " << totalSize << endl;
    cout << "Total op       : " << totalOp << endl;
#endif*/
  } else {  
    /*cout << ((score-m.offset)/m.granularity) << " ";
    cout << pv << " ";
    cout << m.granularity << " ";
#ifdef MEMORYCOUNT
    cout << totalSize << " " << totalOp << " ";
    cout << endl;
#endif*/
  }
}


/**
 * FAST PVALUE
 */

void testFastPvalue (Matrix m, double granularity, double score) {
  
#ifdef VERBOSE
  //cerr << "### FastPvalue (score " << score << ") #########################################" << endl;
#endif

#ifdef MEMORYCOUNT
  qlonglong totalSize = 0;
  qlonglong totalOp = 0;
#endif

  m.computesIntegerMatrix(granularity,true);
  double pvalue = m.fastPvalue(&m,(qlonglong)(score * m.granularity + m.offset));

  if (OPTIONS['h']) {
    /*cout << "Score          : " << score << endl;
    cout << "Pvalue         : " << pvalue << endl;
    cout << "Granularity    : " << m.granularity << endl;
#ifdef MEMORYCOUNT
    cout << "Total map size : " << totalSize << endl;
    cout << "Total op       : " << totalOp << endl;
#endif*/
  } else {  
    /*cout << score << " ";
    cout << pvalue << " ";
    cout << m.granularity << " ";
#ifdef MEMORYCOUNT
    cout << totalSize << " " << totalOp << " ";
    cout << endl;
#endif*/
  }
    
}
  
void testScoreToPvalue (Matrix m, double initialGranularity, double requestedScore, bool forcedGranularity = false, double maxGranularity = 1e-9) {
  
#ifdef VERBOSE
  //cerr << "### ScoreToPvalue (score " << requestedScore << ") #########################################" << endl;
#endif
  
#ifdef MEMORYCOUNT
  qlonglong totalSize = 0;
  qlonglong totalOp = 0;
#endif
  
  qlonglong max;
  qlonglong min;
  double ppv;
  double pv;
  qlonglong score;
  
  for (double granularity = initialGranularity; granularity >= maxGranularity; granularity /= 10) {
#ifdef VERBOSE
    //cerr << "Computing rounded matrix with granularity " << granularity << endl;
#endif
    m.computesIntegerMatrix(granularity);
    max = requestedScore*m.granularity + m.offset + m.errorMax+1;
    min = requestedScore*m.granularity + m.offset - m.errorMax-1;
    score = requestedScore*m.granularity + m.offset;
    
#ifdef VERBOSE
    /*cerr << "Score range : " << m.scoreRange << endl;
    cerr << "Min         : " << min << endl;
    cerr << "Max         : " << max << endl;
    cerr << "Precision   : " << m.granularity << endl;
    cerr << "Error max   : " << m.errorMax << endl;
    cerr << "Computing pvalue for requested score " << requestedScore << " " << score << endl;*/
#endif
    
      
    // computes pvalues for reachable score in range min - max    
#ifdef MEMORYCOUNT
    m.totalMapSize = 0;
    m.totalOp = 0;
#endif
    
    m.lookForPvalue(score,min,max,&ppv,&pv);
    
#ifdef VERBOSE
    /*cerr << "Prev. Pvalue  : " << ppv << endl;
    cerr << "Pvaluex       : " << pv << endl;
    cerr << "Comp. score   : " << score << endl;*/
#endif

#ifdef MEMORYCOUNT
    totalSize += m.totalMapSize;
    totalOp += m.totalOp;
#endif
    
#ifdef VERBOSE
    //cerr << "***********************************************" << endl;
#endif
    
    if (ppv == pv) {
#ifdef VERBOSE
      //cerr << "#####  STOP score computed  #####" << endl;
#endif
      if (!forcedGranularity) {
        break;
      }
    }
    
  }
  
  if (OPTIONS['h']) {
    /*cout << "Score          : " << ((score-m.offset)/m.granularity) << endl;
    cout << "Pvalue         : " << pv << endl;
    cout << "Granularity    : " << m.granularity << endl;
#ifdef MEMORYCOUNT
    cout << "Total map size : " << totalSize << endl;
    cout << "Total op       : " << totalOp << endl;
#endif*/
  } else {  
    /*cout << ((score-m.offset)/m.granularity) << " ";
    cout << pv << " ";
    cout << m.granularity << " ";
#ifdef MEMORYCOUNT
    cout << totalSize << " " << totalOp << " ";
    cout << endl;
#endif*/
    
  }
  
  
}


void testPvalueToScore (Matrix m, double initialGranularity, double requestedPvalue, bool forcedGranularity = false, double maxGranularity = 1e-10, bool sortColumns = false, qlonglong decrgr = 10) {

#ifdef VERBOSE
  //cerr << "### PvalueToScore (pv  " << requestedPvalue << ") #########################################" << endl;
#endif
  
#ifdef MEMORYCOUNT
  qlonglong totalSize;
  qlonglong totalOp;
  totalSize = 0;
  totalOp = 0;
#endif
  
  m.computesIntegerMatrix(initialGranularity);
  qlonglong max = m.maxScore+ceil(m.errorMax+0.5);
  qlonglong min = m.minScore;
  double pv;
  qlonglong score;
  
  for (double granularity = initialGranularity; granularity >= maxGranularity; granularity /= decrgr) {
    
#ifdef VERBOSE    
    //cerr << "Computing rounded matrix with granularity " << granularity << endl;
#endif
    
    m.computesIntegerMatrix(granularity);

#ifdef VERBOSE
    /*cerr << "Score range : " << m.scoreRange << endl;
    cerr << "Min         : " << min << " " << m.minScore << endl;
    cerr << "Max         : " << max << endl;
    cerr << "Precision   : " << m.granularity << endl;
    cerr << "Error max   : " << m.errorMax << endl;
    cerr << "Computing score for requested pvalue " << requestedPvalue << endl;*/
#endif

    double ppv;    
    
#ifdef MEMORYCOUNT
    m.totalMapSize = 0;
    m.totalOp = 0;
#endif
    
    score = m.lookForScore(min,max,requestedPvalue,&pv,&ppv);

#ifdef VERBOSE
    /*cerr << "P-Pvalue      : " << ppv << endl;
    cerr << "Pvalue        : " << pv << endl;
    cerr << "Rounded score : " << score << endl;
    cerr << "Real score    : " << ((score-m.offset)/m.granularity) << endl;*/
#ifdef MEMORYCOUNT
    //cerr << "Memory        : " << m.totalMapSize << " " << totalSize << endl;    
#endif
#endif

#ifdef MEMORYCOUNT
    totalSize += m.totalMapSize;
    totalOp += m.totalOp;
#endif
    
    min = (score - ceil(m.errorMax+0.5)) * decrgr;
    max = (score + ceil(m.errorMax+0.5)) * decrgr;
    
#ifdef VERBOSE
    //cerr << "***********************************************" << endl;
#endif
    if (pv == ppv) {
#ifdef VERBOSE
      //cerr << "#####  STOP Pvalue computed  #####" << endl;
#endif
      if (!forcedGranularity) {        
        break;
      }
    }

  }
  
  if (OPTIONS['h']) {
    /*cout << "Score          : " << ((score-m.offset)/m.granularity) << endl;
    cout << "Pvalue         : " << pv << endl;
    cout << "Granularity    : " << m.granularity << endl;
#ifdef MEMORYCOUNT
    cout << "Total map size : " << totalSize << endl;
    cout << "Total op       : " << totalOp << endl;
#endif*/
  } else {  
    /*if (OPTIONS['i']) {
      cout << score << " ";
    }
    cout << ((score-m.offset)/m.granularity) << " ";
    cout << pv << " ";
    cout << m.granularity << " ";
#ifdef MEMORYCOUNT
    cout << totalSize << " " << totalOp << " ";
#endif
    cout << endl;*/
  }
  
}  //testPvalueToScore


void testDistrib(Matrix m, double granularity, double min, double max) {
  m.computesIntegerMatrix(granularity);
  m.showDistrib(min*m.granularity+m.offset,max*m.granularity+m.offset);  
}


/*void usage (char * const argv[]) {
  cout << "Usage : " << argv[0] << " -a X -t X -g X -c X -m matrix_filename ";
  switch(PROGRAM) {
    case PV2SC:
      cout << " -p pvalue [-w]" << endl;
      cout << "Computes the score threshold associated with a P-value" << endl;
      break;
    case SC2PV:
      cout << "-s threshold [-w]" << endl;
      cout << "Computes the P-value associated with a score threshold" << endl;
      break;
    case FASTPVALUE:
      cout << "-s threshold -G [-w]" << endl;
      cout << "Computes the P-value associated with a score threshold for a given granularity" << endl;
      break;
    case ENUMSC:
      cout << "[-w -G]" << endl;
      cout << "Computes the number of scores of the matrix" << endl;
      break;
    case DISTRIB:
      cout << "-s min_score -S max_score -G granularity" << endl;
      cout << "Computes the distibution of scores between min_score and max_score" << endl;
      break;      
    case LAZY:
      cout << "-p requested pvalue -G granularity" << endl;
      cout << "Computes the score threshold associated with P-value p using the algorithm of Beckstette 2006" << endl;
      break;      
  }
  cout << "  -a -t -c -g : background probabilies" << endl;
  cout << endl;
  cout << "  -m : matrix file" << endl;
#ifndef JASPAR
  cout << "comment on the first line" << endl;
  cout << "A| sc1A sc2A sc3A ..." << endl;
  cout << "C| sc1C sc2C sc3C ..." << endl;
  cout << "G| sc1G sc2G sc3G ..." << endl;
  cout << "T| sc1T sc2T sc3T ..." << endl;
#else
  cout << "sc1A sc2A sc3A ..." << endl;
  cout << "sc1C sc2C sc3C ..." << endl;
  cout << "sc1G sc2G sc3G ..." << endl;
  cout << "sc1T sc2T sc3T ..." << endl;
#endif  
  cout << endl;
  cout << "  -w : the matrix is already a weight matrix, otherwise it is assumed to be a count matrix"<< endl;
  cout << endl;
  cout << "  -p : requested pvalue" << endl;
  cout << endl;
  cout << "  -s : score threshold" << endl;
  cout << endl;
//  cout << "  -G : granularity for integer matrix (a floating number)" << endl;
//  cout << endl;
}

void arguments (int argc, char * const argv[]) {
   // parse options
   char option;
   map<char,bool> opt;
   char options[REQUIRED[PROGRAM].length()+OPTIONAL[PROGRAM].length()+1];
   for (int i = 0; i < REQUIRED[PROGRAM].length(); i++) { options[i] = REQUIRED[PROGRAM][i]; }
   for (int i = 0; i < OPTIONAL[PROGRAM].length(); i++) { options[i+REQUIRED[PROGRAM].length()] = OPTIONAL[PROGRAM][i]; }  
   options[REQUIRED[PROGRAM].length()+OPTIONAL[PROGRAM].length()+1] = '\0';
  while (((option = getopt(argc,argv,options)) != EOF)) {
    if (option == '?') {
      throw new ArgumentException("Bad argument");
    } 
    OPTIONS[option] = optind-1;        
    opt[option] = true;
  }
  for (int i = 0; i < REQUIRED[PROGRAM].length(); i++) {
    if (REQUIRED[PROGRAM][i] != ':' && !opt[REQUIRED[PROGRAM][i]]) throw new ArgumentException("Bad number of args");    
  }
  
}*/



/*int main (int argc, char * const argv[]) {
  
  try {
    arguments(argc,argv);
  } catch (ArgumentException *e) {
    usage(argv);
    exit(1);
  }
  
  
  
  Matrix m(atof(argv[OPTIONS['a']]),atof(argv[OPTIONS['c']]),atof(argv[OPTIONS['g']]),atof(argv[OPTIONS['t']])); 
  
  try {
#ifndef JASPAR
    m.readHorizontalMatrix(argv[OPTIONS['m']]);
#else
    m.readJasparMatrix(argv[OPTIONS['m']]);
#endif   
  } catch (FileException *e) {
    cerr << "Unable to open/read " << argv[OPTIONS['m']] << endl;
    exit(2);
  } catch (ParseException *e) {
    cerr << "The matrix " << argv[1] << " seems to be in a wrong format. The format of the matrix file is : " << endl;
#ifndef JASPAR
    cout << "comment on the first line" << endl;
    cout << "A| sc1A sc2A sc3A ..." << endl;
    cout << "C| sc1C sc2C sc3C ..." << endl;
    cout << "G| sc1G sc2G sc3G ..." << endl;
    cout << "T| sc1T sc2T sc3T ..." << endl;
#else
    cout << "sc1A sc2A sc3A ..." << endl;
    cout << "sc1C sc2C sc3C ..." << endl;
    cout << "sc1G sc2G sc3G ..." << endl;
    cout << "sc1T sc2T sc3T ..." << endl;
#endif  
    cout << endl;
    exit(2);
  }
  
  if (OPTIONS['h']) {
    cout << "Matrix length  : " << m.length << endl;
  } else {
    if (PROGRAM != DISTRIB) {
      cout << m.length << " ";
    }
  }
  
  if (!OPTIONS['w']) {
    m.toLogOddRatio();
  }
  
  
  float start = clock();
  switch (PROGRAM) {
    case PV2SC :
      testPvalueToScore(m,0.1,(atof(argv[OPTIONS['p']])));
      break;
    case SC2PV :
      testScoreToPvalue(m,0.1,atof(argv[OPTIONS['s']]));
      break;
    case ENUMSC :
    {
      qlonglong nbsc = 0;
      if (OPTIONS['G']) {
        m.computesIntegerMatrix(atof(argv[OPTIONS['G']]));
        map<qlonglong,int> t;
        enumScore(&m,0,0,&t);
        nbsc = t.size();
      } else {
        map<double,int> t;
        //        long int sum = 0;      
        enumScoreFloat(&m,0,0,&t);*/
        /*
         map<double,int>::reverse_iterator riter = t.rbegin();
         while (riter->first >= atof(argv[3]) && riter != t.rend()) {
           sum += riter->second;
           riter++;
         }
         */
/*        nbsc = t.size();
      }
      if (OPTIONS['h']) {
        cout << "Number of different scores = " << nbsc << endl;
        cout << "Number of different words  = " << (qlonglong)(pow(4.0,m.length)) << endl;
        
      } else {
        cout << nbsc << " " << (qlonglong)(pow(4.0,m.length)) << endl;
      }
    }
      break;
    case DISTRIB :
      testDistrib(m,atof(argv[OPTIONS['G']]),atof(argv[OPTIONS['s']]),atof(argv[OPTIONS['S']]));
      break;
    case FASTPVALUE:
      testFastPvalue(m,atof(argv[OPTIONS['G']]),atof(argv[OPTIONS['s']]));
      break;
    case LAZY:
      testLazyDistrib(m,atof(argv[OPTIONS['G']]),atof(argv[OPTIONS['p']]));
      break;
  }
  cout << "TIME:" << (clock()-start)/CLOCKS_PER_SEC;

return 0;
}
*/


