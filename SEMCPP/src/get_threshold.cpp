#include "src/iterativeSEM.hpp"
#include "lib/TFM-Pvalue/Matrix.h"
// #include "lib/TFM-Pvalue/TFMpvalue.cpp"
using namespace std;

class Matrix;

const string TEMPFILE = "examples/temp.txt";

//REQUIRES: data is a valid Dataset, in that data is not "missing," with PWM_data filled in
//MODIFIES: data
//EFFECTS: returns score threshold

double get_threshold(Dataset & data, double pval){

	pwm_to_tfm(data);
	Matrix m(0.25, 0.25, 0.25, 0.25);


	ofstream temp_out(TEMPFILE);

#ifdef DEBUG
	assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[1].size());
	assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[2].size());
	assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[3].size());
	assert(data.TFM_data.letter_array[0].size() > 0);
#endif

	for(int i = 0; i < Dataset::TFMdata::LETTER_NUM; ++i){
		for(int j = 0; j < static_cast<int>(data.TFM_data.letter_array[i].size()); ++j){
			temp_out << data.TFM_data.letter_array[i][j] << ' ';
		}
		temp_out << '\n';
	}

	temp_out.close();


//	./TFMpvalue-pv2sc -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m MA0045.pfm -p 1e-5

// my $resultCmd = "bin/TFMpvalue-pv2sc -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m (input) -p $PVAL";
	try{
		m.readJasparMatrix(TEMPFILE);
	}
	catch(...){
		cerr << "problem with readJasparMatrix" << endl;
		exit(1);
	}


	// testPvalueToScore(m, 0.1, { p-value } );

	const double initialGranularity = 0.1;
	const bool forcedGranularity = false;
	const double maxGranularity = 1e-10;
	// const bool sortColumns = false;
	const long decrgr = 10;
	const double requestedPvalue = pval;


#ifdef VERBOSE
  cerr << "### PvalueToScore (pv  " << requestedPvalue << ") #########################################" << endl;
#endif
  
#ifdef MEMORYCOUNT
  long long totalSize;
  long long totalOp;
  totalSize = 0;
  totalOp = 0;
#endif
  
  m.computesIntegerMatrix(initialGranularity);
  long long max = m.maxScore+ceil(m.errorMax+0.5);
  long long min = m.minScore;
  double pv;
  long long score = 0;
  
  for (double granularity = initialGranularity; granularity >= maxGranularity; granularity /= decrgr) {
    
#ifdef VERBOSE    
    cerr << "Computing rounded matrix with granularity " << granularity << endl;
#endif
    
    m.computesIntegerMatrix(granularity);

#ifdef VERBOSE
    cerr << "Score range : " << m.scoreRange << endl;
    cerr << "Min         : " << min << " " << m.minScore << endl;
    cerr << "Max         : " << max << endl;
    cerr << "Precision   : " << m.granularity << endl;
    cerr << "Error max   : " << m.errorMax << endl;
    cerr << "Computing score for requested pvalue " << requestedPvalue << endl;
#endif

    double ppv;    
    
#ifdef MEMORYCOUNT
    m.totalMapSize = 0;
    m.totalOp = 0;
#endif
    
    score = m.lookForScore(min,max,requestedPvalue,&pv,&ppv);

#ifdef VERBOSE
    cerr << "P-Pvalue      : " << ppv << endl;
    cerr << "Pvalue        : " << pv << endl;
    cerr << "Rounded score : " << score << endl;
    cerr << "Real score    : " << ((score-m.offset)/m.granularity) << endl;
#ifdef MEMORYCOUNT
    cerr << "Memory        : " << m.totalMapSize << " " << totalSize << endl;    
#endif
#endif

#ifdef MEMORYCOUNT
    totalSize += m.totalMapSize;
    totalOp += m.totalOp;
#endif
    
    min = (score - ceil(m.errorMax+0.5)) * decrgr;
    max = (score + ceil(m.errorMax+0.5)) * decrgr;
    
#ifdef VERBOSE
    cerr << "***********************************************" << endl;
#endif
    if (pv == ppv) {
#ifdef VERBOSE
      cerr << "#####  STOP Pvalue computed  #####" << endl;
#endif
      if (!forcedGranularity) {        
        break;
      }
    }

  }
  
  if (OPTIONS['h']) {
    cout << "Score          : " << ((score-m.offset)/m.granularity) << endl;
    cout << "Pvalue         : " << pv << endl;
    cout << "Granularity    : " << m.granularity << endl;
#ifdef MEMORYCOUNT
    cout << "Total map size : " << totalSize << endl;
    cout << "Total op       : " << totalOp << endl;
#endif
  } else {  
    if (OPTIONS['i']) {
      cout << score << " ";
    }
    cout << ((score-m.offset)/m.granularity) << " ";
    cout << pv << " ";
    cout << m.granularity << " ";
#ifdef MEMORYCOUNT
    cout << totalSize << " " << totalOp << " ";
#endif
    cout << endl;
  }
  
	return ((score - m.offset ) / m.granularity );

}
