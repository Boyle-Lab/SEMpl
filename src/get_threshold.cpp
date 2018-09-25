#include "src/iterativeSEM.hpp"
#include "lib/TFMPvalue/src/Matrix.h"
using namespace std;

class Matrix;


//REQUIRES: data is a valid Dataset, in that data is not "missing," with PWM_data filled in
//MODIFIES: data
//EFFECTS: returns score threshold

double get_threshold(Dataset & data, double pval){
	pwm_to_tfm(data);
	Matrix m;
	int column;
	const string TEMPFILE = data.base_dir + "temp.txt";
	ofstream temp_out(TEMPFILE);

#ifdef DEBUG
    assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[1].size());
    assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[2].size());
    assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[3].size());
    assert(data.TFM_data.letter_array[0].size() > 0);
#endif

#ifdef VERBOSE
	cerr << "INITIAL MATRIX" << endl;
#endif



    cout << "row size " << data.TFM_data.letter_array.size() << " col size " << data.TFM_data.letter_array[0].size() << endl << flush;
    for(column = 0; column < (int)data.TFM_data.letter_array.size(); ++column){
        for(int row = 0; row < static_cast<int>(data.TFM_data.letter_array[0].size()); ++row){

			temp_out << data.TFM_data.letter_array[column][row] << '\t';
// #ifdef VERBOSE
			cerr << data.TFM_data.letter_array[column][row] << '\t';
// #endif
		}
// #ifdef VERBOSE
		cerr << endl;
// #endif
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

	m.toLogOddRatio();

	const double initialGranularity = 0.1;
	const bool forcedGranularity = false;
	const double maxGranularity = 1e-10;
	// const bool sortColumns = false;
	const long long decrgr = 10;
	const double & requestedPvalue = pval;


  // cerr << "### PvalueToScore (pv  " << requestedPvalue << ") #########################################" << endl;



  m.computesIntegerMatrix(initialGranularity);
  long long max = m.maxScore+ceil(m.errorMax+0.5);
  long long min = m.minScore;
  double pv = 0.0;
  long long score = 0;

  for (double granularity = initialGranularity; granularity >= maxGranularity; granularity /= decrgr) {

    // cerr << "Computing rounded matrix with granularity " << granularity << endl;

    m.computesIntegerMatrix(granularity);

    // cerr << "Score range : " << m.scoreRange << endl;
    // cerr << "Min         : " << min << " " << m.minScore << endl;
    // cerr << "Max         : " << max << endl;
    // cerr << "Precision   : " << m.granularity << endl;
    // cerr << "Error max   : " << m.errorMax << endl;
    // cerr << "Computing score for requested pvalue " << requestedPvalue << endl;

    double ppv = 0.0;
    score = m.lookForScore(min, max, requestedPvalue, &pv, &ppv);

    // cerr << "P-Pvalue      : " << ppv << endl;
    // cerr << "Pvalue        : " << pv << endl;
    // cerr << "Rounded score : " << score << endl;
    // cerr << "Real score    : " << ((score-m.offset)/m.granularity) << endl;

    // cerr << "Memory        : " << m.totalMapSize << endl;

    min = (score - ceil(m.errorMax+0.5)) * decrgr;
    max = (score + ceil(m.errorMax+0.5)) * decrgr;


    // cerr << "***********************************************" << endl;

    if (pv == ppv) {

      // cerr << "#####  STOP Pvalue computed  #####" << endl;

      if (!forcedGranularity) {
        break;
      }
    }

  }

  // if (OPTIONS['h']) {
  //   // cout << "Score          : " << ((score-m.offset)/m.granularity) << endl;
  //   // cout << "Pvalue         : " << pv << endl;
  //   // cout << "Granularity    : " << m.granularity << endl;

  //   // cout << "Total map size : " << totalSize << endl;
  //   // cout << "Total op       : " << totalOp << endl;

  // } else {
  //   if (OPTIONS['column']) {
  //     cout << score << " ";
  //   }
  //   // cout << ((score-m.offset)/m.granularity) << " ";
  //   // cout << pv << " ";
  //   // cout << m.granularity << " ";

  //   // cout << totalSize << " " << totalOp << " ";

  //   // cout << endl;
  // }

	double ans = ((score - m.offset ) / m.granularity );

  // free the memory allocated, not typical Rcpp way
  for(int i=0; i<column; i++){
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

	return ans;
}
