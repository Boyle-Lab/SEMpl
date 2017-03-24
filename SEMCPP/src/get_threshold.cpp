#include "iterativeSEM.hpp"
#include "../lib/TFM-Pvalue/Matrix.h"
using namespace std;

class Matrix;

const string TEMPFILE = "output/temp.txt";

//REQUIRES: data is a valid Dataset, in that data is not "missing," with PWM_data filled in
//MODIFIES: data
//EFFECTS: returns score threshold
//NOTE: I believe pval == 0.0009765625, also the original version takes a file address I believe

double get_threshold(Dataset & data, double pval){

	pwm_to_tfm(data);
	Matrix m(0.25, 0.25, 0.25, 0.25);

	ofstream temp_out(TEMPFILE);

#ifdef DEBUG
	assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[1].size());
	assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[2].size());
	assert(data.TFM_data.letter_array[0].size() == data.TFM_data.letter_array[3].size());
#endif

	for(int i = 0; i < Dataset::TFMdata::LETTER_NUM; i++){
		for(int j = 0; j < static_cast<int>(data.TFM_data.letter_array[i].size()); j++){
			temp_out << data.TFM_data.letter_array[i][j] << ' ';
		}
		temp_out << '\n';
	}



	temp_out.close();

	m.readJasparMatrix(TEMPFILE);

	// implementation of TFMpvalue.cpp below

	double initial_gran = 0.1;
	m.computesIntegerMatrix(initial_gran);
	long long max = m.maxScore + ceil(m.errorMax + 0.5);
	long long min = m.minScore;
	double pv = 0, ppv = 0;
	long long score = 0;
	double final_score = 0.0;

	for(double granularity = initial_gran; granularity >= 1e-10; granularity /= 10){
		m.computesIntegerMatrix(granularity);
		score = m.lookForScore(min, max, pval, &pv, &ppv);
		min = (score - ceil(m.errorMax + 0.5)) * 10;
		max = (score + ceil(m.errorMax + 0.5)) * 10;
	}
#ifdef DEBUG
	assert(score != 0);
	assert(m.granularity != 0.0);
#endif
	final_score = (score - m.offset) / m.granularity;
#ifdef DEBUG
	assert(final_score != 0.0);
#endif
	return final_score;

}
