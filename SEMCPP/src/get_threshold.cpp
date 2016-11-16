#include "../iterativeSEM.hpp"
using namespace std;

void pwm_to_tfm(Dataset & data);

//REQUIRES: data is a valid Dataset, with PWM data filled in
//EFFECTS: return a class object with the data as designed in get_threshold.pl
//NOTE: I believe pval == 0.0009765625, also the original version takes a file address I believe
void get_threshold(Dataset & data, double pval){
	
	pwm_to_tfm(data);

	//UNFINISHED

}
