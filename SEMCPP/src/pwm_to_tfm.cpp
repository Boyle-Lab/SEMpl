#include <vector>
#include <string>
#include "../IterativeSEM.hpp"
using namespace std;

//REQUIRES: pwm data in data is filled in with correct data
//MODIFIES: tfm_format within data
//EFFECTS: puts a, c, g, t into tfm_format
void pwm_to_tfm(Dataset & data){
	
	vector<char> a, c, g, t;

	for(int i = 0; i < Dataset::PWM::MATRIX_SIZE; i++){
		a.push_back(data.PWM_data.matrix_arr[i][1]);
		c.push_back(data.PWM_data.matrix_arr[i][2]);
		g.push_back(data.PWM_data.matrix_arr[i][3]);
		t.push_back(data.PWM_data.matrix_arr[i][4]);
	}
}
