#include <vector>
#include <string>
#include "iterativeSEM.hpp"
using namespace std;

//REQUIRES: pwm data in data is filled in with correct data
//MODIFIES: TFM_data within data
//EFFECTS: puts a, c, g, t data into tfm_format
void pwm_to_tfm(Dataset & data){
	data.TFM_data.letter_array[0].clear();
	data.TFM_data.letter_array[1].clear();
	data.TFM_data.letter_array[2].clear();
	data.TFM_data.letter_array[3].clear();
	for(int i = 0; i < (int)data.PWM_data.matrix_arr.size(); ++i){
		data.TFM_data.letter_array[0].push_back(data.PWM_data.matrix_arr[i][0]);
		data.TFM_data.letter_array[1].push_back(data.PWM_data.matrix_arr[i][1]);
		data.TFM_data.letter_array[2].push_back(data.PWM_data.matrix_arr[i][2]);
		data.TFM_data.letter_array[3].push_back(data.PWM_data.matrix_arr[i][3]);
	}
}
