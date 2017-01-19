#include "iterativeSEM.hpp"
#include <iostream>
using namespace std;

void split(std::string str, std::string splitBy, std::vector<std::string>& tokens);
// if I understand correctly, bowtie places output in the filename specified
// at the end of the command                              // intermediate.dat, as
//                                                        // defined in iterativeSEM.hpp
void bowtie_genome_map(Dataset &data, int length, const string& genome){
  // investigating bowtie
}
