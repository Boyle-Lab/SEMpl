#include "iterativeSEM.cpp"
using namespace std;

void changeBase(Dataset &data, int position, string nucleotide, vector<string> &new_kmer_vec){
  // find new_kmer
  // performs operations on Enumerated_kmer, or kmerHash
  // remove \r from line, not needed here, as stored in memory
  for(auto pair : data.kmerHash){
    if(position > static_cast<int>(pair.first.size())){
      cout << "Position is greater than size!!\n\tEXITING";
      exit(1);
    }
    string firsthalf = pair.first.substr(0, position);
    string secondhalf = pair.first.substr(position + 1);
    if(data.settings.debug){
      cout << "DEBUG: " << pair.second << " has a firsthalf of " << firsthalf
        << " and a second half of " << secondhalf << endl;
    }
    string new_kmer = firsthalf + nucleotide + secondhalf;
    new_kmer_vec.push_back(new_kmer);
  }
}
