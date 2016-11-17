#include "../lib/libBigWig-master/BigWig.h"
#include <string>
using namespace std;

void accumySummary_scale(Dataset &data, string hfile, string cfile, int scale){

	// open file or new object, hfile is bigwig
	bigWigFile_t bwFile = bwOpen(hfile.c_str());
}
