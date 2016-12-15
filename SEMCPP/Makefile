CC=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11

DEBUG_OR_OPTIMIZE = -g

C++SOURCES=iterativeSEM.cpp src/pwm_to_tfm.cpp src/get_threshold.cpp src/accumSummary_scale.cpp src/generatePWMfromSEM.cpp generateSNPEffectMatrix.cpp src/Enumerate_kmer.cpp src/common.cpp src/alignToGenomeWrapper.cpp
BIGWIGLIBRARYLOC=-Llib/libBigWig-master
BIGWIGLIBRARY=-lBigWig
CSOURCES=lib/libBigWig-master/bw{Read,Stats,Values,Write}.c
TFMHEADER=lib/TFM-Pvalue/{Matrix,ParseException,ArgumentException,FileException}.cpp

EXECUTABLE=iterativeSEM.exe

all:
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) $(C++SOURCES) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) $(TFMHEADER) -o $(EXECUTABLE)
debug:
	$(CC) $(C++FLAGS) -g $(C++SOURCES) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) $(TFMHEADER) -o $(EXECUTABLE)

clean:
	rm $(EXECUTABLE)
	rm lib/libBigWig/*.o
