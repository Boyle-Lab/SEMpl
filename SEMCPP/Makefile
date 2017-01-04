CC=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11

DEBUG_OR_OPTIMIZE = -g

C++SOURCES=iterativeSEM.cpp src/pwm_to_tfm.cpp src/get_threshold.cpp src/accumSummary_scale.cpp src/generatePWMfromSEM.cpp generateSNPEffectMatrix.cpp src/Enumerate_kmer.cpp src/common.cpp src/alignToGenomeWrapper.cpp src/checkCache.cpp
BIGWIGLIBRARYLOC=-Llib/libBigWig-master
BIGWIGLIBRARY=-lBigWig
CSOURCES=lib/libBigWig-master/bw{Read,Stats,Values,Write}.c
TFMHEADER=lib/TFM-Pvalue/{Matrix,ParseException,ArgumentException,FileException}.cpp

OBJ=iterativeSEM.o pwm_to_tfm.o get_threshold.o accumSummary_scale.o generatePWMfromSEM.o generateSNPEffectMatrix.o Enumerate_kmer.o common.o alignToGenomeWrapper.o checkCache.o

EXECUTABLE=iterativeSEM.exe

all: $(EXECUTABLE)

iterativeSEM.exe: iterativeSEM.o pwm_to_tfm.o get_threshold.o accumSummary_scale.o generatePWMfromSEM.o generateSNPEffectMatrix.o Enumerate_kmer.o common.o alignToGenomeWrapper.o seq_col_to_fa.o checkCache.o
	$(CC) $(C++FLAGS) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) $(TFMHEADER) -o $(EXECUTABLE) $(OBJ)
 

iterativeSEM.o: iterativeSEM.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c iterativeSEM.cpp
pwm_to_tfm.o: src/pwm_to_tfm.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/pwm_to_tfm.cpp
get_threshold.o: src/get_threshold.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/get_threshold.cpp
generatePWMfromSEM.o: src/generatePWMfromSEM.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/generatePWMfromSEM.cpp
generateSNPEffectMatrix.o: generateSNPEffectMatrix.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c generateSNPEffectMatrix.cpp
Enumerate_kmer.o: src/Enumerate_kmer.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/Enumerate_kmer.cpp
common.o: src/common.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/common.cpp
alignToGenomeWrapper.o: src/alignToGenomeWrapper.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/alignToGenomeWrapper.cpp
accumSummary_scale.o: src/accumSummary_scale.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/accumSummary_scale.cpp
seq_col_to_fa.o: src/seq_col_to_fa.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/seq_col_to_fa.cpp
checkCache.o: src/checkCache.cpp
	$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c src/checkCache.cpp

clean:
	rm *.o
	rm $(EXECUTABLE)
