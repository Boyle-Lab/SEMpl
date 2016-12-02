CC=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11

C++SOURCES=iterativeSEM.cpp src/pwm_to_tfm.cpp src/get_threshold.cpp src/accumSummary_scale.cpp src/generatePWMfromSEM.cpp generateSNPEffectMatrix.cpp
BIGWIGLIBRARYLOC=-Llib/libBigWig-master
BIGWIGLIBRARY=-lBigWig
CSOURCES=lib/libBigWig-master/bw{Read,Stats,Values,Write}.c
TFMHEADER=lib/TFM-Pvalue/*.cpp

EXECUTABLE=iterativeSEM.exe

all: 
	$(CC) $(C++FLAGS) -O1 $(C++SOURCES) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) $(TFMHEADER) -o $(EXECUTABLE)
debug: 
	$(CC) $(C++FLAGS) -g $(C++SOURCES) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) $(TFMHEADER) -o $(EXECUTABLE)

clean:
	rm $(EXECUTABLE)
	rm lib/libBigWig/*.o
