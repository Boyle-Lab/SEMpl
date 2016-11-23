CC=g++

C++FLAGS=-Wall -Werror -pedantic -O1 --std=c++11

C++SOURCES=iterativeSEM.cpp src/pwm_to_tfm.cpp src/get_threshold.cpp src/accumSummary_scale.cpp 
BIGWIGLIBRARYLOC=-Llib/libBigWig-master
BIGWIGLIBRARY=-lBigWig
CSOURCES=lib/libBigWig-master/bw{Read,Stats,Values,Write}.c

EXECUTABLE=iterativeSEM.exe

all:
	$(CC) $(C++FLAGS) $(C++SOURCES) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) -o $(EXECUTABLE)

clean:
	rm $(EXECUTABLE)
	rm lib/libBigWig/*.o
