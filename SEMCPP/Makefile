CC=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11

DEBUG_OR_OPTIMIZE = -g

C++SOURCES=src/iterativeSEM.cpp src/pwm_to_tfm.cpp src/get_threshold.cpp src/accumSummary_scale.cpp src/generatePWMfromSEM.cpp src/generateSNPEffectMatrix.cpp src/Enumerate_kmer.cpp src/common.cpp src/alignToGenomeWrapper.cpp src/checkCache.cpp
C++OBJ=$(subst src/,obj/,$(C++SOURCES))
C++OBJ:=$(C++OBJ:.cpp=.o)
BIGWIGLIBRARY=-lBigWig
BIGWIGLIBRARYLOC=-Llib/libBigWig-master
CSOURCES=bw{Read,Stats,Values,Write}.c
TFMHEADER=lib/TFM-Pvalue/{Matrix,ParseException,ArgumentException,FileException}.cpp
SQLITESOURCE=lib/sqlite3/sqlite3.c

#C++OBJ=obj/iterativeSEM.o obj/pwm_to_tfm.o obj/get_threshold.o obj/accumSummary_scale.o obj/generatePWMfromSEM.o obj/generateSNPEffectMatrix.o obj/Enumerate_kmer.o obj/common.o obj/alignToGenomeWrapper.o obj/checkCache.o

EXECUTABLE=iterativeSEM.exe

all: $(EXECUTABLE)

iterativeSEM.exe: $(C++OBJ)
		$(CC) $(C++FLAGS) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) $(TFMHEADER) -o $(EXECUTABLE) $(C++OBJ)

obj/%.o : src/%.cpp												#pre-requisite   #target
		$(CC) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c $< -o $@
obj/%.o : lib/sqlite3/%.c
		gcc $(DEBUG_OR_OPTIMIZE) -c $< -o $@

clean:
		rm *.o
		rm src/*.o
		rm $(EXECUTABLE)

.PHONY: all
