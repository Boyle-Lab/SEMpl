CXX=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11 -I.

DEBUG_OR_OPTIMIZE = -DDEBUG -g3 -Wno-variadic-macros

#C++SOURCES=src/iterativeSEM.cpp src/pwm_to_tfm.cpp src/get_threshold.cpp src/accumSummary_scale.cpp src/generatePWMfromSEM.cpp src/generateSNPEffectMatrix.cpp src/Enumerate_kmer.cpp src/common.cpp src/alignToGenomeWrapper.cpp src/checkCache.cpp
C++SOURCES=$(wildcard src/*.cpp)

C++OBJ=$(subst src/,obj/,$(C++SOURCES))
C++OBJ:=$(C++OBJ:.cpp=.o)

BIGWIGLIBRARY=-lBigWig
BIGWIGLIBRARYLOC=-Llib/libBigWig-master
BOWTIELIBRARY=-lbowtie
BOWTIELIBRARYLOC=-Llib/bowtie-1.0.0
CSOURCES=bw{Read,Stats,Values,Write}.c
TFMHEADER=lib/TFM-Pvalue/{Matrix,ParseException,ArgumentException,FileException}.cpp
SQLITESOURCE=lib/sqlite3/sqlite3.c
SQLITEOBJ=obj/sqlite3.o


#C++OBJ=obj/iterativeSEM.o obj/pwm_to_tfm.o obj/get_threshold.o obj/accumSummary_scale.o obj/generatePWMfromSEM.o obj/generateSNPEffectMatrix.o obj/Enumerate_kmer.o obj/common.o obj/alignToGenomeWrapper.o obj/checkCache.o

EXECUTABLE=iterativeSEM.exe

all: $(EXECUTABLE)

iterativeSEM.exe: $(C++OBJ) $(SQLITEOBJ)
		$(CXX) $(C++FLAGS) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) $(TFMHEADER) -o $(EXECUTABLE) $(C++OBJ) $(SQLITEOBJ)

obj/%.o : src/%.cpp												#pre-requisite   #target
		$(CXX) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c $< -o $@
obj/%.o : lib/sqlite3/%.c
		gcc $(DEBUG_OR_OPTIMIZE) -c $< -o $@

clean:
		rm -f obj/*.o
		rm -f $(EXECUTABLE)

.PHONY: all
