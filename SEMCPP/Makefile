CXX=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11 -I. -pthread

DEBUG_OR_OPTIMIZE = -O3 
#DEBUG_OR_OPTIMIZE = -g0 -g -DDEBUG 

DISABLE_WARNINGS = -Wno-variadic-macros -Wno-unused-local-typedefs -Wno-enum-compare -Wno-unused-result

C++SOURCES=$(wildcard src/*.cpp)

C++OBJ=$(subst src/,obj/,$(C++SOURCES))
C++OBJ:=$(C++OBJ:.cpp=.o)

BIGWIGLIB=-lBigWig
#BIGWIGLIBRARYLOC=-Llib/libBigWig-master

BOWTIELIB=-lbowtie
#BOWTIELIBRARYLOC=-Llib/bowtie-1.0.0

TFMLIB=-lTFMpvalue
#TFMLIBRARYLOC=-Llib/TFM-Pvalue


LIBLOC=-Llib

CSOURCES=bw{Read,Stats,Values,Write}.c


TFMHEADER=$(wildcard lib/TFM-Pvalue/*.cpp)
TFMOBJ=$(subst lib/TFM-Pvalue/,obj/,$(TFMHEADER))
TFMOBJ:=$(TFMOBJ:.cpp=.o)

SQLITESOURCE=lib/sqlite3/sqlite3.c
SQLITEOBJ=obj/sqlite3.o



EXECUTABLE=iterativeSEM

all: $(EXECUTABLE)

$(EXECUTABLE): $(C++OBJ) $(SQLITEOBJ)
		$(CXX) $(C++FLAGS) $(LIBLOC) -o $(EXECUTABLE) $(C++OBJ) $(SQLITEOBJ) $(TFMLIB) $(BIGWIGLIB) -ldl

obj/%.o : lib/TFM-Pvalue/%.cpp
		$(CXX) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) $(DISABLE_WARNINGS) -c $< -o $@
obj/%.o : src/%.cpp												#pre-requisite   #target
		$(CXX) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) $(DISABLE_WARNINGS) -c $< -o $@
obj/%.o : lib/sqlite3/%.c
		gcc -Wall -Werror -pedantic $(DEBUG_OR_OPTIMIZE) $(DISABLE_WARNINGS) -c $< -o $@

clean:
		rm -f obj/*.o
		rm -f $(EXECUTABLE)

.PHONY: all clean
