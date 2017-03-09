CXX=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11 -I.

DEBUG_OR_OPTIMIZE = -DDEBUG -g3 -Wno-variadic-macros -Wno-unused-local-typedefs -Wno-enum-compare

C++SOURCES=$(wildcard src/*.cpp)

C++OBJ=$(subst src/,obj/,$(C++SOURCES))
C++OBJ:=$(C++OBJ:.cpp=.o)

BIGWIGLIBRARY=-lBigWig
BIGWIGLIBRARYLOC=-Llib/libBigWig-master

BOWTIELIBRARY=-lbowtie
BOWTIELIBRARYLOC=-Llib/bowtie-1.0.0

CSOURCES=bw{Read,Stats,Values,Write}.c

TFMHEADER=$(wildcard lib/TFM-Pvalue/*.cpp)
TFMOBJ=$(subst lib/TFM-Pvalue/,obj/,$(TFMHEADER))
TFMOBJ:=$(TFMOBJ:.cpp=.o)

SQLITESOURCE=lib/sqlite3/sqlite3.c
SQLITEOBJ=obj/sqlite3.o


EXECUTABLE=iterativeSEM

all: $(EXECUTABLE)

$(EXECUTABLE): $(C++OBJ) $(SQLITEOBJ) $(TFMOBJ)
		$(CXX) $(C++FLAGS) $(BIGWIGLIBRARYLOC) $(BIGWIGLIBRARY) -o $(EXECUTABLE) $(C++OBJ) $(SQLITEOBJ) $(TFMOBJ)

obj/%.o : lib/TFM-Pvalue/%.cpp
		$(CXX) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c $< -o $@
obj/%.o : src/%.cpp												#pre-requisite   #target
		$(CXX) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) -c $< -o $@
obj/%.o : lib/sqlite3/%.c
		gcc $(DEBUG_OR_OPTIMIZE) -c $< -o $@

clean:
		rm -f obj/*.o
		rm -f $(EXECUTABLE)

.PHONY: all
