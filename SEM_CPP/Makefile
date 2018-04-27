CXX=g++

C++FLAGS=-Wall -Werror -pedantic --std=c++11 -I. -pthread

#DEBUG_OR_OPTIMIZE = -g -O0 -DDEBUG
DEBUG_OR_OPTIMIZE = -O3 -DDEBUG
#DEBUG_OR_OPTIMIZE = -g2 -DDEBUG

DISABLE_WARNINGS = -Wno-variadic-macros -Wno-unused-local-typedefs -Wno-enum-compare -Wno-unused-result

C++SOURCES=$(wildcard src/*.cpp)

C++OBJ=$(subst src/,obj/,$(C++SOURCES))
C++OBJ:=$(C++OBJ:.cpp=.o)

BIGWIGLIB=-lBigWig
#BIGWIGLIBRARYLOC=-Llib/libBigWig-master

#BOWTIELIB=-lbowtie
#BOWTIELIBRARYLOC=-Llib/bowtie-1.0.0

TFMLIB=-lTFMpvalue
#TFMLIBRARYLOC=-Llib/TFM-Pvalue

LIBLOC=-Llib

STATICLIB=libTFMpvalue.a

CSOURCES=bw{Read,Stats,Values,Write}.c


TFMHEADER=$(wildcard lib/TFM-Pvalue/*.cpp)
TFMOBJ=$(subst lib/TFM-Pvalue/,obj/,$(TFMHEADER))
TFMOBJ:=$(TFMOBJ:.cpp=.o)

SQLITESOURCE=lib/sqlite3/sqlite3.c
SQLITEOBJ=obj/sqlite3.o



EXECUTABLE=iterativeSEM

all: $(EXECUTABLE)

# ./iterativeSEM -PWM examples/MA0114.1.pwm -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak.gz			-big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig -TF_name HNF4A -output results/HNF4A/ -readcache results/HNF4A/HNF4A.cache.db -verbose

test: $(EXECUTABLE)
	export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:$(PWD)/lib"
	./iterativeSEM -PWM examples/MA0114.1.pwm \
			-merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak.gz \
			-big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig \
			-TF_name HNF4A -output results/HNF4A/ \
			-readcache results/HNF4A/HNF4A.cache.db -verbose > 1test.out 2> 1err.out &

ctcf: $(EXECUTABLE)
	export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:$(PWD)/lib"
	./iterativeSEM -PWM /data/data_repo/PWMs/MA0139.1.pwm \
			-merge_file /data/data_repo/ENCODE/H1hESC/stem.narrowPeak.gz \
			-big_wig /data/data_repo/ENCODE/H1hESC/ENCFF000RSD.bigWig \
			-TF_name CTCF -output results/CTCF_H1hESC/ \
			-readcache results/CTCF_H1hESC/CTCF.cache.db -verbose > 1test.out 2> 1err.out


profile: $(EXECUTABLE)
	perf record --call-graph dwarf -e cycles:u nohup ./iterativeSEM -PWM examples/MA0114.1.pwm  -merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak.gz  -big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig   -TF_name HNF4A -output results/HNF4A/ -readcache results/HNF4A/HNF4A.cache.db -verbose > 1test.out 2> 1err.out &


tracetest: $(EXECUTABLE)
	export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:/home/apboyle/SEM_project/SEM_CPP/SEM_CPP/lib"
	strace -f -o /tmp/strace.log -s 200 ./iterativeSEM -PWM examples/MA0114.1.pwm \
			-merge_file examples/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak.gz \
			-big_wig examples/wgEncodeHaibTfbsHepg2Hnf4asc8987V0416101RawRep1.bigWig \
			-TF_name HNF4A -output results/HNF4A/ \
			-readcache results/HNF4A/HNF4A.cache.db -verbose > 1test.out 2> 1err.out &

libexport:
	@cd lib/libBigWig-master && make && mv libBigWig.so ../
	@cd lib/TFM-Pvalue && make SEMCPPobj && mv libTFMpvalue.so ../
	@cd lib/bowtie-1.0.0 && make && mv libbowtie.so ../
	export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:/home/cmorteru/SEM_CPP/SEM_CPP/lib"

$(EXECUTABLE): $(C++OBJ) $(SQLITEOBJ)
		$(CXX) $(C++FLAGS) $(LIBLOC) -o $(EXECUTABLE) $(C++OBJ) $(SQLITEOBJ) $(TFMLIB) $(BIGWIGLIB) -ldl

obj/%.o : lib/TFM-Pvalue/%.cpp
		$(CXX) $(C++FLAGS) -DPROGRAM=0 $(DEBUG_OR_OPTIMIZE) $(DISABLE_WARNINGS) -c $< -o $@
obj/%.o : src/%.cpp												#pre-requisite   #target
		$(CXX) $(C++FLAGS) $(DEBUG_OR_OPTIMIZE) $(DISABLE_WARNINGS) -c $< -o $@
obj/%.o : lib/sqlite3/%.c
		gcc -Wall -Werror -pedantic $(DEBUG_OR_OPTIMIZE) $(DISABLE_WARNINGS) -c $< -o $@



#lib: $(STATICLIB)

#libTFMpvalue.a : $(TFMOBJ)
#	ar rcs libTFMpvalue.a $(TFMOBJ)

clean:
		rm -f obj/*.o
		rm -f $(EXECUTABLE)

.PHONY: all clean lib libexport
