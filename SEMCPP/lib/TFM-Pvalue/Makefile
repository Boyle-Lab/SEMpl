CFLAGS=-O3 -DJASPAR=1
CC=g++
HEADERS=Matrix.h ArgumentException.h FileException.h ParseException.h
SOURCES=TFMpvalue.cpp $(HEADERS:.h=.cpp) 
OBJECTS=$(SOURCES:.cpp=.o)
PROGRAMS=TFMpvalue-pv2sc TFMpvalue-sc2pc TFMpvalue-distrib TFMpvalue-fastpvalue TFMpvalue-lazydistrib

all: TFMpvalue-pv2sc TFMpvalue-sc2pv TFMpvalue-distrib TFMpvalue-fastpvalue TFMpvalue-lazydistrib

SEMCPPobj: $(SOURCES) $(HEADERS)
	${CC} ${CFLAGS} -DPROGRAM=0 $(SOURCES) -o TFMpvalue-pv2sc.o

TFMpvalue-pv2sc: $(SOURCES) $(HEADERS)
	${CC} ${CFLAGS} -DPROGRAM=0 $(SOURCES) -o TFMpvalue-pv2sc

TFMpvalue-sc2pv: $(SOURCES) $(HEADERS)
	${CC} ${CFLAGS} -DPROGRAM=1 $(SOURCES) -o TFMpvalue-sc2pv
	
#enumsc: pvalue.cpp
#	$(CC) $(CFLAGS) -DPROGRAM=2 pvalue.cpp -o enumsc
	
TFMpvalue-distrib: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DPROGRAM=3 $(SOURCES) -o TFMpvalue-distrib

TFMpvalue-fastpvalue: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DPROGRAM=4 $(SOURCES) -o TFMpvalue-fastpvalue

TFMpvalue-lazydistrib: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -DPROGRAM=5 $(SOURCES) -o TFMpvalue-lazydistrib
	
clean:
	rm TFMpvalue-* *.o *~

testpv2sc: TFMpvalue-pv2sc
	./TFMpvalue-pv2sc -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m MA0045.pfm -p 1e-5
	
testsc2pv: TFMpvalue-sc2pv
	./TFMpvalue-sc2pv -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m MA0045.pfm -s 8.77

testdistrib: TFMpvalue-distrib
	./TFMpvalue-distrib -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m MA0045.pfm -s 8 -S 9 -G 1e-5

dist:
	tar -czvf TFMpvalue.tgz $(HEADERS) $(SOURCES) Makefile ReadMe MA0045.pfm

