CFLAGS=-O3 -DJASPAR=1
CC=g++
HEADERS=TFMPvalue/src/Matrix.h TFMPvalue/src/ArgumentException.h TFMPvalue/src/FileException.h TFMPvalue/src/ParseException.h
SOURCES=TFMPvalue/src/TFMpvalue.cpp $(HEADERS:.h=.cpp) 
OBJECTS=$(SOURCES:.cpp=.o)

all: SEMCPPobj

SEMCPPobj: $(SOURCES) $(HEADERS)
	${CC} ${CFLAGS} -DPROGRAM=0 -fPIC -shared $(SOURCES) -o libTFMpvalue.so

clean:
	rm -f TFMpvalue-* *.o *~

.PHONY: clean SEMCPPobj        
