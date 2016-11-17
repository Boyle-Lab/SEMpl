CC=g++

CFLAGS=-Wall -Werror -pedantic -O1 --std=c++11

SOURCES=iterativeSEM.cpp src/pwm_to_tfm.cpp src/get_threshold.cpp src/accumSummary_scale.cpp

EXECUTABLE=iterativeSEM.exe

all:
	$(CC) $(CFLAGS) $(SOURCES) -o $(EXECUTABLE)

clean:
	rm $(EXECUTABLE)
