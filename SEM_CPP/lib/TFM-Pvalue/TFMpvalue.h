#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
//#include <GetOpt.h>
#include <cstdlib>
#include <stdio.h>

#include "Matrix.h"
#include "ArgumentException.h"


using namespace std;

#define VERBOSE
#define PRINTVERBOSE

#define PV2SC 0
#define SC2PV 1
#define ENUMSC 2
#define DISTRIB 3
#define FASTPVALUE 4
#define LAZY 5

#define TOLOGRATIO 1
#define TOINTEGER 1

#ifdef VERBOSE
#define MEMORYCOUNT
#endif

double _beckstette (Matrix m, map<qlonglong, double> **nbOcc, map<qlonglong, double> **pbuf, int pos, qlonglong score, qlonglong d);

double _beckstettePbuf (Matrix m, map<qlonglong, double> **nbOcc, map<qlonglong, double> **pbuf, int pos, qlonglong score, qlonglong d);


