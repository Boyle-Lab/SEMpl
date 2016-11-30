/*
 *  ArgumentException.h
 *  pvalue
 *
 *  Created by Jean-Stéphane Varré on 02/07/07.
 *  Copyright 2007 LIFL-USTL-INRIA. All rights reserved.
 *
 */

#ifndef __ARGUMENTEXCEPTION__
#define __ARGUMENTEXCEPTION__

#include <iostream>

using namespace std;

class ArgumentException { 
public:
  ArgumentException () {} 
  ArgumentException(const char *str) { cerr << str << endl;} 
};

#endif