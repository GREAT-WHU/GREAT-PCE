
/**
*
* @verbatim
	History
	2011-04-26  PV: created
	2012-04-06  JD: extracted type conversion utilities from previous utils.h

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gtypeconv.h
* @brief       Purpose: type conversion utilities
* @author      PV
* @version     1.0.0
* @date        2011-04-26
*
*/

#ifndef GTYPECONV_H
#define GTYPECONV_H

#include <cmath>
#include <string>
#include <iomanip>
#include <iostream>
#include "gexport/ExportLibGnut.h"

using namespace std;

namespace gnut
{

	LibGnut_LIBRARY_EXPORT bool double_eq(const double&, const double&);          // double equivalence according to machine capability
	LibGnut_LIBRARY_EXPORT bool float_eq(const float&,   const float&);           // float equivalence according to machine capability
	LibGnut_LIBRARY_EXPORT double dround(double d);                               // round double 

	LibGnut_LIBRARY_EXPORT string dbl2str(const double&, int p = 3);              // double to string conversion (width/digits by iomanip can be added !!!)
	LibGnut_LIBRARY_EXPORT double str2dbl(const string&);                         // string to double conversion (avoiding blanks)
	LibGnut_LIBRARY_EXPORT double strSci2dbl(const string&);                      // string (Scientific) to double conversion (including convert d,D to E!)
#ifdef STR2DBL
	double str2dbl(const char*);            // faster char* to double conversion (avoiding blanks)
#endif

	LibGnut_LIBRARY_EXPORT string int2str(const int&);                            // integer to string conversion (widtht can be added !!!)
	LibGnut_LIBRARY_EXPORT string int2str(const int& num, const int& width);	   // integer to string conversion (have width)
	LibGnut_LIBRARY_EXPORT int    str2int(const string&);                         // string to integer conversion (avoiding blanks)
	LibGnut_LIBRARY_EXPORT string  bl2str(const bool&);

	LibGnut_LIBRARY_EXPORT size_t substitute(string& line, const string& a, const string& b, bool caseSensitive = true);

	LibGnut_LIBRARY_EXPORT string ltrim(const string&);            // return trimmed string leading spaces
	LibGnut_LIBRARY_EXPORT string rtrim(const string&);            // return trimmed string trailing spaces
	LibGnut_LIBRARY_EXPORT string  trim(const string&);            // returm trimmed string leading and trailing spaces

	LibGnut_LIBRARY_EXPORT double frac(double x);
	LibGnut_LIBRARY_EXPORT string floatWoutZero(double x, int i = 3);

} // namespace

#endif