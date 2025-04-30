
/**
*
* @verbatim
    History
    2011-04-20  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gpoly.h
* @brief      Purpose: implements polynomial approximation
*             directly interpolate
*             get/fit/evaluate polynomials
*.
* @author     JD
* @version    1.0.0
* @date       2011-04-20
*
*/

#ifndef GPOLY_H
#define GPOLY_H


#include "newmat/newmat.h"
#include "newmat/newmatio.h"

#include <vector>
#include <iostream>

using namespace std;

namespace gnut {   
	/** @brief class for t_gpoly. */
	class t_gpoly {
	 public:
	  /** @brief default constructor. */
	   t_gpoly();
	  ~t_gpoly();

	   void reset();

	   int            degree()       const { return _ncoeff-1;}  // get polynomial degree
	   int            ncoeff()       const { return _ncoeff;  }  // get # of coefficients
	   double         rms()          const { return _rms;     }  // get accuracy
	   bool           valid()        const { return _valid;   }  // get validity
	   double         span()         const { return _span;    }  // get x-span
	   double         xref()         const { return _xref;    }  // get x-reference value
   
	   vector<double> polynomials()  const { return _coef;    }  // get polynomials

	   int  interpolate(const vector<double>& X,const vector<double>& Y, const double& x, double& y, double& dy);
	   int  polynomials(const vector<double>& X,const vector<double>& Y);   
	   void evaluate(double x, int I, double& y);

	   int fitpolynom(const vector<double>& X,     // x data  time-difference
			  const vector<double>& Y,     // y data
			  int   N,  double tunit,      // degree of polynom and X-time units [sec] // 
			  const t_gtime& t);           // reference time  //

	 private:    
	   bool            _valid;     // are coefficients valid ?
	   int             _ncoeff;    // polynomial order (n) for n+1 points
	   double          _rms;       // RMS of fitted coefficients  [meters]
	   double          _xref;      // x reference value is always 0.0 !!!!!!
	   double          _span;      // x span
	   vector<double>  _coef;      // coefficients
	};


	double LibGnut_LIBRARY_EXPORT lagrange_interpolate(const vector<double>& X, const vector<double>& Y, double x, const bool& lvel = false);


} // namespace

#endif
