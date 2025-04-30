
/**
*
* @verbatim
    History
    2012-05-02 /PV: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gstochasticmodel.h
* @brief      Purpose: stochastic models
*.
* @author     PV
* @version    1.0.0
* @date       2012-05-02
*
*/

#ifndef STOCHASTIC_H
#define STOCHASTIC_H 
 

#include "gutils/gtime.h"
#include "newmat/newmat.h"

using namespace std;

namespace gnut {   
	/** @brief class for t_stochastic. */
	class LibGnut_LIBRARY_EXPORT t_stochastic
	{
	 public:
	/** @brief default constructor. */
	  t_stochastic();
	  virtual ~t_stochastic() {};
	  /** @brief get Q. */
	  virtual double getQ()
	  {
		return 0.0;
	  }
  

	 protected:
  
	 private:

	};
	/** @brief class for t_randomwalk derive from t_stochastic. */
	class LibGnut_LIBRARY_EXPORT t_randomwalk : public t_stochastic
	{
	 public:
	 /** @brief default constructor. */
	  t_randomwalk();
	  virtual ~t_randomwalk() {};
   
	  virtual double getQ();
	  void setTprev(const t_gtime&);
	  void setTcurr(const t_gtime&);
	  void updateTime(const t_gtime&);
	  void setq(double q);
	  double get_dt();
  
	 protected:
	 private:
	  t_gtime _Tprev; ///< time of prev
	  t_gtime _Tcurr; ///< time of current
	  double _dSig;   ///< dsigma
	};
	/** @brief class for t_whitenoise derive from t_stochastic. */
	class LibGnut_LIBRARY_EXPORT t_whitenoise : public t_stochastic
	{
	 public:
	  /** @brief constructor 1. */
	  t_whitenoise(double);
	  virtual ~t_whitenoise() {};
   
	   virtual double getQ();
	   void setVar(double);
   
	 private:
	   double _var;
	};
	/** @brief class for t_statemode. */
	class LibGnut_LIBRARY_EXPORT t_statemode
	{
	public:
		/** @brief default constructor. */
		t_statemode();
		t_statemode(int order,double dt,double noise);
		virtual ~t_statemode();

		int order;         ///< order
		Matrix M;          ///< M
		SymmetricMatrix P; ///< P

	private:
		static const double _coeff[6]; ///coff
	
	
	};

} // namespace

#endif // STOCHASTIC_H
