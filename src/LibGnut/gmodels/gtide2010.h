
/**
*
* @verbatim
    History
    2012-11-05  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* References:
    [1]  Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
         IERS Technical Note No. 36, BKG (2010) *
*.
* @file       gtide2010.h
* @brief      Purpose: implements tides
*.
* @author     JD
* @version    1.0.0
* @date       2012-11-05
*
*/

#ifndef GTIDE2010_H
#define GTIDE2010_H
 

#include <vector>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gmodels/gephplan.h"
#include "gmodels/gtide.h"

using namespace std;

namespace gnut {
   
    /** @brief class for t_gtide2010 derive from t_gtide. */
    class LibGnut_LIBRARY_EXPORT t_gtide2010 : public t_gtide {

     public:
      /** @brief constructor 1. */
      t_gtide2010(t_glog* l);
      virtual ~t_gtide2010();

      virtual t_gtriple tide_searth(const t_gtime& epo, t_gtriple& xyz);        // solid earth tides
      virtual t_gtriple tide_pole();                                            // pole tides
      virtual t_gtriple load_ocean(const t_gtime& epoch, const string& site, const t_gtriple& xRec);   // ocean tide loading
      virtual t_gtriple load_atmosph();                                         // atmospheric tide loading

     protected:
       // methods for ocean tides
       int _interpolate(const t_gtime& epo, const Matrix& otl_in, const Matrix& doods, Matrix& otl_out);
       int _tidefreq(const t_gtime& epo, const ColumnVector& idd, double& freq, double& phase);
       t_gtime _epo_save;
       ColumnVector _D;
       ColumnVector _DD;
    };

} // namespace

#endif
