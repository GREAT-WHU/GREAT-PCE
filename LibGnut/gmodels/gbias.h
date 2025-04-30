/**
*
* @verbatim
    History
    2012-11-05  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gbias.h
* @brief      Purpose: implements GNSS code biases
*.
* @author     JD
* @version    1.0.0
* @date       2012-11-05
*
*/

#ifndef GBIAS_H
#define GBIAS_H 


#include <vector>

#include "newmat/newmat.h"
#include "gdata/gdata.h"
#include "gdata/gobsgnss.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {   

    /** @brief class for t_gbias. */
    class LibGnut_LIBRARY_EXPORT t_gbias : public t_gdata {

    public:
       /** @brief default constructor. */
       explicit t_gbias();
       virtual ~t_gbias();
   
       void  set(t_gtime beg, t_gtime end, double d, GOBS obs1, GOBS obs2 = X);    // add single differential bias in meters
       void  set(double d, GOBS obs1, GOBS obs2 = X);                              // add single differential bias in meters
       void  ref(GOBS ref);
   
       double bias(bool meter = true);                                          // get signgle differential bias
   
       GOBS gobs()const{ return _gobs; }
       double val()const { return _val; }
       GOBS ref()  const{ return _ref;   }      
   
       void    beg(t_gtime t){ _beg = t; }          // set/get valid from
       t_gtime beg()const{ return _beg;  }
       void    end(t_gtime t){ _end = t; }          // set/get valid until
       t_gtime end()const{ return _end;  }
   
       bool valid(const t_gtime& epo);
   
       private:

       t_gtime        _beg;        // valid from
       t_gtime        _end;        // valid until
       GOBS           _gobs;       // observation
       GOBS           _ref;        // reference
       double         _val;        // code biases are stored in meters
 

    };

} // namespace

#endif
