/**
*
* @verbatim
    History
    2011-01-10  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file      gnav.h
* @brief     implements ephemerides and navigation classes
*            overwritting mode
*            return first/last
*            limit for number of data inclusions per satellites
*.
* @author    JD
* @version   1.0.0
* @date      2011-01-10
*
*/
#ifndef GNAV_H
#define GNAV_H 
 

#include "gdata/geph.h"
#include "gutils/gsys.h"
#include "gutils/gtime.h"

#define MAX_RINEXN_REC     29   // maximum number of RINEXN records for any system !!

#define MAX_NAV_TIMEDIFF  3600*2    // NAV GNS valitity interval [s]
#define MAX_GPS_TIMEDIFF  3600*4.5  // NAV GPS validity interval [s]
#define MAX_GLO_TIMEDIFF  3600		// NAV GLO validity interval [s] 60*17
#define MAX_GAL_TIMEDIFF  3600*3    // NAV GAL validity interval [s]
#define MAX_BDS_TIMEDIFF  3600      // NAV BDS validity interval [s]
#define MAX_SBS_TIMEDIFF   360      // NAV SBS validity interval [s]
#define MAX_IRN_TIMEDIFF  3600*2    // NAV IRN validity interval [s]
#define MAX_QZS_TIMEDIFF  3600      // NAV QZS validity interval [s]

using namespace std;

namespace gnut {
    /** @brief navdata. */
    typedef double t_gnavdata[MAX_RINEXN_REC];
    /** @brief Class for navigation data storing. */
    class  LibGnut_LIBRARY_EXPORT t_gnav : public t_geph {

     public:
      explicit t_gnav();
      virtual ~t_gnav();
  
      static  int nav_validity( GSYS gs );

      virtual int data2nav( string sat, const t_gtime& ep, const t_gnavdata& data ){ return -1; }
      virtual int nav2data(                                      t_gnavdata& data ){ return -1; }

      virtual int iod() const { return -1; }
      virtual int rec() const { return MAX_RINEXN_REC; }

      virtual bool healthy() const;
      virtual string health_str() const { return _health_str(); }
      virtual int chk(set<string>& msg){ return 1; }

      virtual int freq_num() const { return 255; }

     protected:
      virtual bool _healthy() const { return true; }
      virtual string _health_str() const { return ""; }

     private:

    };

} // namespace

#endif
