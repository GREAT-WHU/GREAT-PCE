/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdlib.h>
#include <iostream>

#include "gdata/gnav.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gnav::t_gnav() : t_geph()
{
  id_type(t_gdata::EPH);
  id_group(t_gdata::GRP_EPHEM);
}


// destructor
// ----------
t_gnav::~t_gnav(){}


// get GNSS NAV validity time [s]
int t_gnav::nav_validity( GSYS gs )
{
  switch( gs ){
    case GPS : return MAX_GPS_TIMEDIFF;
    case GLO : return MAX_GLO_TIMEDIFF;
    case GAL : return MAX_GAL_TIMEDIFF;
    case BDS : return MAX_BDS_TIMEDIFF;
    case QZS : return MAX_QZS_TIMEDIFF;
    case SBS : return MAX_SBS_TIMEDIFF;
    case IRN : return MAX_IRN_TIMEDIFF;
    case GNS : return MAX_NAV_TIMEDIFF;
  }
  return MAX_NAV_TIMEDIFF;
}


// healthy check
// ----------
bool t_gnav::healthy() const
{   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  bool tmp = this->_healthy();
   
  _gmutex.unlock(); return tmp;
}


} // namespace
