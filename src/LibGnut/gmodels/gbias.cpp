
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "gmodels/gbias.h"

using namespace std;

namespace gnut {   
#define DIFF_DCB (86400*31)   // 31 days

// constructor
// ----------
t_gbias::t_gbias()
{
  _beg = FIRST_TIME;
  _end = LAST_TIME;
  
  id_type(t_gdata::BIAS);
  id_group(t_gdata::GRP_MODEL);
}


// destructor
// ----------
t_gbias::~t_gbias(){}

// get single code bias
// ----------
// ----------
double t_gbias::bias(bool meter)
{
  if(meter) return _val;
  else return (_val / CLIGHT) * 1e9;  
}

// add bias
// ----------   
void t_gbias::set(t_gtime beg, t_gtime end, double val, GOBS obs1, GOBS obs2)
{
  #ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
  _gmutex.lock();
  
  _beg = beg;
  _end = end;
  
  _gobs = obs1;
  _ref  = obs2;
  
  _val   = val;
  
  _gmutex.unlock(); return;
}

// add bias
// ----------
void t_gbias::set( double val, GOBS obs1, GOBS obs2)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
  _gmutex.lock();
  
  _gobs = obs1;
  _ref  = obs2;
  
  _val   = val;
  
  _gmutex.unlock(); return;

}

// set reference signal
void t_gbias::ref(GOBS ref){
  _gmutex.lock();

  _ref = ref;

  _gmutex.unlock(); return;
}
  
// test validity   
bool t_gbias::valid(const t_gtime& epo)
{
  _gmutex.lock();
  
  bool ret = true;;

  if(epo < _beg || epo > _end) ret = false;
  else ret = true;
  
  _gmutex.unlock();
  return ret;
}

   
} // namespace
