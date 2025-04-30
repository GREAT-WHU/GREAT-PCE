/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdlib.h>
#include <iostream>

#include "gdata/geph.h"
#include "gutils/gtypeconv.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_geph::t_geph()
 : t_gdata(),
   _sat(""),
   _epoch(t_gtime::GPS),
   _validity(true),
   _gio_ptr(0)
{
  _type  = EPHGPS;
}


// destructor
// ----------
t_geph::~t_geph(){
}


// gsys
// ----------
GSYS t_geph::gsys() const
{   
   gtrace("t_geph::gsys");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  GSYS tmp = GNS;
  _gmutex.lock();
  if( _valid() ) tmp = t_gsys::char2gsys( _sat[0] );
  _gmutex.unlock();
  return tmp;
}


// gsat
// ----------
string t_geph::gsat() const
{   
   gtrace("t_geph::gsat");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  string tmp = "";
  _gmutex.lock();
  if( _valid() ) tmp = _sat;
  _gmutex.unlock();
  return tmp;
}


// get parameter value
// ----------
t_timdbl t_geph::param( const NAVDATA& n )
{
  t_timdbl tmp;
  return tmp;
}


// get parameter value
// ----------
int t_geph::param( const NAVDATA& n, double val )
{
  return 0;
}


// get parameter value
// ----------
bool t_geph::param_cyclic( const NAVDATA& n )
{
  if( n == 2 || n == 6 ) return true;
  return false;
}

   
// valid
// ----------
bool t_geph::valid()
{
   gtrace("t_geph::valid");
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  bool tmp = this->_valid();
  _gmutex.unlock();  
  return tmp;
}


// clean data
// ----------
void t_geph::clear()
{ 
   gtrace("t_geph::clear");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  this->_clear();
  _gmutex.unlock();
   return;
}


// clean internal function
// ----------
void t_geph::_clear()
{
   gtrace("t_geph::_clear");   
   
  _sat.clear();
  _epoch = FIRST_TIME;
}


// clean internal function
// ----------
bool t_geph::_valid() const
{   
   gtrace("t_geph::_valid");   
   
  if( !_validity   ||
      _sat.empty() || 
      _sat   == "" ||
      _epoch == FIRST_TIME ) return false;

  return true; 
}

} // namespace
