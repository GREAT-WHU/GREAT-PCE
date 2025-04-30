
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include "gprod/gprod.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gprod::t_gprod(const t_gtime& t, shared_ptr<t_gobj> pt)
 : t_gdata(),
    _epo(t),
    _obj(pt)
{
   gtrace("t_gprod::t_gprod");
   id_group(GRP_PRODUCT);
   id_type(NONE);
}

// destructor
// ----------
t_gprod::~t_gprod()
{
}


// set value + rms
// ----------
int t_gprod::set_val( const string& str, const double& val, const double& rms )
{
  _prod[str] = make_pair(val,rms);
  return 0;
}


// get value
// ----------
int t_gprod::get_val( const string& str, double& val, double& rms )
{
  if( _prod.find(str) != _prod.end() ){
    val = _prod[str].first;
    rms = _prod[str].second;
    return 0;
  }

  return -1;
}

// get value
// ----------
int t_gprod::get_val( const string& str, double& val)
{
  double rms;
  int irc = this->get_val(str, val, rms);

  return irc;
}

// list id
// ----------
set<string> t_gprod::list_id()
{
  set<string> list_id;
  for( itPROD = _prod.begin(); itPROD != _prod.end(); ++itPROD ){
    list_id.insert( itPROD->first );
  }

  return list_id;
}

} // namespace
