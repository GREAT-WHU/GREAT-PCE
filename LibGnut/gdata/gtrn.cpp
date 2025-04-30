/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include "gdata/gtrn.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gtrn::t_gtrn()
  : t_gobj(),
    _channel(255)
    // t_gdata::id_type(TRN) // NEFUNGUJE? 
{
//  cout << "CONSTRUCTOR t_gtrn \n"; cout.flush();
  id_type(TRN);
}

// destructor
// ----------
t_gtrn::~t_gtrn()
{ 
//  boost::mutex::scoped_lock lock(_mutex);

//  _mapchk.clear();
//  cout << "DESTRUCTOR t_gtrn \n"; cout.flush();
}


// add rinex header
// ----------------
void t_gtrn::header(const t_rxnhdr& hdr, string path)
{
  _gmutex.lock();

  t_header_pair pair = make_pair(path,hdr);
  _headers.push_back( pair );
  
  _gmutex.unlock();
}


// get all rinex header
// -------------------
t_gtrn::t_header t_gtrn::headers() const
{   
  return _headers;
}


// get one rinex header
// -------------------
t_rxnhdr t_gtrn::header(string path) const
{
  _gmutex.lock();

  t_rxnhdr rxnhdr;  
  for( auto it = _headers.begin(); it != _headers.end(); ++it ){
    if( it->first.compare(path) == 0 ){ return it->second; }
  }

  return rxnhdr;
}


// set channel
// ----------
void t_gtrn::channel(int chk)
{
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
   _gmutex.lock();
   _channel = chk;
   _gmutex.unlock();

   return;
}


// get channel
// ----------
int t_gtrn::channel() const
{
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
   _gmutex.lock();
   int tmp = _channel;
   _gmutex.unlock();
   
   return tmp;
}


} // namespace
