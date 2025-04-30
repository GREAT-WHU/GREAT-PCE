
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <sstream>
#include <stdlib.h>
#include <algorithm>

#include "gall/gallpcv.h"
 
using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallpcv::t_gallpcv() 
  : t_gdata(),
    _overwrite(false)
{
  id_type( t_gdata::ALLPCV );
}


// destructor
// ----------
t_gallpcv::~t_gallpcv()
{  
//if( _close_with_warning )
//{	
  for( auto itANT = _mappcv.begin(); itANT != _mappcv.end(); ++itANT ){
    string ant = itANT->first;
    for( auto itNUM = _mappcv[ant].begin(); itNUM != _mappcv[ant].end(); ++itNUM ){
      string num = itNUM->first;
      for( auto itPCV = _mappcv[ant][num].begin(); itPCV != _mappcv[ant][num].end(); ++itPCV ){
	 
	shared_ptr<t_gpcv> gpcv = itPCV->second;
	t_gallnote* gallnote = gpcv->gnote();
	if( gallnote ){
  	  vector<t_gnote> notes = gallnote->mesg();
          std::sort( notes.begin(), notes.end() );
          for( auto it = notes.begin(); it != notes.end(); ++it ){
            if( _log ){ _log->comment( 0, it->func(), it->str()); }
	  }
	}
      }
    }
  }
      
}


// return position (virtual)
// ----------
int t_gallpcv::addpcv( shared_ptr<t_gpcv> pcv )
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  string  key = pcv->pcvkey();
  string  typ = pcv->pcvtyp();
  t_gtime beg = pcv->beg();
  t_gtime end = pcv->end();
   
  pcv->gnote(_gnote); // TRANSFER GNOTE

  if( _mappcv.find(key)           == _mappcv.end() ||
      _mappcv[key].find(typ)      == _mappcv[key].end() ||
      _mappcv[key][typ].find(beg) == _mappcv[key][typ].end() ){

     // new instance
     _mappcv[key][typ][beg] = pcv;
#ifdef DEBUG
     cout << "NEW: " << key << " " << typ << " " << endl;
#endif
  }else{
      
    if( _overwrite ){
#ifdef DEBUG
     cout << "DEL: " << key << " " << typ << " " << endl;
#endif       
     _mappcv[key][typ][beg] = pcv;
#ifdef DEBUG       
     cout << "RENEW: " << key << " " << typ << " " << endl;
#endif
    }else{
      if( _log && _log->verb() >= 1 ) _log->comment(1, "gallpcv", "already exists, not overwritten !\n" );
#ifdef DEBUG       
      cout << "OLD: " << key << " " << typ << " " << endl;
#endif
    }
  }

  if( _log && _log->verb() >= 2 ) 
      _log->comment(2,"gallpcv","add PCV " + key + " " + typ 
	                                   + beg.str_ymdhms(" ") + end.str_ymdhms(" "));
   
  _gmutex.unlock(); return 1;
}

// find PCO+PCV model
// ----------
shared_ptr<t_gpcv> t_gallpcv::gpcv( string ant, string ser, const t_gtime& t )
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  shared_ptr<t_gpcv> gpcv = _find( ant, ser, t );
  _gmutex.unlock(); return gpcv;
}


// return list of available antennas
// ----------
vector<string> t_gallpcv::antennas()
{ 
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<string> all_ant;
     
  map<string,t_map_num>::const_iterator itant = _mappcv.begin();   
      
  for( itant = _mappcv.begin(); itant != _mappcv.end(); ++itant ){
  
    string ant = itant->first;

    map<string,t_map_epo>::iterator itnum;
    for( itnum = _mappcv[ant].begin(); itnum != _mappcv[ant].end(); ++itnum ){

      string num = itnum->first;
     
      all_ant.push_back(ant + " " + num);
    }
  }
  _gmutex.unlock(); return all_ant;
}

// return list of available satellites
// ----------
shared_ptr<t_gpcv> t_gallpcv::find( string ant, string ser, const t_gtime& t )
{
   return _find(ant, ser, t);
}
   
// return list of available satellites
// ----------
shared_ptr<t_gpcv> t_gallpcv::_find( string ant, string ser, const t_gtime& t )
{   
  string ant0 = ant;
  string ser0 = "";

#ifdef DEBUG
  cout << "gallpcv A : ant [" + ant + "]  ser [" + ser + "]  ser0 [" + ser0 + "]\n";
#endif 

  // antenna not found in the list
  if( _mappcv.find(ant) == _mappcv.end() ){
     
    if( ant.size() >= 19 )  ant0 = ant.replace(16,4,"NONE"); 

    if( _mappcv.find(ant0) == _mappcv.end() ){
       ostringstream ostr;
       ostr << "Warning: unknown PCO ( antenna " << left << setw(20) << ant << " not found in ATX ) " << t.str_ymdhms();
       if( _log ) _log->comment(2,"gallpcv",ostr.str() );
       return _pcvnull;
    }
  }

  map<string,t_map_epo>::iterator itser;
  if( ser.find("*") == string::npos ){ // SPECIAL SERIAL DEFINED (OR TYPE CALLIBRATION "")
    for( itser = _mappcv[ant0].begin(); itser != _mappcv[ant0].end(); ++itser ){
      
      // antenna (serial number) found in the list!
      if( _mappcv[ant0].find(ser) != _mappcv[ant0].end() ){     
	ser0 = ser; // REPLACE DEFAULT SER-NUMB (INDIVIDUAL CALLIBRATION)
        break;
      }
    }
  }
   
  // try to find the approximation (ser0=DEFAULT, ser=NOT FOUND)
  if( ser0 == "" && _mappcv[ant].find(ser) == _mappcv[ant].end() ){
     ser = "";  // REPLACE ALWAYS WITH TYPE CALLIBRATION !!!     
  }

  // check/return individual and check time validity
  if( ser0 == ser ){ // SERIAL NUMBER FOUND!
    t_map_epo::iterator it = _mappcv[ant][ser0].begin();
    while( it != _mappcv[ant][ser0].end() ){
      if( ((it->second)->beg() < t || (it->second)->beg() == t)  &&
   	  (t < (it->second)->end() || t == (it->second)->end())   ){
	 return it->second;	      
      }
      it++;
    }
  }

#ifdef DEBUG
  cout << "gallpcv F : ant [" + ant + "]  ser [" + ser + "]  ser0 [" + ser0 + "]\n";   
#endif

   ostringstream ostr;
   ostr << "Warning: unknown PCO ( any antenna alternative [" << ant0 << ","<< ser0 << "] found in ATX ) " << t.str_ymdhms();
   if( _log ) _log->comment(1,"gallpcv",ostr.str() );
   
  return _pcvnull;
}

} // namespace
