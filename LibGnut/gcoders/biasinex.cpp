
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <memory> 

#include "gcoders/biasinex.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gbias.h"

using namespace std;

namespace gnut {

// constructor   
t_biasinex::t_biasinex( t_gsetbase* s, string version, int sz, string id )
 : t_sinex( s, version, sz, id )
{   
  gtrace(_class+"::construct");
  _allbias = 0;
}


/* -------- *
 * DECODER
 * -------- */
   
int t_biasinex::_decode_comm()
{
  int idx;
  if( (idx = _line.find("PRN"))  != string::npos ) _mapidx["SAT"]  = make_pair(idx,3);
  if( (idx = _line.find("OBS1")) != string::npos ) _mapidx["OBS1"] = make_pair(idx,4);
  if( (idx = _line.find("OBS2")) != string::npos ) _mapidx["OBS2"] = make_pair(idx,4);

  if( (idx = _line.find("BIAS_START____")) != string::npos ) _mapidx["BEG"] = make_pair(idx,14);
  if( (idx = _line.find("BIAS_END______")) != string::npos ) _mapidx["END"] = make_pair(idx,14);
  
  if( (idx = _line.find("__ESTIMATED_VALUE____")) != string::npos ) _mapidx["EST"] = make_pair(idx,21);
  if( (idx = _line.find("_STD_DEV___")) != string::npos ) _mapidx["STD"] = make_pair(idx,11);
  
  return 1;
}

int t_biasinex::_decode_block()
{ 
  gtrace(_class+"::_decode_block");  

  t_sinex::_decode_block();

  // -------- "BIAS/DESCRIPTION" --------
  if( _block.find("BIAS/DESCRIPTION") != string::npos ){

    if(      _line.find(" OBSERVATION_SAMPLING ") != string::npos ){ if(_log) _log->comment(2,"sinex","Read BIAS/SMP: "+cut_crlf(_line.substr(31))); }
    else if( _line.find(" PARAMETER_SPACING ")    != string::npos ){ if(_log) _log->comment(2,"sinex","Read BIAS/SPC: "+cut_crlf(_line.substr(31))); }
    else if( _line.find(" DETERMINATION_METHOD ") != string::npos ){ if(_log) _log->comment(2,"sinex","Read BIAS/MTD: "+cut_crlf(_line.substr(31))); }
    else if( _line.find(" BIAS_MODE ")            != string::npos ){ if(_log) _log->comment(2,"sinex","Read BIAS/MOD: "+cut_crlf(_line.substr(31))); }
    else if( _line.find(" TIME_SYSTEM ")          != string::npos ){ if(_log) _log->comment(2,"sinex","Read BIAS/TSY: "+cut_crlf(_line.substr(31))); }
  }else if(_block.find("BIAS/SOLUTION") != string::npos ){
    
    string prn = "";
    GOBS gobs1, gobs2;
    gobs1 = gobs2 = X;
    t_gtime beg = FIRST_TIME;
    t_gtime end = LAST_TIME;
    double dcb = 0.0;
//  double std = 0.0;
    
    for(auto it = _mapidx.begin(); it != _mapidx.end(); it++){
      size_t pos = it->second.first;
      size_t len = it->second.second;      
      if(it->first == "SAT") prn = _line.substr(pos, len);
      if(it->first == "OBS1") gobs1 = str2gobs(_line.substr(pos, len));
      if(it->first == "OBS2") gobs2 = str2gobs(_line.substr(pos, len));
      if(it->first == "BEG")  beg.from_str("%Y:%j:%s", _line.substr(pos, len));
      if(it->first == "END" && _line.substr(pos, len)!="0000:000:00000") 
							  end.from_str("%Y:%j:%s", _line.substr(pos, len));
      if(it->first == "EST")  dcb = str2dbl(_line.substr(pos, len));
//      if(it->first == "STD") std = str2dbl(_line.substr(pos, len));
    }
    
#ifdef DEBUG
      cout << "DCB decoding: " << _ac << " " << prn << " " << beg.str_ymdhms() << " " << end.str_ymdhms() << " " << gobs2str(gobs1) << " " << gobs2str(gobs2) << " " << dcb << endl;
      int ooo; cin >> ooo;
#endif    

    shared_ptr<t_gbias> p_bias;
    
    if(_allbias) {
      p_bias = make_shared<t_gbias>();
      p_bias->set(beg, end, dcb* 1e-9 * CLIGHT, gobs1, gobs2);
      _allbias->add( _ac, beg, prn, p_bias);
    }

  }
  
  return 1;
}
      
void t_biasinex::_add_data(string id, t_gdata* pt_data)
{
  gtrace(_class+"::_add_data("+id+")");

  if( pt_data == 0 ) return;

  // ALL OBJECTS
  if( pt_data->id_type() == t_gdata::ALLBIAS ){
    if( ! _allbias ) _allbias = dynamic_cast<t_gallbias*>(pt_data);
    else
      if( _log ) _log->comment(0,_class,"Warning: more gallbias instances. Last used only!");
      else            cerr << _class+" - Warning: more gallbias instances. Last used only!\n";
  }
   
  return;
}
   
} // namespace
