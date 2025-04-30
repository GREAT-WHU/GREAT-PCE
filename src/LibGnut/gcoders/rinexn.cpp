/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <memory>

#include "gdata/gnav.h"
#include "gdata/geph.h"
#include "gdata/gtrn.h"
#include "gdata/grxnhdr.h"
#include "gall/gallobj.h" 
#include "gall/gallnav.h"
#include "gutils/gtime.h"
#include "gcoders/rinexn.h"
#include "gutils/gtypeconv.h"
#include "gutils/gfileconv.h"
#include "gutils/gtetrad.h"
#include "gutils/gsys.h"

using namespace std;

namespace gnut {

/* ----------
 * constructor
 */
t_rinexn::t_rinexn( t_gsetbase* s, string version, int sz )
: t_gcoder( s, version, sz )
{
  _gnsssys     = 'G';
  _rxnhdr_all  = true;
//_check    = true;
  _check_dt = FIRST_TIME;
  
  _gset(s); // HAVE TO BE EXPLICITLY CALLED HERE (AT THE END OF CONSTRUCTOR)
}


/* ----------
 * NAV-RINEX header
 */
int t_rinexn::decode_head( char* buff, int sz, vector<string>& errmsg)
{ 

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();

  if( t_gcoder::_add2buffer( buff, sz) == 0 ){ _mutex.unlock(); return 0; };  

  string line;
  int consume = 0;
  int tmpsize = 0;
  while( ( tmpsize = t_gcoder::_getline( line )) >= 0 ){
    
    consume += tmpsize;
    
    // -------- "RINEX VERSION" --------
    if( line.find("RINEX VERSION",60) != string::npos ){          // first line
      switch ( line[20] ){
        case 'N': _gnsssys = 'G'; break; // Navigation data - according to Rinex specification
        case 'G': _gnsssys = 'R'; break; // GLONASS NAVIGATION - occures sometimes in brdc
        default : { string lg("warning: not rinex navigation data file");
          if( _log )  _log->comment(0,"rinexn", lg);
          else  cerr << lg;
          _mutex.unlock(); return -1;
        }
      }
       
      switch ( line[40] ){
        case 'G': _gnsssys = 'G'; break; // GPS
        case 'R': _gnsssys = 'R'; break; // GLONASS
        case 'E': _gnsssys = 'E'; break; // GALILEO
        case 'J': _gnsssys = 'J'; break; // QZSS
        case 'S': _gnsssys = 'S'; break; // SBAS
        case 'I': _gnsssys = 'I'; break; // IRNSS
        case 'M': _gnsssys = 'M'; break; // MIXED
        case ' ': { 
          if( line[20] == 'N' ) _gnsssys = 'G';
          if( line[20] == 'G' ) _gnsssys = 'R';
          if( _log ) _log->comment(0,"rinexn","warning - RINEXN system not defined, used "+t_gsys::char2str(_gnsssys));
          break; 
        }
        default : { 
          string lg("warning: not supported satellite system "+line.substr(40,1));
          if( _log ) _log->comment(0,"rinexn", lg); 
          else       cerr << lg << endl;
          // _mutex.unlock(); return -1;
        }
      }
       

      _version = trim(line.substr(0,9));

      _rxnhdr.path( _fname );

      if( _log && substitute(_version, " ", "") > 0 ){
        _log->comment(2,"rinexn", "reading VER: " + _version + " SYS: " + string(1,_gnsssys) );
      }

    // -------- "PGM / RUN BY / DATE" --------
    }else if( line.find("PGM / RUN BY / DATE",60) != string::npos ){
      _rxnhdr.program( trim(line.substr( 0,20)) );
      _rxnhdr.runby  ( trim(line.substr(20,20)) );
      t_gtime gtime(t_gtime::UTC);
      if( line.substr(56,3) != "UTC" ) gtime.tsys(t_gtime::LOC);
      
      if(      gtime.from_str("%Y%m%d %H%M%S",    line.substr(40,15)) == 0 ) ;
      else if( gtime.from_str("%Y-%m-%d %H-%M-%S",line.substr(40,20)) == 0 ) ;
      else{    gtime = FIRST_TIME; }
       _rxnhdr.gtime(gtime);
      
      if( _log ) _log->comment(2,"rinexn","PGM / RUN BY / DATE: " + _rxnhdr.program()
                                                            + " " + _rxnhdr.runby()
                                                            + " " + _rxnhdr.gtime().str_ymdhms() );

    // -------- "IONOSPHERIC CORR" --------
    }else if( line.find("IONOSPHERIC CORR",60) != string::npos ){
     
      IONO_CORR   IO = str2iono_corr(line.substr(0,4)); // cout << "IO " << line.substr(0,4);
      t_iono_corr io;

      io.x0 = strSci2dbl(line.substr(5+ 0, 12)); // cout << " [" + line.substr(5+ 0, 12) + "]";
      io.x1 = strSci2dbl(line.substr(5+12, 12)); // cout << " [" + line.substr(5+12, 12) + "]";
      io.x2 = strSci2dbl(line.substr(5+24, 12)); // cout << " [" + line.substr(5+24, 12) + "]";
      io.x3 = strSci2dbl(line.substr(5+36, 12)); // cout << " [" + line.substr(5+36, 12) + "]\n";
                                     
      _rxnhdr.iono_corr( IO, io );

	  map<string, t_gdata*>::iterator it = _data.begin();
	  while (it != _data.end()) {
		  if (it->second->id_type() == t_gdata::ALLNAV ||
			  it->second->id_group() == t_gdata::GRP_EPHEM) {

			  ((t_gallnav*)it->second)->add_iono_corr(IO, io);
		  }
		  it++;
	  }

      if( _log ) _log->comment(2,"rinexn","IONOSPHERIC CORR " + iono_corr2str(IO)
                                                              + dbl2str(io.x0)
                                                              + dbl2str(io.x1)
                                                              + dbl2str(io.x2)
                                                              + dbl2str(io.x3)
                                                        + " " + base_name(_fname)
                              );


    // -------- "ION ALPHA" --------
    }else if( line.find("ION ALPHA",60) != string::npos ){
     
      IONO_CORR   IO = IO_GPSA;
      t_iono_corr io;

      io.x0 = strSci2dbl(line.substr(2+ 0, 12)); // cout << " [" + line.substr(5+ 0, 12) + "]";
      io.x1 = strSci2dbl(line.substr(2+12, 12)); // cout << " [" + line.substr(5+12, 12) + "]";
      io.x2 = strSci2dbl(line.substr(2+24, 12)); // cout << " [" + line.substr(5+24, 12) + "]";
      io.x3 = strSci2dbl(line.substr(2+36, 12)); // cout << " [" + line.substr(5+36, 12) + "]\n";
                                     
      _rxnhdr.iono_corr( IO, io );
	  map<string, t_gdata*>::iterator it = _data.begin();
	  while (it != _data.end()) {
		  if (it->second->id_type() == t_gdata::ALLNAV ||
			  it->second->id_group() == t_gdata::GRP_EPHEM) {

			  ((t_gallnav*)it->second)->add_iono_corr(IO, io);
		  }
		  it++;
	  }

      if( _log ) _log->comment(2,"rinexn","ION ALPHA " + iono_corr2str(IO)
                                                       + dbl2str(io.x0)
                                                       + dbl2str(io.x1)
                                                       + dbl2str(io.x2)
                                                       + dbl2str(io.x3)
                                                 + " " + base_name(_fname)
                              );

    // -------- "ION BETA" --------
    }else if( line.find("ION BETA",60) != string::npos ){
     
      IONO_CORR   IO = IO_GPSB;
      t_iono_corr io;

      io.x0 = strSci2dbl(line.substr(2+ 0, 12)); // cout << " [" + line.substr(5+ 0, 12) + "]";
      io.x1 = strSci2dbl(line.substr(2+12, 12)); // cout << " [" + line.substr(5+12, 12) + "]";
      io.x2 = strSci2dbl(line.substr(2+24, 12)); // cout << " [" + line.substr(5+24, 12) + "]";
      io.x3 = strSci2dbl(line.substr(2+36, 12)); // cout << " [" + line.substr(5+36, 12) + "]\n";
                                     
      _rxnhdr.iono_corr( IO, io );
	  map<string, t_gdata*>::iterator it = _data.begin();
	  while (it != _data.end()) {
		  if (it->second->id_type() == t_gdata::ALLNAV ||
			  it->second->id_group() == t_gdata::GRP_EPHEM) {

			  ((t_gallnav*)it->second)->add_iono_corr(IO, io);
		  }
		  it++;
	  }

      if( _log ) _log->comment(2,"rinexn","ION BETA "  + iono_corr2str(IO)
                                                       + dbl2str(io.x0)
                                                       + dbl2str(io.x1)
                                                       + dbl2str(io.x2)
                                                       + dbl2str(io.x3)
                                                 + " " + base_name(_fname)
                              );

    // -------- "TIME SYSTEM CORR" --------
    }else if( line.find("TIME SYSTEM CORR",60) != string::npos ){

      TSYS_CORR   TS = str2tsys_corr(line.substr(0,4)); // cout << "TS " << line.substr(0,4);
      t_tsys_corr ts;
     
      if( line.substr(5,2) != "  " ){ // ELIMINATE INCORRECT HEADERS FROM SOME RECEIVERS!

        ts.a0 = strSci2dbl(line.substr(5,  17)); // cout << " [" + line.substr( 5, 17) + "]";
        ts.a1 = strSci2dbl(line.substr(22, 16)); // cout << " [" + line.substr(22, 16) + "]";
        ts.T  =    str2int(line.substr(38,  7)); // cout << " [" + line.substr(38,  7) + "]";
        ts.W  =    str2int(line.substr(45,  5)); // cout << " [" + line.substr(45,  5) + "]\n";
                                     
        _rxnhdr.tsys_corr( TS, ts );

        if( _log ) _log->comment(2,"rinexn","TIME SYSTEM CORR " + tsys_corr2str(TS)
                                                                + dbl2str(ts.a0*1e9,6)
                                                                + dbl2str(ts.a1*1e9,6)
                                                          + " " + int2str(ts.T)
                                                          + " " + int2str(ts.W)
                                                          + " " + base_name(_fname)
                                 );
        }

    // -------- "DELTA-UTC" --------
    }else if( line.find("DELTA-UTC: A0,A1,T,W",60) != string::npos ){

      TSYS_CORR   TS = TS_GPUT;
      t_tsys_corr ts;
     
      ts.a0 = strSci2dbl(line.substr(3,  19)); // cout << " [" + line.substr( 3, 19) + "]";
      ts.a1 = strSci2dbl(line.substr(22, 19)); // cout << " [" + line.substr(22, 16) + "]";
      ts.T  =    str2int(line.substr(41,  9)); // cout << " [" + line.substr(38,  7) + "]";
      ts.W  =    str2int(line.substr(50,  9)); // cout << " [" + line.substr(45,  5) + "]\n";

      _rxnhdr.tsys_corr( TS, ts );

      if( _log ) _log->comment(2,"rinexn","DELTA-UTC: A0,A1,T,W" + tsys_corr2str(TS)
                                                                 + dbl2str(ts.a0*1e9,6)
                                                                 + dbl2str(ts.a1*1e9,6)
                                                           + " " + int2str(ts.T)
                                                           + " " + int2str(ts.W)
                                                           + " " + base_name(_fname)
                               );

    // -------- "CORR TO SYSTEM" --------
    }else if( line.find("CORR TO SYSTEM",60) != string::npos ){      // <= RINEX v2.10
      
      TSYS_CORR   TS = TS_GLUT;
      t_tsys_corr ts;

      ts.a0 = strSci2dbl(line.substr(21, 19));
      ts.a1 = 0.0;
      ts.T  = 0;
      ts.W  = 0;

      _rxnhdr.tsys_corr( TS, ts );

      if( _log ) _log->comment(2,"rinexn","TIME SYSTEM CORR " + tsys_corr2str(TS)
                                                              + dbl2str(ts.a0*1e9,6)
                                                              + dbl2str(ts.a1*1e9,6)
                                                        + " " + int2str(ts.T)
                                                        + " " + int2str(ts.W)
                                                        + " " + base_name(_fname)
                               );

      if( _log ) _log->comment(2,"rinexn","reading CORR TO SYSTEM");


    // -------- "D-UTC" --------
    }else if( line.find("D-UTC A0,A1,T,W,S,U",60) != string::npos ){ // == RINEX v2.11

      TSYS_CORR   TS = TS_SBUT;
      t_tsys_corr ts;

      ts.a0 = strSci2dbl(line.substr(0,  19)); // cout << " [" + line.substr( 0, 19) + "]";
      ts.a1 = strSci2dbl(line.substr(20, 19)); // cout << " [" + line.substr(20, 19) + "]";
      ts.T  =    str2int(line.substr(40,  7)); // cout << " [" + line.substr(40,  7) + "]";
      ts.W  =    str2int(line.substr(48,  5)); // cout << " [" + line.substr(48,  5) + "]\n";

      _rxnhdr.tsys_corr( TS, ts );

      if( _log ) _log->comment(2,"rinexn","D-UTC: A0,A1,T,W,S,U" + tsys_corr2str(TS)
                                                                 + dbl2str(ts.a0*1e9,6)
                                                                 + dbl2str(ts.a1*1e9,6)
                                                           + " " + int2str(ts.T)
                                                           + " " + int2str(ts.W)
                                                           + " " + base_name(_fname)
                               );

    // -------- "LEAP SECONDS" --------
    }else if( line.find("LEAP SECONDS",60) != string::npos ){
      _rxnhdr.leapsec( str2int(line.substr(0,6)) );
      if( _log ) _log->comment(2,"rinexn","reading LEAP SECONDS");


    // -------- "COMMENT" --------
    }else if( line.find("COMMENT",60) != string::npos ){
      _comment.push_back(line.substr(0,60));
      if( _log ) _log->comment(2,"rinexn","reading COMMENT");


    // -------- "END OF HEADER" --------
    }else if( line.find("END OF HEADER",60) != string::npos ){
      _fill_head();
      if( _log ) _log->comment(2,"rinexn","reading END OF HEADER ");
      t_gcoder::_consume(tmpsize);
      _mutex.unlock(); return -1;
    }
    t_gcoder::_consume(tmpsize);
  }

  _mutex.unlock(); return consume;
}

   
/* ----------
 * NAV-RINEX body
 */
int t_rinexn::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();
   
  if( t_gcoder::_add2buffer(buff, sz) == 0 ){ _mutex.unlock(); return 0; };

#ifdef DEBUG
  cout << " BUFFER : \n" << _buffer << "\n size = " << sz  << " END OF BUFFER \n\n"; cout.flush();
#endif
  
  t_gtime epoch;
  int b, e, s, l, i;
  int maxrec = MAX_RINEXN_REC;  // implicite
  string timstr;
  t_gnavdata data = { 0 };

  // RINEX v2.xx   
  if( _version[0] == '2' ){        b = 0; e = 22; l = 19; s = 3; // timstr = "%2s %2d %02d %02d %02d %02d %02d";
    // RINEX v3.xx
  }else if( _version[0] == '3' ){  b = 0; e = 23; l = 19; s = 4; // timstr = "%3s %4d %02d %02d %02d %02d %02d";
    // RINEX
  }else{                           b = 0; e = 23; l = 19; s = 4; // timstr = "%3s %4d %02d %02d %02d %02d %02d";
  }

  string line;
  int consume = 0;
  int tmpsize = 0;
  int recsize = 0;

  while( ( tmpsize = t_gcoder::_getline( line, 0 ) ) >= 0 ){
    
#ifdef DEBUG
    cout << "0: " << line << endl; cout.flush();
#endif

    consume += tmpsize;
    recsize += tmpsize;
    string epostr = line.substr(b,e);
    
    istringstream istr(line.substr(b,e) );
    istr.clear();

    string prn;
    int yr = 0,mn = 0,dd = 0,hr = 0,mi = 0;
    int svn = 0;
    int irc = 0;
    float sec=0.0;
    //char tmpbuff[82];   // RINEX 2 (80) + RINEX 3 (81)

	/*lvhb changed in 20200321 according to novatel decoding empherise*/
	char tmpbuff[83];   // RINEX 2 (80) + RINEX 3 (81)
    int min_sz = 23; // minimum size to be decoded (timestamp)
     
    switch( _version[0] ){
      case '2'  : min_sz = 22; break;  // RINEX 2
      case '3'  : min_sz = 23; break;  // RINEX 3
    }
       
    if( line.size() > 82 || _decode_buffer.size() <= min_sz ){ // avoid decoding such case
//      cout << "skip [" << size() << ":" << tmpsize <<  ":" << recsize <<  ":" << consume << "]:" << line; cout.flush();
      t_gcoder::_consume(tmpsize);
      recsize = consume = 0; break; // read buffer 
    }

    strncpy(tmpbuff, line.c_str(), min_sz );
	tmpbuff[line.size()] = '\0';
    
    if( _version[0] == '2' ){          // don't use '%2i' in scan, but '%2d' instead !
      irc = sscanf(tmpbuff, "%2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %5f", &svn,&yr,&mn,&dd,&hr,&mi,&sec);

      if( irc < 7 ){ // not success - remove this data
//        cout << "fail [" << irc << ":" << size() << ":" << tmpsize <<  ":" << recsize <<  ":" << consume << "]:[" << tmpbuff << "]\n"; cout.flush();
        t_gcoder::_consume(tmpsize);
        recsize = consume = 0; continue;
      }
      prn = t_gsys::eval_sat(svn, t_gsys::char2gsys(_gnsssys));
//        cout << "LINE [" << irc << "]:[" << tmpbuff << "]\n"; cout.flush();

    }else{ 
      istringstream ss(tmpbuff);
      try {
          ss >> prn >> yr >> mn >> dd >> hr >> mi >> sec;
      }
      catch (std::ios_base::failure e) {
          t_gcoder::_consume(tmpsize);
          recsize = consume = 0; continue;
      }
    }
     
//    cout << "LINE IS OK:" << prn << ":" << line; cout.flush();
    shared_ptr<t_gnav> geph = make_shared<t_gnav>();
       
    // !!! SAT musi byt preveden na PRN (pro RINEX ver < 3) --> neni implementovano!
    if(       prn[0] == 'G' ){   maxrec = MAX_RINEXN_REC_GPS;  geph = make_shared<t_gnavgps>();  geph->glog(_log); epoch.tsys(t_gtime::GPS);
    }else if( prn[0] == 'R' ){   maxrec = MAX_RINEXN_REC_GLO;  geph = make_shared<t_gnavglo>();  geph->glog(_log); epoch.tsys(t_gtime::UTC);
    }else if( prn[0] == 'E' ){   maxrec = MAX_RINEXN_REC_GAL;  geph = make_shared<t_gnavgal>();  geph->glog(_log); epoch.tsys(t_gtime::GAL);
    }else if( prn[0] == 'J' ){   maxrec = MAX_RINEXN_REC_QZS;  geph = make_shared<t_gnavqzs>();  geph->glog(_log); epoch.tsys(t_gtime::GPS);   // QZSS is equal to GPS time
    }else if( prn[0] == 'S' ){   maxrec = MAX_RINEXN_REC_SBS;  geph = make_shared<t_gnavsbs>();  geph->glog(_log); epoch.tsys(t_gtime::GPS);   // SBAS is equal to GPS time
    }else if( prn[0] == 'C' ){   maxrec = MAX_RINEXN_REC_BDS;  geph = make_shared<t_gnavbds>();  geph->glog(_log); epoch.tsys(t_gtime::BDS);
    }else if( prn[0] == 'I' ){   maxrec = MAX_RINEXN_REC_IRN;  geph = make_shared<t_gnavirn>();  geph->glog(_log); epoch.tsys(t_gtime::GPS); 
    }else{
      string lg("Warning: not supported satellite satellite system: "+prn); mesg(GWARNING,lg);
      if( _log )  _log->comment(0,"rinexn", lg); 
      else cerr << lg << endl;
      t_gcoder::_consume(tmpsize);
      recsize = consume = 0; continue;
    }

    if (prn[0] == 'R' && str2dbl(_version) > 3.04) maxrec += 4;
    
    epoch.from_ymdhms(yr, mn, dd, hr, mi, sec);

    if( fabs(epoch - _check_dt) > 7*86400 && _check_dt != FIRST_TIME ){
      string lg(prn+" strange epoch ["+epoch.str_ymdhms()+"] or corrupted file ["+base_name(_fname)+"]"); mesg(GWARNING,lg);
      t_gcoder::_consume(tmpsize);
      recsize = consume = 0; continue;
    }
     
     
    if( tmpsize < 57+s ) break;

    data[0] = strSci2dbl(line.substr( 19+s, l ));
    data[1] = strSci2dbl(line.substr( 38+s, l ));
    data[2] = strSci2dbl(line.substr( 57+s, l ));

    i = 2;
    while( i < MAX_RINEXN_REC ){

      // incomplete record
      if( ( tmpsize = t_gcoder::_getline( line, recsize ) ) < 0 ){ break; }

      consume += tmpsize;
      recsize += tmpsize;
      if( ++i < maxrec ){ if( tmpsize>   s ) data[i] = strSci2dbl(line.substr(    s, l )); } //cout << "1 => " << data[i] << "\n"; }
      if( ++i < maxrec ){ if( tmpsize>19+s ) data[i] = strSci2dbl(line.substr( 19+s, l )); } //cout << "2 => " << data[i] << "\n"; }
      if( ++i < maxrec ){ if( tmpsize>38+s ) data[i] = strSci2dbl(line.substr( 38+s, l )); } //cout << "3 => " << data[i] << "\n"; }
      if( ++i < maxrec ){ if( tmpsize>57+s ) data[i] = strSci2dbl(line.substr( 57+s, l )); } //cout << "4 => " << data[i] << "\n"; }

      
      // is record complete and filter-out GNSS systems
      if( geph && i+1 >= maxrec ){
        t_gcoder::_consume(recsize);
        recsize = 0;

        // filter GNSS and SAT
        if( !_filter_gnss(prn) ){
          if( _log ) _log->comment( 4, "rinexn", "skip "+prn);
          break;
        }

        if( epoch < _beg - MAX_NAV_TIMEDIFF || epoch > _end + MAX_NAV_TIMEDIFF){
          if( _log ) _log->comment( 4, "rinexn", "skip "+prn+" "+epoch.str_ymdhms());
          break;
        }

        geph->data2nav( prn, epoch, data );
        geph->gio(_gio_ptr.lock());
	 
	// reset check_dt (already filtered sat)
        if( geph->healthy() ){ _check_dt = geph->epoch(); };

        if( !_filter_gnav(geph, prn) ){
          if( _log ) _log->comment( 4, "rinexn", "skip "+prn+" nav type ["+gnavtype2str(geph->gnavtype())+"]");
          break;
        }

        // collect 1-line messages
        // mesg(GWARNING,geph->linefmt());
        if( _log ) _log->comment(-3,geph->linefmt()+" "+base_name(_fname) );
 
        // fill data
        map<string,t_gdata*>::iterator it = _data.begin();
        while( it != _data.end() ){
          
          if( it->second->id_type()  == t_gdata::ALLNAV ||
              it->second->id_group() == t_gdata::GRP_EPHEM){

            if (prn[0] != 'E' || ((int)data[20] & 0x05) != 0x00)
                ((t_gallnav*)it->second)->add( geph );
            
          }else if( it->second->id_type() == t_gdata::ALLOBJ ){
            t_gallobj* all_obj = (t_gallobj*)it->second;
            shared_ptr<t_gobj>  one_obj = all_obj->obj(geph->sat());
            if( one_obj != 0 ) one_obj->glog(_log);	       	          
            if( geph->id_type() == t_gdata::EPHGLO ){
              int ch = dynamic_pointer_cast<t_gnavglo>(geph)->channel();
              if( one_obj != 0 ) one_obj->channel(ch);
            }
          }
          it++;
        }
        cnt++;
        break;
      }
    }
    if( recsize != 0 ) break; // break the initialization loop if read not finished correctly
  }
  _mutex.unlock(); return consume;
}


/* -------- *
 * ENCODER
 * -------- */

/* ----------
 * NAV-RINEX header
 */
int t_rinexn::encode_head( char* buff, int sz, vector<string>& errmsg)
{ 

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();

  t_gtime tt(t_gtime::UTC);
  t_gtime beg(t_gtime::UTC);
  t_gtime end(t_gtime::UTC);

  // check if data exists
  int count = 0;
  for( auto it = _data.begin(); it != _data.end(); ++it ){
    if( it->second->id_type()  == t_gdata::ALLNAV ||
        it->second->id_group() == t_gdata::GRP_EPHEM )
    {
      count += dynamic_cast<t_gallnav*>(it->second)->vec_nav().size();
      beg    = dynamic_cast<t_gallnav*>(it->second)->beg_gnav("");
      end    = dynamic_cast<t_gallnav*>(it->second)->end_gnav("");
    }
  }

  if( _ss_position == 0 && count > 0 ){ // only if data
      
    if( _version < "3.00" ){
       
      // RINEX 2.x (GNSS-specific)
      if( _sys.size() != 1 ){
        if( _log ) _log->comment(0,"rinexn2","Warning: NAV no encoded, more than one GNSS system identified!");
        _mutex.unlock(); return -1;
      }
      set<string>::iterator itSYS = _sys.begin();
      switch( t_gsys::str2gsys(*itSYS) ){
        case GPS : _ss << "     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE" << endl; break;
        case GLO : _ss << "     2.11           G: GLO NAV DATA                         RINEX VERSION / TYPE" << endl; break;
        case GAL : _ss << "     2.12           N: GNSS NAV DATA    E: GAL NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case BDS : _ss << "     2.12           N: GNSS NAV DATA    C: BDS NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case QZS : _ss << "     2.12           N: GNSS NAV DATA    Q: QZS NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case SBS : _ss << "     2.12           N: GNSS NAV DATA    S: SBS NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case IRN : _ss << "     2.12           N: GNSS NAV DATA    I: IRN NAV DATA     RINEX VERSION / TYPE" << endl; break;
        default :  _ss << "     2.12           N: GPS NAV DATA                         RINEX VERSION / TYPE" << endl; break;
      }
      
    }else{ // RINEX 3.x (GNSS-mixed)
                   _ss << "     3.03           N: GNSS NAV DATA    M (MIXED)           RINEX VERSION / TYPE" << endl;
    }

    //  << "G-Nut/Anubis        GOP/RIGTC           20120618 000000 UTC PGM / RUN BY / DATE"  << endl
    _ss << left << setw(20) << "G-Nut/"+ _pgm 
                            << setw(20) << "GOP/RIGTC"
                            << tt.str("%Y%m%d %H%M%S UTC")                 << " PGM / RUN BY / DATE " << endl
        << "Multi-GNSS (GPS/GLO/GAL/BDS/QZS/SBS/IRN) navigation data   "   << " COMMENT             " << endl
        << "merged from all available IGS-MGEX (and other) files       "   << " COMMENT             " << endl
        << "Contact: gnss@pecny.cz                                     "   << " COMMENT             " << endl;
      
    if( _version >= "3.00" ){
       
      t_rxnhdr rxnhdr = _consolidate_header();

      set<IONO_CORR> iono_list = rxnhdr.iono_corr();
      for( auto it = iono_list.begin(); it != iono_list.end(); ++it ){
        t_iono_corr io = rxnhdr.iono_corr(*it);
        _ss << scientific << setw(4) << iono_corr2str(*it)
            << " " << right
            << setw(12) << setprecision(4) << io.x0
            << setw(12) << setprecision(4) << io.x1
            << setw(12) << setprecision(4) << io.x2
            << setw(12) << setprecision(4) << io.x3
            << left 
            << setw(7)  << " " << setw(20) << "IONOSPHERIC CORR"
            << endl;
      }

      set<TSYS_CORR> tsys_list = rxnhdr.tsys_corr();
      for( auto it = tsys_list.begin(); it != tsys_list.end(); ++it ){
        t_tsys_corr ts = rxnhdr.tsys_corr(*it);
        _ss << scientific << setw(4) << tsys_corr2str(*it)
            << " " << right
            << setw(17) << setprecision(10) << ts.a0
            << setw(16) << setprecision( 9) << ts.a1
            << setw( 7) << setprecision( 0) << ts.T
            << setw( 5) << setprecision( 0) << ts.W
            << left
            << setw(10) << " " << setw(20) << "TIME SYSTEM CORR"
            << endl;
      }
    }

    // LEAP SECONDS + END  (use BEG+END clean_outer to get correct LEAP SEC!)
    if( fabs(end-beg) < 2*86400 ){
      _ss << right << setw(6) << end.leapsec() - TAI_GPS
          <<       "                                                      LEAP SECONDS        " << endl;
    }else{
      cerr << "Warning: RINEXN: NAV BEG/END too far, LEAP SECONDS omitted" << endl;
    }

    _ss << "                                                            END OF HEADER       " << endl;
  }

  int size = _fill_buffer( buff, sz );

  _mutex.unlock(); return size;
}


/* ----------
 * NAV-RINEX data
 */
int t_rinexn::encode_data( char* buff, int sz, int& cnt, vector<string>& errmsg )
{ 

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();

  // RINEX 2.x (GNSS-specific)
  if( _version < "3.00" && _sys.size() != 1 ){
    if( _log ) _log->comment(0,"rinexn2","Warning: NAV no encoded, more than one GNSS system identified!");
    _mutex.unlock(); return -1;
  }

  map<string,int> count_mesg;

  // fill data
  if( _ss_position == 0 ){

    map<string,t_gdata*>::iterator it = _data.begin();
    for( it = _data.begin(); it != _data.end(); ++it ){

      if( it->second->id_type()  == t_gdata::ALLNAV ||
          it->second->id_group() == t_gdata::GRP_EPHEM ){

        vector<shared_ptr<t_geph>> vec_nav = dynamic_cast<t_gallnav*>(it->second)->vec_nav();
        vector<shared_ptr<t_geph>>::iterator itNAV;
        for( itNAV = vec_nav.begin(); itNAV != vec_nav.end(); ++itNAV )
        {
          shared_ptr<t_gnav> gnav = dynamic_pointer_cast<t_gnav>(*itNAV);

          if( gnav->epoch() < _beg || gnav->epoch() > _end ){ continue; } // filter BEG/END

#ifdef DEBUG
  cout << "NAV: " << gnav->sat() << " " << gnav->epoch().str_ymdhms(" ")
	                         << " " << gnav->valid() << endl;
#endif         
          count_mesg[gnav->sat()] += 1;

          t_gnavdata data;
          if( gnav->nav2data( data ) == 0 ){

            if( _version < "3.00" ){
              _ss << setw(2) << str2int(gnav->sat().substr(1,2))
              << " " << gnav->epoch().str("%y %m %d %H %M") // y-digit year, float seconds
              << " " << fixed << setprecision(1) << setw(4) << gnav->epoch().secs() + gnav->epoch().dsec();
            }else{
              _ss << gnav->sat() + " " + gnav->epoch().str("%Y %m %d %H %M %S");  // 4-digit year !!!
            }
            _ss  << fixed << scientific << setprecision(12); // setfill('0')
	     
            for( int i=0; i<gnav->rec(); ++i ){
              if( (i-3)%4 == 0 ){
                
                if( _version < "3.00" )
                _ss << endl << setw(3) << ""; // RINEX 2: 3-CH space before first column first
                else _ss << endl << setw(4) << ""; // RINEX 3: 4-CH space before first column first
                
              }
              _ss << setw(19) << data[i];
            }
            _ss << endl;
          }
        }
      }
    }
  }

  int size = _fill_buffer( buff, sz );

  for( auto it = count_mesg.begin(); it != count_mesg.end(); ++it ){
    if( _log ) _log->comment(1,"rinexn","#NAV "+it->first+": "+int2str(it->second));
  }

  _mutex.unlock(); return size;
}


// fill header information
// ----------
int t_rinexn::_fill_head()
{
  gtrace("t_rinexn::_fill_head");

  int cnt = 0;

  for( auto itDAT = _data.begin(); itDAT != _data.end(); ++itDAT )
  {
    if( itDAT->second->id_type() == t_gdata::ALLOBJ )
    {
      t_gallobj* all_obj = (t_gallobj*)itDAT->second;

      shared_ptr<t_gtrn> obj = dynamic_pointer_cast<t_gtrn>(all_obj->obj(ID_ALLSAT));

      if( obj == 0 ){
        obj = make_shared<t_gtrn>();
        obj->id(ID_ALLSAT);                   
//      obj->name(_fname);
        obj->name(ID_ALLSAT);
        obj->glog(_log);
        obj->header(_rxnhdr, _fname);
        all_obj->add(obj);

        if( _log ) _log->comment(1,"rinex", "Object created: "  +obj->id());
      }else{ 
        if( _log ) _log->comment(1,"rinex", "Object completed: "+obj->id());
        obj->header(_rxnhdr, _fname);        
      }

      ++cnt;
    }
  }

//  cout << "End fill_head\n";

  return cnt;
}

   
// filter out navigation mess. types
// ---------   
bool t_rinexn::_filter_gnav(shared_ptr<t_gnav> geph, const string& prn)
{
  gtrace("t_rinexn::_filter_gnav");
   
  bool ret = true; 
   
  GSYS gs = t_gsys::char2gsys(prn[0]);

   
  // currently only GALILEO can be filtered!
  if( gs == GAL ){

    GNAVTYPE gnav;
    if( _nav[gs].find(gnavtype2str(INAV)) != _nav[gs].end() ){
          gnav = dynamic_pointer_cast<t_gnavgal>(geph)->gnavtype(false); }  // any-source fit in INAV request
    else{ gnav = dynamic_pointer_cast<t_gnavgal>(geph)->gnavtype(true);  }  // exact fit for NAV source

    if( _nav[gs].size() == 0 ||
        _nav[gs].find(gnavtype2str(gnav)) != _nav[gs].end()){ ret = true;  }
    else{                                                     ret = false; }
  }

  return ret; // ALL OTHER SYSTEMS
}


// consolidate header records
// ---------   
t_rxnhdr t_rinexn::_consolidate_header()
{
  gtrace("t_rinexn::_consolidate_header");

  t_rxnhdr rxnhdr;
  
  map<IONO_CORR, map<t_gtetrad, int>> io_sort; // use just 4 records, no support for BDS sat-specific IONO so far!
  map<TSYS_CORR, map<t_gtetrad, int>> ts_sort;

  for( auto itDAT = _data.begin(); itDAT != _data.end(); ++itDAT )
  {
    if( itDAT->second->id_type() != t_gdata::ALLOBJ ){ continue; }

    t_gallobj* all_obj = (t_gallobj*)itDAT->second;

    shared_ptr<t_gtrn> obj = dynamic_pointer_cast<t_gtrn>(all_obj->obj(ID_ALLSAT));

    if( obj != 0 ){
     
      t_gtrn::t_header headers = obj->headers();
      if( _log ) _log->comment(1,"rinexn","Consolidate RINEXN headers: " + int2str(headers.size()));;

      for( auto itHDR = headers.begin(); itHDR != headers.end(); ++itHDR ){
        
        if( !_rxnhdr_all && itHDR != headers.begin() ) break; // just just first HDR

        if( _log ) _log->comment(3,"rinexn","Header: " + itHDR->first );

        t_rxnhdr hdr = itHDR->second;

        // IONO corrections
        set<IONO_CORR> iono_list = hdr.iono_corr();
        for( auto itION = iono_list.begin(); itION != iono_list.end(); ++itION ){
          t_iono_corr io = hdr.iono_corr(*itION);
          t_gtetrad data(io.x0,io.x1,io.x2,io.x3);

          io_sort[*itION][data]++; // just a short cut assuming initialized with zero!
        }

        // TSYS corrections
        set<TSYS_CORR> tsys_list = hdr.tsys_corr();
        for( auto itSYS = tsys_list.begin(); itSYS != tsys_list.end(); ++itSYS ){
          t_tsys_corr ts = hdr.tsys_corr(*itSYS);
          t_gtetrad data(ts.a0,ts.a1,ts.T,ts.W);

          ts_sort[*itSYS][data]++; // just a short cut assuming initialized with zero!
        }
      }
    }else{
      if( _log ) _log->comment(1,"rinexn","Header not found\n");
    }
  }

  t_gtetrad tr_null(0.0,0.0,0.0,0.0);

  // consolidate
  for( auto itION = io_sort.begin(); itION != io_sort.end(); ++itION )
  {
    IONO_CORR io  = itION->first;
         auto itS = io_sort[io].begin();
    for( auto it  = io_sort[io].begin(); it != io_sort[io].end(); ++it )
    { if( it->second > itS->second ) itS = it; }

    if( itS->first == tr_null ){ } // remove empty records
    else{
      t_iono_corr io_corr;
      io_corr.x0 = itS->first[0];
      io_corr.x1 = itS->first[1];
      io_corr.x2 = itS->first[2];
      io_corr.x3 = itS->first[3];
      rxnhdr.iono_corr(io, io_corr);
    }
  }

  for( auto itSYS = ts_sort.begin(); itSYS != ts_sort.end(); ++itSYS)
  {
    TSYS_CORR ts  = itSYS->first;
         auto itS = ts_sort[ts].begin();
    for( auto it  = ts_sort[ts].begin(); it != ts_sort[ts].end(); ++it )
    { if( it->second > itS->second ) itS = it; }

    if( itS->first == tr_null ){ } // remove empty records
    else{
      t_tsys_corr ts_corr;
      ts_corr.a0 = itS->first[0];
      ts_corr.a1 = itS->first[1];
      ts_corr.T  = itS->first[2];
      ts_corr.W  = itS->first[3];
      rxnhdr.tsys_corr(ts, ts_corr);
    }
  }

  return rxnhdr;
}

} // namespace
