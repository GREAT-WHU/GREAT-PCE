/**
*
* @verbatim
     (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

     (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
        Research Institute of Geodesy, Topography and Cartography
        Ondrejov 244, 251 65, Czech Republic
  @endverbatim
*
* @file        gsetout.cpp
* @brief       implements output setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gset/gsetout.h"
#include "gutils/gfileconv.h"

using namespace std;
using namespace pugi;

namespace gnut {


// Convertor for OUT formats
// ----------
OFMT t_gsetout::str2ofmt(const string& s)
{ 
  string tmp = s;
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  if( tmp == "OUT"     )  return XXX_OUT;
  if( tmp == "LOG"     )  return LOG_OUT;
  if( tmp == "PPP"     )  return PPP_OUT;
  if( tmp == "FLT"     )  return FLT_OUT;
  if (tmp == "RECCLK")	 return RECCLK_OUT;
  if (tmp == "SATCLK")	 return SATCLK_OUT;
  if (tmp == "CLK")      return CLK_OUT;
  if (tmp == "EPODIR") return EPODIR_OUT;
  return OFMT(-1);
}

// Convertor for OUT formats
// ----------
string t_gsetout::ofmt2str(const OFMT& f)
{ 
  switch ( f ){        
   case XXX_OUT:     return "OUT";  
   case LOG_OUT:     return "LOG";
   case PPP_OUT:     return "PPP";
   case FLT_OUT:     return "FLT";
   case RECCLK_OUT:	 return "RECCLK";
   case SATCLK_OUT:	 return "SATCLK";
   case CLK_OUT:     return "CLK";
   case EPODIR_OUT: return "EPODIR";
   default:          return "UNDEF";
  }
  return "UNDEF";
}


// Constructor
// ----------
t_gsetout::t_gsetout()
  : t_gsetbase(),
     _append(false),
     _verb(0),
     _upd(0),
     _len(0),
     _smp(0)
{
   _set.insert(XMLKEY_OUT);
}


// Destructor
// ----------
t_gsetout::~t_gsetout()
{}


// Get formats output size
// ----------
int t_gsetout::output_size(const string& fmt)
{
  _gmutex.lock();
   
  int tmp = _outputs(fmt).size();

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
int t_gsetout::verb()
{
  _gmutex.lock();

  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("verb").as_int();

  _gmutex.unlock(); return tmp;
}


// Retrun value
// -----------
bool t_gsetout::append()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("append").as_bool();

  _gmutex.unlock(); return tmp;
}

   
// Get string outputs
// ----------
string t_gsetout::outputs(const string& fmt)
{ 
  _gmutex.lock();
  
  string tmp = _outputs(fmt);

  _gmutex.unlock(); return tmp;
}

int t_gsetout::sp3_obslimit()
{
	_gmutex.lock();
	int obs;
	obs = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).attribute("obslimit").as_int();
	
	_gmutex.unlock(); return obs;
}
   
// get string output version
// ---------
string t_gsetout::version(const string& fmt)
{
  _gmutex.lock();
   
  string ver = DEFAULT_FILE_VER;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
  if( ! fmt.empty() &&
      ! node.attribute("ver").empty() )
  {
    ver = node.attribute("ver").as_string();
  }

  _gmutex.unlock(); return ver;
}


   
// get period [min] for output file content
// ----------
int t_gsetout::out_length(const string& fmt)
{
  _gmutex.lock();
   
  int len = DEFAULT_FILE_LEN;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());

  if( ! fmt.empty() &&
      ! node.attribute("len").empty() )
  {
    len = node.attribute("len").as_int();
  }

  _gmutex.unlock(); return len;
}

   
// get sample [sec] for output file data sampling
// ---------
float t_gsetout::out_sample(const string& fmt)
{
  _gmutex.lock();
   
  float smp = DEFAULT_FILE_SMP;
  xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).child(fmt.c_str());
  if( ! fmt.empty() &&
      ! node.attribute("smp").empty() )
  {
    smp = node.attribute("smp").as_float();
  }
   
  _gmutex.unlock(); return smp;
}

// Get output formats
// ----------
set<string> t_gsetout::oformats()
{
  return _oformats();
}


// Get output formats
// ----------
set<string> t_gsetout::_oformats()
{
  set<string> tmp;
  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).first_child(); node; node = node.next_sibling() ){
    tmp.insert( node.name() );
  }
  return tmp;
}


// Get formats outputs
// ----------
string t_gsetout::_outputs(const string& fmt)
{
  string str;

  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).first_child(); node; node = node.next_sibling() ){
    if( node.name() == fmt ){
      istringstream is( node.child_value() );
      while( is >> str && !is.fail() ){
        if( str.find("://") == string::npos ) str = GFILE_PREFIX + str;
        return str;
      }       
    }
  }
  return "";
}

// settings check
// ----------
void t_gsetout::check()
{
  this->_upd_glog(); // update log!
   
  _gmutex.lock();
   
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_OUT);

  // check existence of attributes
  _default_attr(node,"append", _append);

  // check supported input formats (see OFMT enum !)
  set<string> ofmt = _oformats();
  set<string>::const_iterator itFMT = ofmt.begin();
  while( itFMT != ofmt.end() ){
    string fmt = *itFMT;
    OFMT  ofmt = str2ofmt( fmt );
    if( ofmt < 0 ){
      _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).remove_child( node.child(fmt.c_str()) );
      if( _log ) _log->comment(0, "Warning: " + fmt + " out format not implemented [gsetout::check()]!");
      else                cerr << "Warning: " + fmt + " out format not implemented [gsetout::check()]!\n";
      itFMT++;
      continue;
    }
  
    // check application-specific output format
    if( _OFMT_supported.find( ofmt ) == _OFMT_supported.end() ){
      _doc.child(XMLKEY_ROOT).child(XMLKEY_OUT).remove_child( node.child(fmt.c_str()) );
      if( _log ) _log->comment(0, "Warning: " + fmt + " out format not supported by this application!");
      else                cerr << "Warning: " + fmt + " out format not supported by this application!\n";       
    }     
    itFMT++;            
  }  
   
  _gmutex.unlock(); return;
   
  this->_upd_glog(); // update log!   
}


// settings help
// ----------
void t_gsetout::help()
{
  _gmutex.lock();
       
  cerr << " <outputs append=\""  << _append << "\" verb=\"" << _verb << "\" >\n"
       << "   <flt> file://dir/name </flt>    \t\t <!-- filter output encoder -->\n"
       << " </outputs>\n";
   
  cerr << "\t<!-- outputs description:\n"
       << "\t <encoder> path </encoder>\n"
       << "\t ... \n"
       << "\t where path contains [file,tcp,ntrip]:// depending on the application\n"
       << "\t -->\n\n";

   _gmutex.unlock(); return;
}


// Read glog from settings
void t_gsetout::_upd_glog()
{
  if( _log ){ // only if exists
    string logstr( this->outputs("log") );

    if( !logstr.empty() ){
      _log->mask(logstr);
      _log->append( this->append() );
      _log->verb( this->verb() );
    }
  //  int verb = this->verb();
  //  if( verb > _log->verb() ) _log->verb( verb ); // only if higher verbosity
  }   
}


} // namespace
