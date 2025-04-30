/**
*
* @verbatim
	History
	2011-10-15  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file      glog.cpp
* @brief
*.  Purpose: implements glog class derived from giof
* @author   JD
* @version  1.0.0
* @date     2011-01-15
*
*/

#include <iostream>
#include <algorithm>

#include "gio/glog.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include <stdarg.h>
using namespace std;

namespace gnut {

// constructor
// ---------
t_glog::t_glog( string mask )
  : t_giof( mask ),
    _time(false),
    _verb(0),
    _size(CACHE_LINES),
	_compare(true)
{}


// destructor
// ---------
t_glog::~t_glog(){}


// time stamp
// ---------
void t_glog::time_stamp(bool b)
{
  _log_gmutex.lock();
  _time = b;
  _log_gmutex.unlock();
}

// time stamp
// ---------
bool t_glog::time_stamp() const
{
  return _time;
}

// cache size
// ---------
void t_glog::cache_size(int i)
{
  _log_gmutex.lock();
  _size = i;
  _log_gmutex.unlock();
}

// get time
// ---------
int t_glog::cache_size() const
{
  return _size;
}

// verbosity
// ---------
void t_glog::verb(int i)
{
  _log_gmutex.lock();
  _verb = i;
  _log_gmutex.unlock();
}

// verbosity
// ---------
int t_glog::verb() const
{
  return _verb;
}

// clear cache
// ---------
void t_glog::clear()
{
  _log_gmutex.lock();
  _cache.clear();
  _log_gmutex.unlock();
}

// compare with _cache or not
// ---------
void t_glog::compare(bool compare)
{
	_log_gmutex.lock();
	_compare = compare;
	_log_gmutex.unlock();
}


int t_glog::_lv2int(const LOG_LV& level)
{
	switch (level)
	{
	case LOG_LV::LOG_ERROR: return 0;
	case LOG_LV::LOG_WARN: return 1;
	case LOG_LV::LOG_INFO: return 0;
	case LOG_LV::LOG_DEBUG: return 1;
	case LOG_LV::LOG_CERR: return 1;
	default:
		break;
	}
	return 999;
}

string t_glog::_lv2str(const LOG_LV& level)
{
	switch (level)
	{
	case LOG_LV::LOG_ERROR: return "ERROR";
	case LOG_LV::LOG_WARN:  return "WARN ";
	case LOG_LV::LOG_INFO:  return "INFO ";
	case LOG_LV::LOG_DEBUG: return "DEBUG";
	case LOG_LV::LOG_CERR:  return "CERR";
	default:
		break;
	}
	return "DEF";
}

string t_glog::_format(const char* fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	int len = vsnprintf(nullptr, 0, fmt, ap);
	va_end(ap);
	std::string buf(len + 1, '\0');
	va_start(ap, fmt);
	vsnprintf(&buf[0], buf.size(), fmt, ap);
	va_end(ap);
	buf.pop_back();
	return buf;
}

   
// comment
// ---------
void t_glog::comment(int l, const string& str)
{ 
	if (fabs(float(l)) > _verb) return;
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_log_gmutex);
#endif
  _log_gmutex.lock();

   string text, text1;
   // negative will not print date !
   if( l < 0 || _time == false ){
     text  = str + "\n";
     text1 = str + "\n"; 
   }else{
     text1 = text = t_gtime::current_time(t_gtime::LOC).str_ymdhms() + " " + str + "\n";

//   text1.replace(9,1," ");     // cache tens of seconds
//   text1.replace(8,2,"  ");    // cache full minutes
     text1.replace(14,5,"     ");    // cache full minutes (PV corrected)
//   text1.replace(6,4,"    ");  // cache tens of minutes
//   text1.replace(5,5,"     "); // cache full hours
   }

   // print repeating message maximally every minut
   if( _compare && find(_cache.begin(), _cache.end(), text1) == _cache.end() ){
       
     this->write( text.c_str(), text.size());
     this->flush();

     // maintain cache
     _cache.push_back(text1);
     if( (int)_cache.size() > _size ) _cache.erase(_cache.begin());
   }

   if (!_compare)
   {
	   this->write(text.c_str(), text.size());
	   this->flush(); 

	   // maintain cache
	   _cache.push_back(text1);
	   if ((int)_cache.size() > _size) _cache.erase(_cache.begin());
   }

   // -------------------------------

   _log_gmutex.unlock();
   return;
}


// comment & identificator
// ---------
void t_glog::comment(int l, const string& ide, const string& str)
{
  comment(l,"[" + ide + ":" + int2str(l) + "] " + str);
}

void t_glog::logInfo(const string& str)
{
    t_glog::comment(LOG_LV::LOG_INFO, "", "", str);
}

void t_glog::logInfo(const string& class_id, const string& func_id, const string& str)
{
    t_glog::comment(LOG_LV::LOG_INFO, class_id, func_id, str);
}

void t_glog::logError(const string& class_id, const string& func_id, const string& str)
{
    t_glog::comment(LOG_LV::LOG_ERROR, class_id, func_id, str);
}

void t_glog::logDebug(const string& class_id, const string& func_id, const string& str)
{
    t_glog::comment(LOG_LV::LOG_DEBUG, class_id, func_id, str);
}

void t_glog::comment(const LOG_LV& level, const string& class_id, const string& func_id, const string& str)
{
	string tmp = "[" + 
    _format("%20s", class_id.substr(0, 18).c_str()) + "::" + 
    _format("%20s",  func_id.substr(0, 18).c_str()) + "::" +
    _lv2str(level) + "] " + str ;
	comment(_lv2int(level),  tmp);
}

} // namespace
