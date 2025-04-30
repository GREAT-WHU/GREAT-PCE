/**
*
* @verbatim
	History
	2011-10-15  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file      glog.h
* @brief
*.  Purpose: implements glog class derived from giof
* @author   JD
* @version  1.0.0
* @date     2011-01-15
*
*/
#ifndef GLOG_H
#define GLOG_H

#include <string>
#include <vector>

#include "gio/giof.h"

#define CACHE_LINES 300

using namespace std;

namespace gnut
{
	/** @brief class for t_glog. */
	class LibGnut_LIBRARY_EXPORT t_glog : public t_giof
	{
	public:
		enum class LOG_LV : int { LOG_ERROR = 0, LOG_WARN = 1, LOG_INFO = 2, LOG_DEBUG = 3, LOG_CERR = 4, LOG_DEF = 999 };
		/** @brief constructor. */
		t_glog(string mask = "");
		virtual ~t_glog();

		void comment(int l, const string& str);
		void comment(int l, const string& ide, const string& str);
		void logInfo(const string& str);
		void logInfo(const string& class_id, const string& func_id, const string& str);
		void logError(const string& class_id, const string& func_id, const string& str);
		void logDebug(const string& class_id, const string& func_id, const string& str);
		void comment(const LOG_LV& level, const string& class_id, const string& func_id, const string& str);

		void time_stamp(bool b);
		bool time_stamp() const;

		void cache_size(int i);
		int  cache_size() const;

		void verb(int i);
		int  verb() const;

		void clear();

		void compare(bool compare);

	protected:
		bool            _compare;          // whether compare with _cache, default true
		bool            _time;             // time stamp
		int             _verb;             // verbosity
		int             _size;             // cache size
		vector<string>  _cache;            // cache for messages
		t_gmutex        _log_gmutex;       // special mutex for comments

#ifdef BMUTEX
		boost::mutex    _log_mutex;
#endif

	private:
		int            _lv2int(const LOG_LV& level);
		string         _lv2str(const LOG_LV& level);
		string         _format(const char* fmt, ...);
	};

} // namespace

#endif
