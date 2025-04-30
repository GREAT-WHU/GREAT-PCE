/**
*
* @verbatim
    History
    2011-10-15  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file      giof.h
* @brief
*.  Purpose: class derived from fstream
*           - file name/mask added
*           - support file hunting according to date/time
*           - automated open/close handling
* @author   JD
* @version  1.0.0
* @date     2011-01-15
*
*/

#ifndef GIOF_H
#define GIOF_H


#include <fstream>
#include <string>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gexport/ExportLibGnut.h"

using namespace std;

namespace gnut
{
	/** @brief class for fstream. */
	class LibGnut_LIBRARY_EXPORT t_giof : public fstream
	{

	public:
		/** @brief constructor 1. */
		t_giof(string mask = "");
		virtual ~t_giof();

		int    mask(string mask);                     // set file mask
		string mask()const { return _mask; }             // get file mask
		string name()const { return _name; }             // get last filename

		int irc()const { return _irc; };                 // get irc status

		void loop(bool l) { _loop = l; }               // set loop read
		bool loop()const { return _loop; }               // get loop read

		void toff(int i) { _toff = i; }                // set time offset [min] for the file name
		int  toff()const { return _toff; }               // get time offset [min]

		virtual int write(const char* buff, int size);        // writting
		int write(const string& s);

		int read(char* buff, int size);               // reading   
		void append(const bool& b = true);               // append mode [false/true]

		void tsys(t_gtime::t_tsys);                   // set/get time system for replacement
		t_gtime::t_tsys tsys();

		string md5sum();                                // get md5sum

	protected:
		string _replace();                              // replace mask to name

		int                _irc;                        // irc status OK=0, Warning>0, Error<0
		string             _mask;                       // original name mask
		string             _name;                       // actual (evaluated) name
		ios::openmode      _omode;                      // output open mode
		bool               _repl;                       // replace if time-specific
		int                _toff;                       // if replace, time offset [min] for the file name
		bool               _loop;                       // loop read
		t_gtime::t_tsys    _tsys;                       // time system for replacement
		t_gmutex           _gmutex;

#ifdef BMUTEX
		boost::mutex       _mutex;                      // mutual exlusion
#endif
	private:

	};

} // namespace

#endif
