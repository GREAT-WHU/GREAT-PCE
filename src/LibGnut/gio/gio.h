/**
*
* @verbatim
    History
    2011-01-10  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gio.h
* @brief      Purpose: implements gtime class (day and precise time)
*.
* @author     JD
* @version    1.0.0
* @date       2011-01-10
*
*/

#ifndef GIO_H
#define GIO_H


#include <stdio.h>
#include <fstream>
#include <string>

#include "gio/glog.h"
#include "gutils/gmutex.h"
#include "gcoders/gcoder.h"

#define BUF_SIZE 1024

using namespace std;

namespace gnut
{

	class t_gcoder;
	/** @brief class for t_gfile. */
	class LibGnut_LIBRARY_EXPORT t_gio
	{

	public:
		/** @brief default constructor. */
		t_gio();
		virtual ~t_gio();
		/** @brief start reading (could be run in a separate thread). */
		virtual void run_read();
		/** @brief start writing (could be run in a separate thread). */
		virtual void run_write();

		virtual void coder(t_gcoder* coder) { _coder = coder; }
		/**
		* @brief init write.
		* @return
			@retval >0    success
		*/
		virtual int init_write() { return _opened = _init_common(); }
		/**
		* @brief init read.
		* @return
			@retval >0    success
		*/
		virtual int init_read() { return _opened = _init_common(); }
		/**
	    * @brief stop write.
	    * @return
	    *    @retval =0    success
	    */
		virtual int stop_write() { _stop_common(); return 0; }
		/**
		* @brief stop read.
		* @return
		*    @retval =0    success
		*/
		virtual int stop_read() { _stop_common(); return 0; }
		/**
		* @brief set file://dir/name.
		* @param[in]    str        file path
		* @return
		*    @retval =-1 false
		*    @retval =1    true
		*/
		virtual int    path(string str);
		/**
		* @brief get file://dir/name.
		* @return        file path
		*/
		virtual string path()const { return _path; }
		/**
		* @brief set local i/o file.
		* @param[in]    f    file name
		*/
		void   file(const char* f) { _giof.mask(f); }
		/**
		* @brief get local i/o file.
		* @return        file name
		*/
		string file() { return _giof.mask(); }
		/**
		* @brief set verbosity.
		* @param[in]    i    verbosity
		*/
		void   verb(int i) { _verb = i; }
		/**
		* @brief get verbosity.
		* @return        verbosity
		*/
		int    verb() { return _verb; }
		/** @brief stop. */
		void   stop() { _stop = 1; }
		/**
		* @brief get size.
		* @return        size
		*/
		size_t size() { return _size; }
		/**
		* @brief get running.
		* @return        running status
		*/
		int   running() { return _running; }
		/** @brief get opened. */
		int    opened() { return _opened; }
		/** @brief get opened. */
		int connected() { return _opened; }
		/**
		* @brief get glog pointer.
		* @return        glog pointer
		*/
		void glog(t_glog* l) { _log = l; }
		t_glog* glog() { return _log; }
		/** @brief override opreator < == <<. */
		bool            operator<(const t_gio& n) const;
		bool            operator==(const t_gio& n) const;
		friend ostream& operator<<(ostream& os, const t_gio& n);

	protected:
		virtual int _gio_write(const char* buff, int size) = 0;
		virtual int _gio_read(char* buff, int size) = 0;

		virtual int _locf_write(const char* buff, int size);
		virtual int _locf_read(char* buff, int size);

		virtual int _init_common();
		virtual int _stop_common();

		t_glog*         _log;              // log pointer
		int             _fd;               // file descriptor
		size_t          _size;             // buffer size
		string          _path;             // URL-like path (e.g. file:///home/honza/file)
		t_giof          _giof;             // local file
		int             _count;            // record counter
		int             _verb;             // verbosity
		int             _stop;             // require a stop at run() loop
		int             _opened;           // 1: opened/connected
		int             _running;          // running
		t_gcoder*       _coder;            // decoder/encoder
		t_gmutex        _gmutex;

#ifdef BMUTEX
		boost::mutex    _mutex;            // mutual exlusion
#endif

	private:

	};

} // namespace

#endif
