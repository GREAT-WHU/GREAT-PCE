/**
*
* @verbatim
    History
    2011-01-10  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       grec.h
* @brief
*.
* @author     JD
* @version    1.0.0
* @date       2011-01-10
*
*/

#ifndef GREC_H
#define GREC_H


#include <stdio.h>
#include <string>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gexport/ExportLibGnut.h"
#include "gdata/gobj.h"
#include "gdata/gdata.h"
#include "gdata/grnxhdr.h"

using namespace std;

namespace gnut {
	/** @brief class for grec. */
	class LibGnut_LIBRARY_EXPORT t_grec : public t_gobj
	{

	public:
		t_grec();
		virtual ~t_grec();
		/** @brief map rec. */
		typedef map<t_gtime, string>      t_maprec;
		/** @brief map header. */
		typedef map<t_gtime, t_rnxhdr>    t_maphdr;

		virtual void addhdr(const t_rnxhdr& hdr, const t_gtime& epo, string path);
		t_maphdr gethdr();
		t_rnxhdr gethdr(const t_gtime& epo);

		t_maprec get_maprec() { return _maprec; }

		void rec(string rec, const t_gtime& beg, const t_gtime& end = LAST_TIME);
		string rec(const t_gtime& t) const;       // set/get receiver  

		void rec_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const;

		virtual bool isrec() { return true; }
		virtual bool istrn() { return false; }

		virtual void compare(shared_ptr<t_grec> grec, const t_gtime& tt, string source);

		virtual vector<t_gtime> rec_id() const;

		void fill_rnxhdr(const t_rnxhdr& rnxhdr);

	protected:
		/** @brief fill data members form rinex header. */
		void _fill_rnxhdr(const t_rnxhdr& rnxhdr);
		/** @brief get one rinex headr. */
		t_rnxhdr _gethdr(const t_gtime& epo);
		/** @brief set receiver name. */
		void   _rec(string rec, const t_gtime& beg, const t_gtime& end = LAST_TIME);
		/** @brief get receiver name (>=t). */
		string _rec(const t_gtime& t) const;

		t_maprec        _maprec;                  // map of receviers
		t_maphdr        _maphdr;                  // map of rinex header information

	private:
	};

} // namespace

#endif
