
/**
* @verbatim
	History
	2012-09-26  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gobs.h
* @brief       Purpose: definition of GNSS observation types
*.
* @author      JD
* @version     1.0.0
* @date        2012-09-26
*
*/

#ifndef GOBS_H
#define GOBS_H

#include <map>
#include <list>
#include <vector>
#include <string>
#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gutils/gnss.h"
#include "gexport/ExportLibGnut.h"

using namespace std;

namespace gnut {

	/** @brief class for t_gattr. */
	class LibGnut_LIBRARY_EXPORT t_gattr {

	public:
		/** @brief constructor 1. */
		t_gattr() { _gattr = ATTR; };
		t_gattr(GOBSATTR a) { _gattr = a; };
		~t_gattr() {};
		/** @brief set attr. */
		virtual void attr(const GOBSATTR& a);
		/** @brief get attr. */
		virtual GOBSATTR attr()const;
		/** @brief override operator ==. */
		virtual bool operator==(const t_gattr& g)const;
		virtual bool valid()const;

	protected:
		GOBSATTR _gattr; ///< gnss attr
	};


	/** @brief class for t_gband derive from t_gattr. */
	class LibGnut_LIBRARY_EXPORT t_gband : public t_gattr {

	public:
		/** @brief constructor 1. */
		t_gband() :t_gattr() { _gband = BAND; };
		t_gband(GOBSBAND b, GOBSATTR a) :t_gattr(a) { _gband = b; };
		virtual ~t_gband() {};

		virtual void band(const GOBSBAND& g);              // set band
		virtual GOBSBAND band()const;                      // get band

		virtual void gattr(const t_gattr& g);              // set t_gattr
		virtual t_gattr gattr()const;                      // get t_gattr

		virtual bool operator==(const t_gband& g)const;
		virtual bool valid()const;

	protected:
		GOBSBAND _gband; ///< gnss band
	};


	/** @brief class for t_gobs derive from t_gattr. */
	class LibGnut_LIBRARY_EXPORT t_gobs : public t_gband {

	public:
		t_gobs() :t_gband() { _gtype = TYPE; };
		t_gobs(GOBSTYPE t, GOBSBAND b, GOBSATTR a) :t_gband(b, a) { _gtype = t; };
		t_gobs(GOBS g) { gobs(g); };
		virtual ~t_gobs() {};

		virtual void type(const GOBSTYPE& t);               // set type (only, inherit band&attr!)  
		virtual GOBSTYPE type()const;                       // get type

		virtual void gband(const t_gband& g);               // set attr
		virtual t_gband gband()const;                       // get attr 

		int gobs(const GOBS& g);                            // set type (only! inherit)
		int gobs(const string& s);

		GOBS gobs()const;                                   // get gobs enum
		GOBS gobs2CH(GSYS gs)const;                         // get 2char gobs (only C1/C2 and P1/P2 and L1/L2)
		GOBS gobs3CH()const;                                // get 3char gobs ( )

		void gobs2to3(GSYS gs);							  // change obs from 2 to 3

		bool operator==(const t_gobs& g)const;

		bool valid()const;
		bool is_code()const;
		bool is_phase()const;
		bool is_doppler()const;

	protected:

		GOBSTYPE          _gtype; ///< gtype
	};

} // namespace

#endif // GOBS_H

