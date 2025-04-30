/**
*
* @verbatim
History
-1.0 JD  2011-02-14  creat the file.
-1.1 JD  2018-08-13  update the file
*
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     gallnav.h
* @brief    container for all navigation systems
*
*
* @author   PV
* @version  1.1.0
* @date     2018-08-13
*
*/
#ifndef GALLNAV_H
#define GALLNAV_H


#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include "gexport/ExportLibGnut.h"
#include "gdata/gdata.h"
#include "gdata/geph.h"
#include "gdata/gnavgps.h"
#include "gdata/gnavglo.h"
#include "gdata/gnavgal.h"
#include "gdata/gnavqzs.h"
#include "gdata/gnavsbs.h"
#include "gdata/gnavbds.h"
#include "gdata/gnavirn.h"
#include "gdata/grxnhdr.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtetrad.h"

#define MAX_GPS_PRN 32
#define MAX_GLO_PRN 24
#define MAX_GAL_PRN 30
#define NAV_BUF   1024


using namespace std;

namespace gnut
{
	/**
	 *@brief Class for navigation system setting derive from t_gdata
	 */
	class LibGnut_LIBRARY_EXPORT t_gallnav : public t_gdata
	{

	public:
		/** @brief constructor. */
		t_gallnav();
		virtual ~t_gallnav();

		typedef multimap<t_gtime, shared_ptr<t_geph> > t_map_ref; // all data for a single satellite
		typedef map<string, t_map_ref>                 t_map_sat; // all data for all satellites

		void chk_health(bool b) { _chk_health = b; }        // set nav healthy status control
		void chk_navig(bool b) { _chk_navig = b; }          // set nav internal quality control
		void chk_tot(bool b) { _chk_tot = b; }             // test if tot < t
		/**
		* @brief get satellite health.
		*
		* @param[in]  sat    satellite
		* @param[in]  t        time
		* @return    status
		*/
		virtual bool health(string sat, const t_gtime& t);       // get satellite health
		/**
		 * @brief
		 *
		 * @param sat
		 * @param t
		 * @param xyz
		 * @param var
		 * @param vel
		 * @param chk_mask
		 * @return int
		 */
		virtual int nav(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true); // [m]
		/**
		 * @brief
		 *
		 * @param sat
		 * @param t
		 * @param xyz
		 * @param var
		 * @param vel
		 * @param chk_mask
		 * @return int
		 */
		virtual int pos(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true); // [m]
		/**
		* @brief return clock corrections.
		*
		* @param[in]  sat        satellite
		* @param[in]  t            time
		* @param[in]  clk        clock offset
		* @param[in]  var
		* @param[in]  dclk        difference of clock offset
		* @param[in]  chk_mask
		* @return    irc
		*/
		virtual int clk(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL, bool chk_mask = true); // [s]  
		/**
		* @brief print function.
		*/
		virtual void print(string sat, const t_gtime& t);
		/**
		* @brief clean invalid messages.
		*/
		virtual void clean_invalid();
		/**
		* @brief clean duplicit messages.
		*/
		virtual void clean_duplicit();
		/**
		* @brief erase data.
		*/
		void erase_time(const t_gtime& t);
		/**
		* @brief clean function.
		*
		* @param[in]  beg        begin time
		* @param[in]  end        end time
		* @return    void
		*/
		virtual void clean_outer(const t_gtime& beg = FIRST_TIME,
			const t_gtime& end = LAST_TIME);

		virtual t_gtime beg_gnav(string prn = "");                            // get first t_gnav epoch
		virtual t_gtime end_gnav(string prn = "");                            // get last  t_gnav epoch
		// IMPROVE beg_time/end_time to distinguish GALLNAV/GALLPREC - t_satview !
		virtual t_gtime beg_time(string prn = "") { return beg_gnav(prn); }    // get first t_gnav epoch
		virtual t_gtime end_time(string prn = "") { return end_gnav(prn); }    // get first t_gnav epoch
		// =========================================================
		/**
        * @brief get all satellites.
        *
        * @return    all satellites
        */
		virtual set<string> satellites();
		/**
		* @brief get all systems.
		*
		* @return    all systems
		*/
		virtual set<GSYS>   systems();
		/**
		* @brief add single navigation message.
		*
		* @param[in]  nav    navigation system
		* @return    0
		*/
		virtual int add(shared_ptr<t_gnav> nav);
		/**
		* @brief get number of epochs.
		*
		* @param[in]  prn    satellite prn
		* @return    number of epochs
		*/
		virtual unsigned int nepochs(const string& prn);
		/**
		* @brief get list of nav epochs.
		*
		* @param[in]  prn    satellite prn
		* @return    list of nav epochs
		*/
		virtual set<t_gtime> epochs(string prn = "");
		/**
		* @brief get list of nav messages.
		*
		* @param[in]  prn    satellite prn
		* @return    list of nav messages
		*/
		virtual vector<shared_ptr<t_geph> > vec_nav(string prn = "");
		/**
		* @brief get list of calculated crd.
		*
		* @param[in]  prn    satellite prn
		* @param[in]  beg    begin time
		* @return    list of calculated crd
		*/
		virtual map<string, t_gtriple>  map_xyz(set<string> prns, const t_gtime& beg);
		/**
		* @brief get list of calculated clk.
		*
		* @param[in]  prn    satellite prn
		* @param[in]  beg    begin time
		* @return    list of calculated clk
		*/
		virtual map<string, double>     map_clk(set<string> prns, const t_gtime& beg);

		/**
		* @brief multimap mode
		*
		* @return    multimap
		*/
		virtual void multi(bool b) { _multimap = b; }
		virtual bool multi() { return _multimap; }
		/**
		* @brief overwrite mode
		*
		* @return    overwrite
		*/
		virtual void overwrite(bool b) { _overwrite = b; }
		virtual bool overwrite() { return _overwrite; }

		/**
		* @brief position/clock reference point.
		*
		* @return    com
		*/
		virtual void com(bool b) { _com = b; }
		virtual bool com() { return _com; }
		/**
		* @brief offset for satdata.
		*
		* @return    offset
		*/
		virtual void offset(int i) { _offset = i; }
		virtual int  offset() { return _offset; }

		/**
		* @brief get number of satellites.
		*
		* @param[in]  gs        navigation system
		* @return    number of satellites
		*/
		virtual int nsat(GSYS gs);                                         // get number of satellites
		virtual int intv(GSYS gs);                                         // get interval between messages
		virtual int have(GSYS gs, const t_gtime& beg, const t_gtime& end); // get existing number of messages
		virtual int expt(GSYS gs, const t_gtime& beg, const t_gtime& end); // get expected number of messages
		virtual int excl(GSYS gs, const t_gtime& beg, const t_gtime& end); // get excluded number of messages

		virtual int consolidate(double cfdi = 0.0);                       // consolidate (confident interval, 0 = auto)
		virtual int consolidate_others();                                  // consolidate healthy status & biases

		/**
		* @brief find appropriate t_geph element (interface only).
		*
		* @param[in]  sat        satellite
		* @param[in]  t            time
		* @param[in]  chk_mask
		* @return    tmp
		*/
		virtual shared_ptr<t_geph> find(string sat, const t_gtime& t, bool chk_mask = true);   // find appropriate t_geph element (interface only)
		/**
		* @brief find appropriate t_geph elements (interface only)
		*
		* @param[in]  sat        satellite
		* @param[in]  t            time
		* @return    tmp
		*/
		virtual vector<shared_ptr<t_geph> > find_mult(string sat, const t_gtime& t);
		/**
		* @brief return frequency number of GLONASS.
		*
		* @return    _glo_freq_num
		*/
		map<string, int> glo_freq_num() { return _glo_freq_num; }
		/**
		* @brief get ionosphere correction.
		*
		* @param[in]  c        ionosphere correction
		* @return    io
		*/
		t_iono_corr get_iono_corr(const IONO_CORR c) const;
		/**
		* @brief add ionosphere correction.
		*
		* @param[in]  c        ionosphere correction
		* @param[in]  io        ionosphere correction
		* @return    void
		*/
		void add_iono_corr(const IONO_CORR c, t_iono_corr io);

	protected:
		/**
		* @brief find appropriate t_geph element.
		*
		* @param[in]  sat        satellite
		* @param[in]  t            time
		* @param[in]  chk_mask
		* @return    null
		*/
		virtual shared_ptr<t_geph> _find(string sat, const t_gtime& t, bool chk_mask = true);
		/**
		* @brief find vector of appropriate t_geph elements
		*
		* @param[in]  sat        satellite
		* @param[in]  t            time
		* @return    vector of appropriate t_geph elements
		*/
		vector<shared_ptr<t_geph> > _find_mult(string sat, const t_gtime& t);
		/**
		* @brief find appropriate t_geph element.
		*
		* @param[in]  sat        satellite
		* @param[in]  iod        issue of data
		* @param[in]  t            time
		* @param[in]  chk_mask
		* @return
		*/																
		virtual shared_ptr<t_geph> _find(string sat, int iod, const t_gtime& t);

		bool               _com;         // position/clock reference point (com = true; apc = false);
		int                _offset;      // offset for RTCM corrections
		int                _nepoch;      // maximum number of epochs (0 = keep all)
		t_map_sat          _mapsat;      // map over all satellites (positions,..?)
		bool               _multimap;    // use multimap for redundant records
		bool               _overwrite;   // overwrite mode (for derived classes with MAP)
		bool               _chk_health;  // check satellite health (navigation)
		bool               _chk_navig;   // check navigation messages (internal)
		shared_ptr<t_geph> _null;        // null pointer
		bool _chk_tot;                   // check tot
		map<string, int> _glo_freq_num;  // frequency number of GLONASS
		t_map_iono       _brdc_iono_cor; // brdc iono correction

	};

} // namespace

#endif
