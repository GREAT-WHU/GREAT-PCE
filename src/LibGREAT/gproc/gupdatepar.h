/**
 * @file         gupdatepar.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GUPDATEPAR_H
#define GUPDATEPAR_H

#include "gall/gallpar.h"
#include "gutils/gtime.h"
#include "gset/gsetbase.h"
#include "gutils/gcycleslip.h"
#include "gproc/glsqmatrix.h"
#include "gmodels/gstochasticmodel.h"

using namespace gnut;

namespace great
{
	/**
	* @brief  class for update parameters information
	*/
	class LibGREAT_LIBRARY_EXPORT t_gupdateparinfo
	{
	public:
		/** @brief default constructor */
		t_gupdateparinfo();
		~t_gupdateparinfo();

		bool exist(int id);
		/** @brief add remove info */
		void add(int id);
		void add(t_gpar newpar);
		/** @brief add update info */
		void add(int id,t_gpar newpar);
		/** @brief add update info with state equ */
		void add(vector<pair<int,t_gpar> > update_par, t_glsqEquationMatrix state_equ);
		/** @brief get remove info */
		void get(vector<int>& remove_id);
		/** @brief get new parlist info */
		void get(vector<t_gpar>& newparlist);
		/** @brief get state equ info */
		void get(vector<t_gpar>& update_newpar,t_glsqEquationMatrix& equ);
		/** @brief get all remove info */
		void get(vector<int>& remove_id,vector<t_gpar>& newparlist,vector<t_gpar>& equ_parlist, t_glsqEquationMatrix& equ);
	


	private:
		set<int>       _remove_id;   // remove par id
		vector<t_gpar> _new_parlist; // new parameters list
		vector<t_gpar> _equ_parlist; // equation list
		t_glsqEquationMatrix _state_equ; // state equation
	};

	/**
	* @brief  class for update parameters
	*/
	class LibGREAT_LIBRARY_EXPORT t_gupdatepar
	{
	public:
		/** @brief default constructor */
		t_gupdatepar();
		virtual ~t_gupdatepar();

		t_gupdatepar(const t_gupdatepar& Other);
		/** @brief set interval */
		void set_interval(double interval);
		/** @brief set time */
		void set_time(const t_gtime& beg, const t_gtime& end, const double& interval) {
			_beg_time = beg;
			_end_time = end;
			_intv = interval;
		}
		/** @brief set iono sigma */
		void set_sig_ion(double sigion);
		/** @brief set cycleslip */
		void set_cycleslip(shared_ptr<t_gcycleslip> cycleslip);
		/** @brief set amb update */
		void set_amb_update_way(bool way);
		/** @brief set system bias */
		void set_sysbias(const SYSBIASMODEL& s, const double& d) { _sys_bias = s; _isb_intv = d; }
		/** @brief set par state mode */
		void set_par_state_mode(par_type type, int order, double dt, double noise);
		/** @brief set band index */
		void set_band_index(const map< GSYS, map<FREQ_SEQ, GOBSBAND> >& band_index) { _band_index = band_index; }
		/** @brief set reference clock */
		void set_ref_clk(const string& s) { _ref_clk = s; }
		/** @brief set bds2 ISB */
		void set_bds2_isb(const bool& b) { _use_bds2_isb = b; }
		/** @brief update amb parameters */
		void update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		/** @brief update ISB parameters */
		void update_isb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		void update_isb_pars_new(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);
		/** @brief get all update parameters */
		virtual t_gupdateparinfo get_all_update_parameters(const t_gtime& epoch,t_gallpar& allpars,const vector<t_gsatdata>& obsdata);

		

	protected:
		/** @brief update amb parameters */
		virtual void _update_amb_pars(const t_gtime& epoch,t_gallpar& allpars,const vector<t_gsatdata>& obsdata,t_gupdateparinfo& update_info) = 0;
		/** @brief update process parameters */
		void _update_process_pars(const t_gtime& epoch, t_gallpar& allpars,t_gupdateparinfo& update_info);
		void _update_process_par(const t_gtime& epoch, int update_id, t_gallpar& allpars, t_gupdateparinfo& update_info);
		/** @brief update state parameters */
		void _update_state_par(const t_gtime& epoch,int update_id, t_gallpar& allpars, t_gupdateparinfo& update_info);
		void _udpate_gps_rec_ifb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info);

	protected:
		
		double _intv    = 0.0;          // interval
		t_gtime _beg_time = FIRST_TIME; // begin time
		t_gtime _end_time = LAST_TIME;  // end time
		double _sig_ion = 9000;         // iono sigma

		int _new_par = 0;
		shared_ptr<t_gcycleslip>   _cycleslip;
		map<par_type, t_statemode> _state_mode;
		bool _lite_update_amb = false;
		map< GSYS, map<FREQ_SEQ, GOBSBAND> >  _band_index;

		SYSBIASMODEL _sys_bias = SYSBIASMODEL::AUTO_WHIT;
		double _isb_intv = 0.0;            // IBS interval
		string _ref_clk;                   // reference clock
		map<string, string> _site_ref_sys;
		map<string, string> _site_bds_ref;
		map<string, string> _sys_ref_site;
		bool _use_bds2_isb = false;

		map<string, double> _bd2_isb;
	};



}


#endif /* GUPDATEPAR_H */