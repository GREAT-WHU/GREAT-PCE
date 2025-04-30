/**
*
* @verbatim
History
-1.0 ZhengHJ  2019-09-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		  gallrecover.h
* @brief	  Storage all recover info
* @author     ZhengHJ, Wuhan University
* @version	  1.0.0
* @date		  2019-09-25
*
*/

#ifndef GALLRECOVER_H
#define GALLRECOVER_H

#include "gexport/ExportLibGnut.h"
#include "gdata/grecoverdata.h"
#include "gall/gallsta.h"
#include "gdata/gpoleut1.h"
#include "gall/gallpar.h"
#include <gall/gallprec.h>

namespace great
{

	typedef map<t_gtime, vector<t_grecover_equation*> > t_map_time_equ; ///< time equ
	typedef map<string, t_map_time_equ>	t_map_sat_equ;                  ///< sat equ
	typedef map<string, t_map_sat_equ> t_map_site_equ;                  ///< time equ

	typedef map<t_gtime, vector<t_grecover_par*> > t_map_time_par;      ///< time par
	typedef map<par_type, vector<t_grecover_par*> > t_map_type_par;     ///< type par

	LibGnut_LIBRARY_EXPORT typedef tuple<string, string, string, double, double, t_gtime, t_gtime> t_tuple_ion; // type, rec, sat, value, sigma,beg,end
	/**
	*@brief Class for t_gallrecover
	*/
	class LibGnut_LIBRARY_EXPORT t_gallrecover: public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gallrecover();
		~t_gallrecover();

		void add_allrecover(const t_gallrecover& other);
		/** @brief add recover equation. */
		void add_recover_equation(const t_grecover_equation& recover_equ);
		/** @brief add recover parameter. */
		void add_recover_par(const t_grecover_par& recover_par);
		/** @brief get clk data. */
		void get_clkdata(t_gallprec& clkdata,t_gallprec::clk_type type = t_gallprec::UNDEF);
		/** @brief get iono data. */
		void get_iondata(vector< t_tuple_ion>& ion_data);
		/** @brief get site data. */
		void get_stadata(t_gallsta& stadata);

		vector<t_grecover_equation> get_first_phase_recover_equation(string site, string sat, string freq = "LC");
		bool get_first_phase_recover_equation(string site, string sat, vector<t_grecover_equation>& equ, string freq = "LC");
		/** @brief get recover parameter. */
		vector<t_grecover_par> get_recover_par(par_type parType);
		/** @brief get begin time. */
		t_gtime get_beg_time() const;
		/** @brief get end time. */
		t_gtime get_end_time() const;
		/** @brief get equ begin time. */
		t_gtime get_equ_beg_time() const;
		/** @brief get interval. */
		double get_interval() const;
		/** @brief get sigma0. */
		double get_sigma0() const;
		/** @brief set interval. */
		void set_interval(double intv);
		/** @brief set sigma0. */
		void set_sigma0(double sigma0);

		// get all record
		const vector<t_grecover_data*>& get_all_recover_data() const;
		double get_recover_data_value(par_type parType);
		// get map_par
		const t_map_time_par& get_map_time_par() const;
		const t_map_time_equ& get_map_time_equ() const;
		const t_map_site_equ& get_map_site_equ() const;
		/** @brief get sat list. */
		set<string> get_sat_list() const;
		/** @brief get site list. */
		set<string> get_site_list() const;
		/** @brief get time list. */
		set<t_gtime> get_time_list() const;
		vector<t_gtime> get_all_obs_time_list() const;
		/** @brief get all parameters. */
		t_gallpar get_all_pars();
		map<t_gtime, int> get_sat_number_per_epo() const;

	protected:
		t_grcover_head _recover_head;
		vector<t_grecover_data*> _recover_data;

	private:

		void _add_common_data(t_grecover_data*  data);

		// Mainly for index for time for and sat but no for storaging 
		t_map_time_equ _time_equmap;
		t_map_site_equ _site_sat_time_equmap;
		t_map_time_par _time_parmap;
		map<t_gtime, vector<t_grecover_par> > _time_parmap_pce;

		t_map_type_par _type_parmap;
	};


}


#endif // !GALLRECOVER_H