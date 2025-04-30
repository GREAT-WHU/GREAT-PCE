/**
 * @file         gpcelsqIF.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GPCELSQIF_H
#define GPCELSQIF_H

#include "gmodels/glsqprocIF.h"
#include "gexport/ExportLibGREAT.h"
#include "gcoders/sp3.h"

namespace great
{

	/**
	* @brief main control class for PCE
	*/
	class LibGREAT_LIBRARY_EXPORT t_gpcelsqIF :public t_glsqproc
	{
	public:
		using t_all_res = map<pair<string, string>, double>;
		using t_B_mat = vector<pair<int, double>>;
		/**@brief constructor */
		t_gpcelsqIF(t_gsetbase* set, t_gallproc* data, t_glog* log);
		virtual ~t_gpcelsqIF();
		/** @brief init process data and parameter */
		bool InitProc(t_gallproc* data, const t_gtime& beg, const t_gtime& end);
		/** @brief process batch */
		bool ProcessBatch(t_gallproc* data, const t_gtime&beg, const t_gtime& end) override;
		/** @brief generate product */
		bool GenerateProduct() override;

	protected:
		/** @brief init lsq product data */
		bool _initLsqProdData(t_gallprod* data) override;
		/** @brief init lsq process pars */
		bool _initLsqProcPars(t_glsq* lsq)      override;
		/** @brief init lsq process data */
		bool _initLsqProcData(t_gallproc* data) override;
		/** @brief process one receiver and one satellite data */
		bool _processOneRecOneSat(const t_gtime& crt_epoch, const std::string& rec,
			const std::string& sat,
			t_gsatdata& crd_obs) override;
		/** @brief generate solved clock product */
		void _print_clk_solved();
		/** @brief init one epoch data */
		bool _initOneEpoch();
		/** @brief process one epoch data */
		bool _processOneEpoch(const t_gtime& crt_epoch, vector<t_gsatdata>& crt_obs) override;
		/** @brief process one receiver data */
		bool _processOneRec_thread_safe(const t_gtime& crt_epoch, const string& crt_rec, vector<t_gsatdata>& crt_obs, t_glsqEquationMatrix& equ_result);
		/** @brief solve one epoch equation */
		bool _solveEpoch(vector<t_gsatdata>& crt_obs);
		/** @brief whether reference clock has observations in current epoch */
		bool _ref_clk_valid(t_glsq* lsq);
		/** @brief check and reset the reference clock */
		bool _check_ref_clk(t_glsq* lsq);
		/** @brief get current satellite clocks */
		void _get_clk_crt(t_glsq* lsq, bool update_std=true);
		/** @brief get current obs num */
		void _get_obs_crt_num();
		/** @brief find satellite with maximum observations */
		string _sat_obs_max(t_glsq* lsq);

		shared_ptr<t_gqualitycontrol>  _quality_control = nullptr;
		t_gio* _gioout = nullptr;

		string _ref_clk;     //reference clock
		string _ref_clk_crt;

		double _crt_mjd = 0; // cunrrent mjd time
		chrono::high_resolution_clock::time_point _beg_epo_time;
		double _prepare_time = 0.0;
		
		map<string, t_glsqEquationMatrix> _map_all_equ; // map of all equation
		set<string> _bad_site;

		// store solved clock values 
		map<string, double> _clk_crt;
		map<string, double> _clk_std_crt;
		double _crt_sigma = 0.0;
		bool _epo_solved = false;
	};
}

#endif /*  GPCELSQIF_H */