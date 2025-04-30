/**
 * @file         glsqproc.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gmodels/glsqproc.h"
#include "gset/gsetpar.h"
#include "gset/gsetamb.h"
#include "gset/gsetout.h"
#include "gset/gsetturboedit.h"
#include "gio/gfile.h"
#include "gutils/ginfolog.h"
#include "gutils/gcycleslip.h"
#include "gutils/gturboedit.h"
#include "gutils/gtime.h"
#include "gcoders/recover.h"
#include "gcoders/rinexc.h"
#include "gcoders/poleut1.h"
#include "gcoders/sp3.h"
#include "gutils/gfileconv.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gmodels/glsqprocIF.h"
#include "gproc/gupdateparIF.h"
#include "gutils/ginfolog.h"
#include "gproc/gupdateparALL.h"
#include "gcoders/biabernese.h"
#include "gutils/gstring.h"
#include <algorithm>
#include <thread>
#include <sstream>


#if  defined _WIN32 || defined _WIN64
#include <io.h>
#include <direct.h> 
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#include <unistd.h>
#include <sys/stat.h>
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

namespace great
{
	t_glsqproc::t_glsqproc()
	{
	}

	t_glsqproc::t_glsqproc(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		_glog(log),
		_gset(set),
		_gall_proc(data),
		_recclk_threshold(1e-8)
	{
		// check log
		if (_glog == nullptr)
		{
			_myLog = _create_log(_glog);
		}

		
		// get settings from xml
		_init_xml_settings(set);
		_init_lsq_settings(set, data, log);

	}

	t_glsqproc::~t_glsqproc()
	{
		if (_lsq)
		{
			delete _lsq;
			_lsq = nullptr;
		}

		if (_myLog)
		{
			_delete_log(_glog);
		}
	}

	void t_glsqproc::add_coder(const vector<t_gcoder*>& coder)
	{
		_gcoder = coder;
	}

	bool t_glsqproc::ProcessBatch(t_gallproc* data, const t_gtime& beg, const t_gtime& end)
	{
		return true;
	}

	bool t_glsqproc::GenerateProduct()
	{
		return false;
	}

	bool t_glsqproc::_init_rec_trop_pars(t_glsq* lsq)
	{
		if (_rec_list.empty() || !_lsq || !_gset) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{
			// init the ztd
			string rec = *it;
			t_gpar par_ztd(rec, par_type::TRP, ++ipar_num, "");
			shared_ptr<t_gobj> _grec = _gall_obj->obj(rec);
			if (_grec->id_type() == t_gdata::REC_LEO)
			{
				//it++;
				continue;
			}
			double sig_ztd = dynamic_cast<t_gsetpar*>(_gset)->sigZtd(rec);
			double sig_TropPd = dynamic_cast<t_gsetpar*>(_gset)->sigTropPd(rec);
			par_ztd.apriori(sig_ztd);
			par_ztd.value(0.0);
			par_ztd.setTime(_beg_time, _beg_time + _pwc_intv * 60.0);
			lsq->add_parameter(par_ztd);
			if (_ztd_model == ZTDMODEL::PWC)
			{
				lsq->add_par_state_equ(par_ztd.parType, 1, _pwc_intv / 60.0, sig_TropPd);
			}
			if (_ztd_model == ZTDMODEL::STO)
			{
				lsq->add_par_state_equ(par_ztd.parType, 1, _obs_intv / 3600.0, sig_TropPd);
			}

		}
		return true;
	}

	bool t_glsqproc::_init_rec_clk_pars(t_glsq* lsq)
	{

		if (_rec_list.empty() || !lsq || !_gset)
		{
			return false;
		}

		string site_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		auto it_beg = _rec_list.begin();
		auto it_end = _rec_list.end();
		for (auto it = it_beg; it != it_end; it++)
		{

			int ipar_num = lsq->_x_solve.parNumber();

			string site = *it;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigRecCLK(site);
			double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();

			if (site == site_ref)
			{
				if (sig_refclk == 0.0)
				{
					//it++;
					continue;
				}
				else
				{
					sigclk = sig_refclk;
				}
			}

			//intit the Clock
			t_gpar par_clk(site, par_type::CLK, ++ipar_num, "");
			par_clk.value(0.0);
			par_clk.setTime(_beg_time, _beg_time);
			par_clk.apriori(sigclk);
			lsq->add_parameter(par_clk);
		}
		return true;
	}

	bool t_glsqproc::_init_sat_clk_pars(t_glsq* lsq)
	{
		if (!lsq || !_gset || _sat_list.empty()) return false;

		int ipar_num = _lsq->_x_solve.parNumber();
		auto it_beg = _sat_list.begin();
		auto it_end = _sat_list.end();
		string sat_ref = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		double sig_refclk = dynamic_cast<t_gsetproc*>(_gset)->sig_ref_clk();
		
		for (auto iter = it_beg; iter != it_end; iter++)
		{
			string satID = *iter;
			double sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			if (satID == sat_ref && sig_refclk == 0.0)
			{
				continue;
			}
			else if (satID == sat_ref && sig_refclk != 0.0)
			{
				sigclk = sig_refclk;
			}
			else
			{
				sigclk = dynamic_cast<t_gsetpar*>(_gset)->sigSatCLK(satID);
			}

			t_gpar par_satclk("", par_type::CLK_SAT, ++ipar_num, satID);
			par_satclk.setTime(_beg_time, _beg_time);
			par_satclk.apriori(sigclk);
			par_satclk.value(0.0);
			lsq->add_parameter(par_satclk);
		}
		return true;
	}
	
	bool t_glsqproc::_initLsqProdData(t_gallprod* data)
	{
		return false;
	}

	bool t_glsqproc::_initLsqProcPars(t_glsq* lsq)
	{
		return false;
	}

	bool t_glsqproc::_initLsqProcData(t_gallproc* data)
	{
		return false;
	}
	
	bool t_glsqproc::_select_obs(const t_gtime& epoch, vector<t_gsatdata>& all_obs)
	{
		if (!_slip12)
		{
			_glog->logError("t_glsqproc", "_select_obs", "_slip 12 is nullptr");
			return false;
		}

		// ===================================================================================================================
		// set tb log in satdata
		for (auto& iter : all_obs)
		{
			const string& obs_rec = iter.site();
			const string& obs_sat = iter.sat();
			const GSYS&   obs_sys = iter.gsys();

			if (_frequency < 2) return false;
			bool obs_12 = _slip12->use_of_obs(obs_rec, obs_sat, epoch) && (_freq_number[obs_sys] >= 2); iter.tb12(obs_12);
		}
		return true;
	}

	bool t_glsqproc::_select_rec_obs(const string& rec, const t_gtime& epoch, vector<t_gsatdata>& all_obs, vector<t_gsatdata>& rec_obs)
	{		
		for (auto& iter : all_obs)
		{
			const string& obs_rec = iter.site();
			const string& obs_sat = iter.sat();

			if (obs_rec != rec) continue;

			bool obs_valid = true;
			if (_lite_turboedit)
			{
				set<GOBSBAND> band_avail = iter.band_avail();
				bool freq1_ok = false;
				for (auto itband : band_avail)
				{
					FREQ_SEQ freq2 = _freq_index[iter.gsys()][itband];
					switch (freq2)
					{
					case FREQ_SEQ::FREQ_1:
						freq1_ok = true;
						break;
					case FREQ_SEQ::FREQ_2:
						iter.tb12(true);
						break;
					default:
						break;
					}
				}
				if (!freq1_ok)
				{
					iter.tb12(false);
				}
				
			}

			if (iter.tb12())
			{
				rec_obs.push_back(iter);
			}
		}

		// ===================================================================================================================
		if (rec_obs.empty())
		{
			_glog->logDebug("t_glsqproc", "_select_rec_obs", "This epoch have no useful data : " + rec);
			return false;
		}

		return true;
	}

	bool t_glsqproc::_processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs)
	{
		// ========================================================================================================================================
		// check obs size()	
		if (crt_obs.empty()) {
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("crt obs is empty"));
			return false;
		}
		
		// ========================================================================================================================================
		for (auto iter = crt_obs.begin(); iter != crt_obs.end();) {
			if (iter->sys() == "R" && iter->channel() >= DEF_CHANNEL) {
				if (_glofrq_num.find(iter->sat()) != _glofrq_num.end()) {
					iter->channel(_glofrq_num.at(iter->sat()));
				}
				else {
					_glog->logInfo("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful frequency id for: " + iter->sat()));
					iter = crt_obs.erase(iter);
					continue;
				}
			}
			++iter;
		}
		chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
		// ========================================================================================================================================
		bool select_obs = _select_obs(crt_epoch, crt_obs);
		if (!select_obs)
		{
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("select_obs failed"));
			return false;
		}
		
		map<string, vector<t_gsatdata> > map_site_obs;
		vector<t_gsatdata> crt_obs_new;
		for (const auto& crt_rec : _rec_list) {
			vector<t_gsatdata> crt_rec_obs;
			// get all the obs for the rec
			bool select_valid = _select_rec_obs(crt_rec, crt_epoch, crt_obs, crt_rec_obs);
			if (!select_valid) {
				_glog->logInfo("t_gpcelsqIF", "_processOneEpoch", crt_epoch.str_mjdsod("no useful data : " + crt_rec));
				continue;
			}
			if (_crd_est != CONSTRPAR::KIN) _quality->processOneEpoch(crt_epoch, crt_rec, _rec_crds[crt_rec], crt_rec_obs);
			map_site_obs[crt_rec] = crt_rec_obs;
			crt_obs_new.insert(crt_obs_new.end(), crt_rec_obs.begin(), crt_rec_obs.end());
		}

		bool updata_valid = _lsq->update_parameter(crt_epoch, crt_obs_new, _matrix_remove);
		if (!updata_valid || crt_obs_new.empty())
		{
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("update_parameter failed"));
			return false;
		}
		chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
		_remove_par_msec += chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count();

		beg_t = chrono::high_resolution_clock::now();
		// ========================================================================================================================================
		vector<string> vec_sites;
		for (const auto& item : map_site_obs) {
			vec_sites.push_back(item.first);
		}
		for (int site_i = 0; site_i < vec_sites.size(); site_i++)
		{
			string crt_rec = vec_sites[site_i];
			_crt_equ.clear_allequ();

			// Process this site and get the equations
			bool proc_site = _processOneRec(crt_epoch, crt_rec, map_site_obs[crt_rec]);
			if (!proc_site)
			{
				_glog->logInfo("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful equations : " + crt_rec));
				continue;
			}

			// ========================================================================================================================================
			// add the new equations
			try
			{
				if (_frequency == 2) _crt_equ.set_newamb(FREQ_1, FREQ_2, dynamic_cast<t_gturboedit*>(_slip12.get()));

				_lsq->add_equation(_crt_equ, crt_epoch, _write_equ);
			}
			catch (...)
			{
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("add_equation throw error"));
				return false;
			}
		}

		end_t = chrono::high_resolution_clock::now();
		_cmb_equ_msec += chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count();
		return true;
	}

	bool t_glsqproc::_processOneEpoch_thread_safe(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs)
	{
		if (crt_obs.empty()) {
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("crt obs is empty"));
			return false;
		}

		for (auto iter = crt_obs.begin(); iter != crt_obs.end();) {
			if (iter->sys() == "R" && iter->channel() >= DEF_CHANNEL) {
				if (_glofrq_num.find(iter->sat()) != _glofrq_num.end()) { 
					iter->channel(_glofrq_num.at(iter->sat())); 
				}
				else {
					_glog->logInfo("t_glsqproc", "_processOneEpoch_thread_safe", crt_epoch.str_mjdsod("no useful frequency id for: " + iter->sat()));
					iter = crt_obs.erase(iter);
					continue;
				}
			}
			++iter;
		}

		chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
		// ========================================================================================================================================
		bool select_obs = _select_obs(crt_epoch, crt_obs);
		if (!select_obs) {
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("select_obs failed"));
			return false;
		}

		map<string, vector<t_gsatdata> > map_site_obs;
		vector<t_gsatdata> crt_obs_new;
		for (const auto& crt_rec : _rec_list) {
			vector<t_gsatdata> crt_rec_obs;
			// get all the obs for the rec
			bool select_valid = _select_rec_obs(crt_rec, crt_epoch, crt_obs, crt_rec_obs);
			if (!select_valid) {
				_glog->logInfo("t_gpcelsqIF", "_processOneEpoch", crt_epoch.str_mjdsod("no useful data : " + crt_rec));
				continue;
			}
			if (_crd_est != CONSTRPAR::KIN) _quality->processOneEpoch(crt_epoch, crt_rec, _rec_crds[crt_rec], crt_rec_obs);
			map_site_obs[crt_rec] = crt_rec_obs;
			crt_obs_new.insert(crt_obs_new.end(), crt_rec_obs.begin(), crt_rec_obs.end());
		}

		bool updata_valid = _lsq->update_parameter(crt_epoch, crt_obs_new, _matrix_remove);
		if (!updata_valid || crt_obs_new.empty())
		{
			_glog->logError("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("update_parameter failed"));
			return false;
		}
		chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
		_remove_par_msec += chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count();

		beg_t = chrono::high_resolution_clock::now();
		vector<string> vec_sites;
		for (const auto& item : map_site_obs) {
			vec_sites.push_back(item.first);
		}

		t_gmutex add_mtx;
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
		for (int site_i = 0; site_i < vec_sites.size(); site_i++)
		{
			t_glsqEquationMatrix equ_temp;
			string crt_rec = vec_sites[site_i];

			// Process this site and get the equations
			bool proc_site = _processOneRec_thread_safe(crt_epoch, crt_rec, map_site_obs[crt_rec], equ_temp);
			if (!proc_site)
			{
				_glog->logInfo("t_glsqproc", "_processOneEpoch", crt_epoch.str_mjdsod("no useful equations : " + crt_rec));
				continue;
			}

			if (_frequency == 2) equ_temp.set_newamb(FREQ_1, FREQ_2, dynamic_cast<t_gturboedit*>(_slip12.get()));

			// add the new equations
			add_mtx.lock();
			_lsq->add_equation(equ_temp, crt_epoch, _write_equ);
			add_mtx.unlock();
		}

		end_t = chrono::high_resolution_clock::now();
		_cmb_equ_msec += chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count();
		return true;
	}
	

	bool t_glsqproc::_processOneRec(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs)
	{
		// ========================================================================================================================================
		if (crt_obs.empty() || crt_rec.empty())
		{
			_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("have no useful data : " + crt_rec));
			return false;
		}

		// save as mark
		_crt_rec = crt_rec;

		// ========================================================================================================================================
		// loop for Combine equation
		vector<double> codeOmc;
		bool while_valid = false;
		do
		{
			// process one rec one sat and get equ
			// ========================================================================================================
			_crt_equ.clear_allequ();

			std::vector<t_gsatdata>::iterator it_beg = crt_obs.begin();
			std::vector<t_gsatdata>::iterator it_end = crt_obs.end();

			for (auto it = it_beg; it != it_end; it++)
			{
				_crt_sat = it->sat();
				if (_sat_list.empty() || _sat_list.find(_crt_sat) == _sat_list.end()) continue;

				bool proc_valid = _processOneRecOneSat(crt_epoch, _crt_rec, _crt_sat, *it);
				if (!proc_valid)
				{
					_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("can not process the rec : " + _crt_rec + " and the sat : " + _crt_sat + " "));
					continue;
				}
			}

			codeOmc = _crt_equ.get_codeomc();
			int idx = _lsq->_x_solve.getParam(crt_rec, par_type::CLK, "");
			double recclk = 0.0;
			double sigma = 0.0;
			int k = Check_Range(codeOmc, recclk, sigma);
			
			if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > _recclk_threshold ||
				(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15)))
			{
				recclk += _lsq->_x_solve[idx].value();
				_lsq->_x_solve[idx].value(recclk);
				_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

				codeOmc.clear();
				while_valid = true;
			}
			else
			{
				while_valid = false;
			}

			// =======================================================================================================
		} while (while_valid);

		if (!_init_sat_clk())
		{
			_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("_init_sat_clk failed" + _crt_sat));
			return false;
		}

		return true;
	}


	bool t_glsqproc::_processOneRec_thread_safe(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs,t_glsqEquationMatrix& equ_result)
	{
		if (crt_obs.empty() || crt_rec.empty())
		{
			_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("have no useful data : " + crt_rec));
			return false;
		}

		// loop for Combine equation
		vector<double> codeOmc;
		bool while_valid = false;
		do
		{
			// process one rec one sat and get equ
			// ========================================================================================================
			equ_result.clear_allequ();

			string crt_sat("NONE");
			for (auto it:crt_obs) 
			{
				crt_sat = it.sat();
				if (_sat_list.empty() || _sat_list.find(crt_sat) == _sat_list.end())
				{
					continue;
				}
				t_gtime epoch = crt_epoch;
				if(!_base_model->cmb_equ(epoch, _lsq->_x_solve, it , equ_result))
				{
					_glog->logInfo("t_glsqproc", "_processOneRec", crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}
			}

			codeOmc = equ_result.get_codeomc();

			//reset recclk omc, check range
			// =======================================================================================================
			int idx = _lsq->_x_solve.getParam(crt_rec, par_type::CLK, "");
			double recclk = 0.0;
			double sigma = 0.0;
			int k = Check_Range(codeOmc, recclk, sigma);
			
			if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > _recclk_threshold/*1e-8*/ ||
				(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15)))
			{
				recclk += _lsq->_x_solve[idx].value();
				_lsq->_x_solve[idx].value(recclk);
				_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

				codeOmc.clear();
				while_valid = true;
			}
			else
			{
				while_valid = false;
			}

			// =======================================================================================================
		} while (while_valid);

		return true;				
	
	}

	bool t_glsqproc::_processOneRecOneSat(const t_gtime& crt_epoch, const std::string& rec, const std::string& sat, t_gsatdata& crt_obs)
	{

		t_gtime epoch = crt_epoch;
		bool code_valid = _base_model->cmb_equ(epoch, _lsq->_x_solve, crt_obs, _crt_equ);
		if (!code_valid)
		{
			_glog->logInfo("t_glsqproc", "_processOneRecOneSat", crt_epoch.str_mjdsod("can not combine equ for code in "));
			return false;
		}
		return true;
	}

	bool t_glsqproc::_init_sat_clk_IF(vector<t_gsatdata>& crt_obs, ColumnVector& l, t_glsq* lsq)
	{
		if (!_lsq || crt_obs.empty()) return false;
		for (int i = 0; i < crt_obs.size(); i++)
		{
			string sat = crt_obs[i].sat();
			int idx_satclk = _lsq->_x_solve.getParam("", par_type::CLK_SAT, sat);
			if (idx_satclk >= 0 && _lsq->_x_solve[idx_satclk].value() == 0.0)
			{
				try
				{
					// satclk delay = - code_omc
					double satclk = -(l(2 * i + 1) + 1.0);
					// code omc
					l(2 * i + 1) = l(2 * i + 1) + satclk;
					// phase omc
					l(2 * i + 2) = l(2 * i + 2) + satclk;
					// clk par init
					_lsq->_x_solve[idx_satclk].value(satclk);

				}
				catch (exception e)
				{
					write_log_info(_glog, 1, "ERROR", e.what());
					return false;
				}
			}
		}

		return true;
	}

	bool t_glsqproc::_init_sat_clk()
	{
		set<string> sat_list = _crt_equ.get_satlist(_crt_rec);
		for (auto sat_iter = sat_list.begin(); sat_iter != sat_list.end(); sat_iter++)
		{
			// get satclk par idx
			int idx_satclk = _lsq->_x_solve.getParam("", par_type::CLK_SAT, *sat_iter);
			if (idx_satclk >= 0 && double_eq(_lsq->_x_solve[idx_satclk].value(),0.0))
			{
				try
				{
					// get all phase code observ equations;
					vector<int> equ_list = _crt_equ.find_equ(_crt_rec, *sat_iter);
					double satclk = 0.0;
					// get satclk delay initial value 
					for (int idx : equ_list)
					{
						if (_crt_equ.get_obscombtype(idx).is_code())
						{
							satclk = -(_crt_equ.l[idx] + 1.0);
							break;
						}
					}
					// clk par init
					_lsq->_x_solve[idx_satclk].value(satclk);
					_bias_model->update_obj_clk(*sat_iter, _crt_time, satclk/CLIGHT);

					for (int idx : equ_list)
					{
						_crt_equ.l[idx] += satclk;
					}
				}
				catch (exception e)
				{
					write_log_info(_glog, 1, "ERROR", e.what());
					return false;
				}
			}
		}
		return true;
	}

	bool t_glsqproc::_init_xml_settings(t_gsetbase* set)
	{
		// for all values with time
		_beg_time = dynamic_cast<t_gsetgen*>(set)->beg();
		_end_time = dynamic_cast<t_gsetgen*>(set)->end();
		_obs_intv = dynamic_cast<t_gsetgen*>(set)->sampling();

		// for all values with mode
		_obs_mode = dynamic_cast<t_gsetproc*>(set)->obs_combin();
		_lsq_mode = dynamic_cast<t_gsetproc*>(set)->lsq_mode();
		_write_equ = dynamic_cast<t_gsetproc*>(set)->write_equ();

		_slip_model = dynamic_cast<t_gsetproc*>(_gset)->slip_model();
		_iono_order = dynamic_cast<t_gsetproc*>(_gset)->iono_order();
		_ztd_model = dynamic_cast<t_gsetproc*>(_gset)->ztd_model(&_pwc_intv);
		_crd_est = dynamic_cast<t_gsetproc*>(_gset)->crd_est();
		_sysbias_model = dynamic_cast<t_gsetproc*>(_gset)->sysbias_model(_isb_pwc_intv);
		_ifb_model = dynamic_cast<t_gsetproc*>(_gset)->ifb_model();
		_frequency = dynamic_cast<t_gsetproc*>(set)->frequency();
		_quality = make_shared<t_gqualitycontrol>(_gset, nullptr);

		// for all settings with sat and rec
		_rec_list = dynamic_cast<t_gsetgen*>(set)->rec_all();
		_sat_list = dynamic_cast<t_gsetgnss*>(set)->sat();
		_sys_list = dynamic_cast<t_gsetgen*>(set)->sys();

		_lite_turboedit = dynamic_cast<t_gsetturboedit*>(set)->liteMode();
		if (_slip_model == SLIPMODEL::DEF_DETECT_MODEL) _lite_turboedit = true;

		_num_threads = dynamic_cast<t_gsetproc*>(set)->num_threads();
		_cmb_equ_multi_thread = dynamic_cast<t_gsetproc*>(set)->cmb_equ_multi_thread();
#ifdef USE_OPENMP
		omp_set_num_threads(_num_threads);
#endif
		_matrix_remove = dynamic_cast<t_gsetproc*>(set)->matrix_remove();

		_maxres_norm = dynamic_cast<t_gsetproc*>(set)->max_res_norm();
		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(set)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(gnut::QZS);

		_freq_number[gnut::GPS] = _freq_index[gnut::GPS].size();
		_freq_number[gnut::GAL] = _freq_index[gnut::GAL].size();
		_freq_number[gnut::GLO] = _freq_index[gnut::GLO].size();
		_freq_number[gnut::BDS] = _freq_index[gnut::BDS].size();
		_freq_number[gnut::QZS] = _freq_index[gnut::QZS].size();

		std::set<string> tmp_sat_list;
		for (const auto& sat : _sat_list)
		{
			auto gsys = t_gsys::sat2gsys(sat);
			auto ssys = t_gsys::gsys2str(gsys);
			if (_sys_list.count(ssys) == 0) continue;
			tmp_sat_list.insert(sat);
		}
		_sat_list = tmp_sat_list;

		if (_sat_list.empty())
		{
			for (auto sys_iter = _sys_list.begin(); sys_iter != _sys_list.end(); sys_iter++)
			{
				auto sats_temp = dynamic_cast<t_gsetgnss*>(set)->sat(t_gsys::str2gsys(*sys_iter));
				_sat_list.insert(sats_temp.begin(), sats_temp.end());
			}
		}

		std::set<string> sat_rm = dynamic_cast<t_gsetgen*>(set)->sat_rm();
		for (const auto& sat : sat_rm) { _sat_list.erase(sat); }


		// Get site coordinates from XML/rinex O files
		std::set<string> rec_rm;
		for (const auto& rec : _rec_list)
		{
			_rec_crds[rec] = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(rec);
			_rec_stds[rec] = dynamic_cast<t_gsetrec*>(_gset)->get_std_xyz(rec);
			if (_rec_crds[rec].zero())
			{
				rec_rm.insert(rec);
			}
		}

		return true;
	}

	bool t_glsqproc::_init_lsq_settings(t_gsetbase* set, t_gallproc* data, t_glog* log)
	{
		// =========================================================================================================
		if (_frequency >= 2) _slip12 = shared_ptr<t_gcycleslip>(new t_gturboedit(set, log));

		// =========================================================================================================
		// Init Paramter in Constructor, make sure lsq is new class
		if (_lsq) { delete _lsq; _lsq = nullptr;} _lsq = new t_glsq(_gset);
		
		_bias_model = shared_ptr<t_gbiasmodel>(new t_gprecisebias(data, log, set));
		
		if (_frequency == 2 && _obs_mode == OBSCOMBIN::IONO_FREE)
		{
			shared_ptr<t_gupdatepar> updateIF(new t_gupdateparIF(_slip12, _band_index));
			updateIF->set_amb_update_way(_lite_turboedit);
			updateIF->set_time(_beg_time, _end_time, _obs_intv);
			updateIF->set_sysbias(_sysbias_model, _isb_pwc_intv);
			updateIF->set_ref_clk(dynamic_cast<t_gsetproc*>(set)->ref_clk());
			_lsq->set_update_par(updateIF);

			_base_model = shared_ptr<t_gbasemodel>(new t_gcombIF(set, log, _bias_model, data));
		}
		else
		{
			throw logic_error("Support IF[DOUBL]");
		}


		return true;
	}

	bool t_glsqproc::_extract_resfile(t_gallrecover& recover_data,string resfile)
	{

		// encode recover file
		t_gfile gout; 
		t_resfile grecover_coder(_gset);
		string recover_path = dynamic_cast<t_gsetout*>(_gset)->outputs("recover");
		if (recover_path.empty())
			recover_path = dynamic_cast<t_gsetout*>(_gset)->outputs("rcv");
		if (recover_path.empty())
			recover_path = resfile;
		if (recover_path.empty())
			recover_path = "resfile_temp_$(date)";

		gnut::substitute(recover_path, "$(date)", _beg_time.str_yyyydoy(), false);
		gnut::substitute(recover_path, "$(rec)", _crt_rec, false);
		//set I/O
		gout.path(recover_path); 
		gout.glog(_glog);
		//set coder
		grecover_coder.glog(_glog);
		grecover_coder.path(recover_path);
		grecover_coder.add_data("ID1", &recover_data);

		gout.coder(&grecover_coder);
		// write
		gout.run_write();

		return true;
	}

	bool t_glsqproc::_extract_clkfile(t_gallrecover& recover_data,t_gallprec::clk_type type)
	{
		// get allprec data
		t_gallprec clk_data;

		recover_data.get_clkdata(clk_data,type);

		// encoder data
		t_gfile gout;
		string clk_path = "";
		if (type == t_gallprec::AS) 
		{
			if(clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("satclk");
			if(clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("cas");
		}
		else if (type == t_gallprec::AR) 
		{
		    if (clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("recclk");
			if (clk_path.empty()) clk_path = dynamic_cast<t_gsetout*>(_gset)->outputs("car");
		}
		
		if (clk_path.empty()) 
		{
			clk_path = "clkfiletemp_$(date)";
		}
		gnut::substitute(clk_path, "$(date)",_beg_time.str_yyyydoy(), false);
		gnut::substitute(clk_path, "$(rec)", _crt_rec, false);
		
		gout.path(clk_path);
		gout.glog(_glog);

		t_rinexc clk_coder(_gset);
		clk_coder.add_data("ID1", &clk_data);
		gout.coder(&clk_coder);

		gout.run_write();
		return true;
	}

	bool t_glsqproc::_init_rec_crd_pars(t_glsq* lsq, bool isPPP)
	{
		int ipar_num = _lsq->_x_solve.parNumber();

		//
		set<string> rec_tmp = _gall_obj->objects_set(t_gdata::REC);
		set<string> rec_final;
		set_intersection(rec_tmp.begin(), rec_tmp.end(), _rec_list.begin(), _rec_list.end(), insert_iterator<set<string>>(rec_final, rec_final.begin()));

		_rec_list = rec_final;
		if (_rec_list.empty())
		{
			write_log_info(_glog, 0, "ERROR : The rec list is empty !");
			return false;
		}

		for (auto it = _rec_list.begin(); it != _rec_list.end(); )
		{
			string site = *it;
			shared_ptr<t_gobj> _grec = _gall_obj->obj(site);


			t_gtriple rec_crd;
			t_gtriple rec_std;

			rec_crd = _rec_crds[site];
			rec_std = _rec_stds[site];

			if (rec_crd.zero()) {
				_grec->get_recent_crd(_beg_time, 100, rec_crd, rec_std);
			}

			//HwangShih: change rec_std using XML
			auto crd_est = dynamic_cast<t_gsetproc*>(_gset)->crd_est();
			if (crd_est == CONSTRPAR::EST)
			{
				rec_std = dynamic_cast<t_gsetpar*>(_gset)->sigRecPosTriple("XXXX");
			}

			if (rec_std.zero() || isPPP)
			{
				rec_std = dynamic_cast<t_gsetpar*>(_gset)->sigRecPosTriple(site);
			}

			if (_crd_est == CONSTRPAR::FIX)
			{
				// only need sinex file 
				t_gtriple snx_crd, snx_std;
				if (!_grec->get_recent_crd(_beg_time,0.1,snx_crd,snx_std)) 
				{					
					if (!_siteres.empty())
					{
						if (_siteres.find(str2upper(site)) != _siteres.end())
						{
							snx_crd[0] = _siteres[str2upper(site)]->get_recover_data_value(par_type::CRD_X);
							snx_crd[1] = _siteres[str2upper(site)]->get_recover_data_value(par_type::CRD_Y);
							snx_crd[2] = _siteres[str2upper(site)]->get_recover_data_value(par_type::CRD_Z);
							if (double_eq(snx_crd[0], 0.0) || double_eq(snx_crd[1], 0.0) || double_eq(snx_crd[2], 0.0))
							{
								_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + "'s crd is 0 in resfile, Erase the site!");
								it = _rec_list.erase(it);
							}
							else
							{
								snx_std = t_gtriple(0.0001, 0.0001, 0.0001);
								_grec->crd(snx_crd, snx_std, FIRST_TIME, LAST_TIME, true);
								it++;
							}
						}
						else
						{
							_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + " has no crd in resfile, Erase the site!");
							it = _rec_list.erase(it);
						}
					}
					else
					{
						snx_crd = dynamic_cast<t_gsetrec*>(_gset)->get_crd_xyz(site);
						snx_std = dynamic_cast<t_gsetrec*>(_gset)->get_std_xyz(site);
						if (snx_std[0] > 0.1 || snx_std[1] > 0.1 || snx_std[2] > 0.1)
						{
							_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + "'s std precision is low, Erase the site!");
							it = _rec_list.erase(it);
							continue;
						}
						if (double_eq(snx_crd[0], 0.0) || double_eq(snx_crd[1], 0.0) || double_eq(snx_crd[2], 0.0))
						{
							_glog->comment(t_glog::LOG_LV::LOG_INFO, "glsqproc", "+init_rec_crd_pars", site + "'s crd precision is absent, Erase the site!");
							it = _rec_list.erase(it);
						}
						else
						{
							snx_std = t_gtriple(0.0001, 0.0001, 0.0001);
							_grec->crd(snx_crd, snx_std, FIRST_TIME, LAST_TIME, true);
							it++;
						}
					}
				}
				else 
				{
					_grec->crd(snx_crd, snx_std,FIRST_TIME,LAST_TIME,true);
					it++;
				}
				continue;
			}

			if (rec_crd.zero() && !_rec_crds[site].zero()) rec_crd = _rec_crds[site];
			else if (!rec_crd.zero() && _rec_crds[site].zero()) _rec_crds[site] = rec_crd;

			// init xyz 
			t_gpar par_x(site, par_type::CRD_X, ++ipar_num, "");
			t_gpar par_y(site, par_type::CRD_Y, ++ipar_num, "");
			t_gpar par_z(site, par_type::CRD_Z, ++ipar_num, "");

			// set the parameter tim
			par_x.value(rec_crd[0]);
			par_y.value(rec_crd[1]);
			par_z.value(rec_crd[2]);

			// set the inital value
			t_gtime endT = _end_time;
			if (_crd_est == CONSTRPAR::KIN || _crd_est == CONSTRPAR::SIMU_KIN) 
			{
				endT = _beg_time;
			}
			par_x.setTime(_beg_time, endT);
			par_y.setTime(_beg_time, endT);
			par_z.setTime(_beg_time, endT);

			// set the apriori  [crd read from rinex]
			par_x.apriori(rec_std[0]);
			par_y.apriori(rec_std[1]);
			par_z.apriori(rec_std[2]);

			//  add parameters
			lsq->add_parameter(par_x);
			lsq->add_parameter(par_y);
			lsq->add_parameter(par_z);

			it++;
		}
		return true;
	}

	bool t_glsqproc::_create_log(t_glog* log)
	{
		if (!log)
		{
			this->_glog = new t_glog();
			this->_glog->mask("podlsqALL.log");
			this->_glog->cache_size(99);
			this->_glog->tsys(t_gtime::GPS);
			this->_glog->time_stamp(true);
			this->_glog->verb(1);
			return true;
		}
		return false;
	}

	bool t_glsqproc::_delete_log(t_glog* log)
	{
		if (log)
		{
			delete log;
			log = nullptr;
			return true;
		}
		return false;
	}
}