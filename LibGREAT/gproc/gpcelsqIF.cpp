/**
 * @file         gpcelsqIF.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gproc/gpcelsqIF.h"
#include "gset/gsetpar.h"
#include "gset/gsetout.h"
#include "gset/gsetgen.h"
#include "gutils/ginfolog.h"
#include "gproc/gupdateparIF.h"

#include "gproc/ginverse_Eigen.h"
#include <gcoders/recover.h>
#include <gio/gfile.h>
#include "gutils/gstring.h"
#include <math.h>
#include <thread>
#include <chrono>

namespace great
{

	t_gpcelsqIF::t_gpcelsqIF(t_gsetbase* set, t_gallproc* data, t_glog* log) :
		t_glsqproc(set, data, log)
	{
		_ref_clk = dynamic_cast<t_gsetproc*>(_gset)->ref_clk();
		_quality_control = make_shared< t_gqualitycontrol>(_gset, nullptr);

		if (_lsq->mode() == LSQMODE::EPO)
		{
			string satclk_file = "clk_" + _beg_time.str_yyyydoy() + "_epo";
			_satclkfile = new t_giof(satclk_file);
			_satclkfile->tsys(t_gtime::GPS);
			_satclkfile->append(false);
		}
	}

	t_gpcelsqIF::~t_gpcelsqIF()
	{
		if (_gioout) { delete _gioout; _gioout; }
		if (_clk_log) { delete _clk_log; _clk_log = nullptr; }
		if (_satclkfile) { delete _satclkfile; _satclkfile = nullptr; }
	}

	bool t_gpcelsqIF::InitProc(t_gallproc* data, const t_gtime& beg, const t_gtime& end)
	{
		if (!data ||
			(beg > end))
		{
			write_log_info(_glog, 0, "ERROR", "beg > end or no proc data");
			return false;
		}

		// get the time
		_beg_time = (_beg_time < beg) ? beg : _beg_time;
		_end_time = (_end_time > end) ? end : _end_time;
		_crt_time = beg;
		_crt_mjd = _crt_time.mjd();

		bool init_data = _initLsqProcData(data);
		if (!init_data)
		{
			_glog->comment(0, "ERROR : can not init data");
			return false;
		}

		bool init_pars = _initLsqProcPars(_lsq);
		if (!init_pars)
		{
			_glog->comment(0, "ERROR : can not init pars");
			return false;
		}

		_quality_control->setNav(_gall_nav);

#ifdef USE_OPENMP
		dynamic_cast<t_gprecisebias*>(_bias_model.get())->set_multi_thread(_rec_list);
		omp_set_num_threads(_num_threads);
#endif

		return true;
	}

	bool t_gpcelsqIF::ProcessBatch(t_gallproc* data, const t_gtime& beg, const t_gtime& end)
	{
		if (!InitProc(data, beg, end)) return false;

		chrono::high_resolution_clock::time_point end_epo_time;

		while (_crt_time <= _end_time) {
			_beg_epo_time = chrono::high_resolution_clock::now();
			cout << _crt_time.str_ymdhms("-------------------") << " -------------------" << endl;
			write_log_info(_glog, 0, "NOTE", _crt_time.str_ymdhms("Processing epoch "));

			_initOneEpoch();
			/* get observation */
			vector<t_gsatdata> crt_obs = _gall_obs->obs(_rec_list, _crt_time);
			bool epoch_valid = false;
			int num_site = _gall_obs->end_obs_sites(_crt_time, _rec_list);

			epoch_valid = _processOneEpoch(_crt_time, crt_obs);

			if (!epoch_valid) {
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_gpcelsqIF", "ProcessBatch", "Processing failed in epoch");
			}

			if (_lsq->mode() == LSQMODE::EPO) {
				_solveEpoch(crt_obs);
			}

			write_log_info(_glog, 1, "NOTE", _crt_time.str_ymdhms("End Processing epoch "));
			end_epo_time = chrono::high_resolution_clock::now();
			double compute_time = chrono::duration_cast<chrono::milliseconds>(end_epo_time - _beg_epo_time).count() / 1000.0;
			if (_lsq->mode() == LSQMODE::EPO) {
				cout << _crt_time.str_ymdhms("Finish epoch") << ": " << setw(8) << setprecision(3) << fixed << _prepare_time << setw(8) << setprecision(3) << fixed << compute_time
					<< setw(8) << setprecision(3) << fixed << _prepare_time + compute_time << " sec" << ", nrec = " << setw(3) << _map_all_equ.size()
					<< ", nobs = " << setw(5) << _obs_crt_num << ", npar = " << setw(5) << _lsq->_x_solve.parNumber() << ", sigma0 = " << setw(9) << setprecision(5) << _crt_sigma << endl;
			}
			else {
				cout << _crt_time.str_ymdhms("Finish epoch") << ": " << setw(8) << setprecision(3) << fixed << _prepare_time << setw(8) << setprecision(3) << fixed << compute_time
					<< setw(8) << setprecision(3) << fixed << _prepare_time + compute_time << " sec" << ", npar = " << setw(5) << _lsq->_x_solve.parNumber() << endl;
			}
			
			_crt_time = _crt_time + _obs_intv;
		}

		writeLogInfo(_glog, 0, "NOTE", "###REMOVE_PAR " + dbl2str(_remove_par_msec / 1000.0) + " sec.");
		writeLogInfo(_glog, 0, "NOTE", "###COMBINE_EQU " + dbl2str(_cmb_equ_msec / 1000.0) + " sec.");

		chrono::high_resolution_clock::time_point beg_t = chrono::high_resolution_clock::now();
		try
		{
			_check_ref_clk(_lsq);
			_lsq->solve_NEQ();
		}
		catch (exception e)
		{
			write_log_info(_glog, 0, e.what(), "Solve Equation Fail!");
			return false;
		}
		chrono::high_resolution_clock::time_point end_t = chrono::high_resolution_clock::now();
		writeLogInfo(_glog, 0, "NOTE", "###SOLVE_LSQ " + dbl2str(chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count() / 1000.0) + " sec.");
		
		beg_t = chrono::high_resolution_clock::now();
		try
		{
			_lsq->get_result_parameter(*_gall_rcv);
			_lsq->recover_parameters(*_gall_rcv);
			_extract_resfile(*_gall_rcv);
		}
		catch (exception e)
		{
			write_log_info(_glog, 1, e.what(), "Recover Parameter Fail!");
			return false;
		}
		end_t = chrono::high_resolution_clock::now();
		writeLogInfo(_glog, 0, "NOTE", "###RECOVER " + dbl2str(chrono::duration_cast<chrono::milliseconds>(end_t - beg_t).count() / 1000.0) + " sec.");

		return true;
	}

	bool t_gpcelsqIF::GenerateProduct()
	{
		// update clk file
		_extract_clkfile(*_gall_rcv, t_gallprec::AS);
		_extract_clkfile(*_gall_rcv, t_gallprec::AR);
		return true;
	}

	bool t_gpcelsqIF::_initLsqProdData(t_gallprod* data)
	{
		return false;
	}

	bool t_gpcelsqIF::_initLsqProcPars(t_glsq* lsq)
	{
		string class_id = "t_gpcelsqIF";
		string funct_id = "_initLsqProcPars";

		if (_rec_list.empty() || _sat_list.empty())
		{
			this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "no rec or sat");
			return false;
		}

		// check reference clock
		if (!_ref_clk.empty() && _sat_list.find(_ref_clk) == _sat_list.end() && _rec_list.find(_ref_clk) == _rec_list.end()) {
			_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "reference clock not in satellite list and site list"); 
			return false; 
		}
		
		// Crd pars
		bool rec_crd_valid = _init_rec_crd_pars(lsq);
		bool sat_clk_valid = _init_sat_clk_pars(lsq);
		bool rec_clk_valid = _init_rec_clk_pars(lsq);
		bool trop_valid = _init_rec_trop_pars(lsq);

		if (!rec_crd_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "rec_crd_valid  is false"); return false; }
		if (!sat_clk_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "sat_clk_valid  is false"); return false; }
		if (!rec_clk_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "rec_clk_valid  is false"); return false; }
		if (!trop_valid) { this->_glog->comment(t_glog::LOG_LV::LOG_ERROR, class_id, funct_id, "trop_valid  is false"); return false; }

		write_log_info(_glog, 1, "INIT", "init the pars estimated in PCE_LSQ_IF.");
		return true;
	}

	bool t_gpcelsqIF::_initLsqProcData(t_gallproc* data)
	{
		gtrace("t_gpcelsqIF::_initLsqProcData");

		// prepare obs data
		_gall_bias = dynamic_cast<t_gallbias*>((*data)[t_gdata::ALLBIAS]);
		_gall_obs = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS]);
		_gall_obj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gall_rcv = dynamic_cast<t_gallrecover*>((*data)[t_gdata::ALLRECOVER]);
		_gall_nav = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gall_nav->overwrite(true);

		if (!_gall_bias ||
			!_gall_obs ||
			!_gall_obj ||
			!_gall_rcv)
		{
			write_log_info(_glog, 0, "ERROR : can not init lsq process data");
			return false;
		}

		// set glofrq number 
		_glofrq_num = _gall_obs->glo_freq_num();
		if (_glofrq_num.empty() && _gall_nav)
		{
			_glofrq_num = _gall_nav->glo_freq_num();
		}

		write_log_info(_glog, 4, "INIT : process data in lsq is prepared.");
		return true;
	}

	bool t_gpcelsqIF::_processOneRecOneSat(const t_gtime& crt_epoch,
		const std::string& rec, const std::string& sat, t_gsatdata& crt_obs)
	{
		t_gtime epoch = crt_epoch;
		bool proc_valid = _base_model->cmb_equ(epoch, _lsq->_x_solve, crt_obs, _crt_equ);
		string class_id = "t_gpcelsqIF";
		string funct_id = "_processOneRecOneSat";

		if (!proc_valid)
		{
			_glog->comment(t_glog::LOG_LV::LOG_WARN, class_id, funct_id,
				crt_epoch.str_mjdsod("can not combine equ for code in " + sat + rec + "at "));
			return false;
		}
		return true;
	}

	bool t_gpcelsqIF::_processOneEpoch(const t_gtime& crt_epoch, std::vector<t_gsatdata>& crt_obs)
	{
		string class_id = "t_gpcelsqIF";
		string funct_id = "_processOneEpoch";
		_map_all_equ.clear();

		for (auto iter = crt_obs.begin(); iter != crt_obs.end();) {
			if (iter->sys() == "R" && iter->channel() >= DEF_CHANNEL) {
				if (_glofrq_num.find(iter->sat()) != _glofrq_num.end()) {
					iter->channel(_glofrq_num.at(iter->sat()));
				}
				else {
					_glog->logInfo("t_gpcelsqIF", "_processOneEpoch", crt_epoch.str_mjdsod("no useful frequency id for: " + iter->sat()));
					iter = crt_obs.erase(iter);
					continue;
				}
			}
			++iter;
		}

		bool select_obs = _select_obs(crt_epoch, crt_obs);
		if (!select_obs || !_slip12 || crt_obs.empty()) {
			_glog->logError("t_gpcelsqIF", "_processOneEpoch", crt_epoch.str_mjdsod("select_obs failed"));
			return false;
		}

		// prepare Obs
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
			if (_crd_est != CONSTRPAR::KIN) _quality_control->processOneEpoch(crt_epoch, crt_rec, _rec_crds[crt_rec], crt_rec_obs);

			//number of sat less than 4 continue
			if (crt_rec_obs.size() < 4) {
				continue;
			}
			map_site_obs[crt_rec] = crt_rec_obs;
			crt_obs_new.insert(crt_obs_new.end(), crt_rec_obs.begin(), crt_rec_obs.end());
		}

		chrono::high_resolution_clock::time_point beg_time = chrono::high_resolution_clock::now();
		bool updata_valid = _lsq->update_parameter(crt_epoch, crt_obs_new, _matrix_remove, _write_equ);
		chrono::high_resolution_clock::time_point end_time = chrono::high_resolution_clock::now();
		_remove_par_msec += chrono::duration_cast<chrono::milliseconds>(end_time - beg_time).count();

		if (!updata_valid || crt_obs_new.empty()) {
			write_log_info(_glog, 1, crt_epoch.str(), ": no observation found ");
			return false;
		}

		vector<string> vec_sites;
		for (const auto& item : map_site_obs) {
			vec_sites.push_back(item.first);
		}
		if (vec_sites.empty()) return false;

		beg_time = chrono::high_resolution_clock::now();
		t_gmutex add_mtx;
		// multi-thread  Process REC [Use OpenMP]
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
		for (int site_i = 0; site_i < vec_sites.size(); site_i++) {
			t_glsqEquationMatrix equ_temp;
			string rec_temp = vec_sites[site_i];
			// Process this site and get the equations
			bool proc_site = _processOneRec_thread_safe(crt_epoch, rec_temp, map_site_obs[rec_temp], equ_temp);
			if (!proc_site) {
				cout << crt_epoch.str_ymdhms(rec_temp + " has no equations ", false, false) << endl;
				_glog->logInfo("t_gpcelsqIF", "_processOneEpoch", crt_epoch.str_mjdsod("no useful equations : " + rec_temp));
				continue;
			}

			// add the new equations
			add_mtx.lock();
			if (_lsq->mode() == LSQMODE::EPO) {
				if (equ_temp.num_equ() > 0) _map_all_equ[rec_temp] = equ_temp;
			}
			else {
				_lsq->add_equation(equ_temp, crt_epoch, _write_equ);
				_obs_crt_num = _lsq->get_equ_obs_total_num();
			}
			add_mtx.unlock();
		}
		end_time = chrono::high_resolution_clock::now();
		_cmb_equ_msec += chrono::duration_cast<chrono::milliseconds>(end_time - beg_time).count();
		_glog->logDebug("t_gpcelsqIF", "_processOneEpoch", "Finish form equations");
		return true;
	}

	bool t_gpcelsqIF::_processOneRec_thread_safe(const t_gtime& crt_epoch, const std::string& crt_rec, std::vector<t_gsatdata>& crt_obs, t_glsqEquationMatrix& equ_result)
	{
		const string class_id = "t_gpcelsqIF";
		const string funct_id = "_processOneRec";

		if (crt_obs.empty() || crt_rec.empty()) {
			if (_glog) _glog->comment(2, class_id + ": " + funct_id, "This epoch have no useful data: " + crt_rec);
			return false;
		}

		vector<double> codeOmc;

		// loop for Combine equation
		bool while_valid = false;
		do {
			// process one rec one sat and get equ
			equ_result.clear_allequ();
			string crt_sat;
			for (auto& it : crt_obs) {
				crt_sat = it.sat();
				if (_sat_list.find(crt_sat) == _sat_list.end()) continue;
				t_gtime epo(crt_epoch);
				if (!_base_model->cmb_equ(epo, _lsq->_x_solve, it, equ_result)) {
					_glog->comment(t_glog::LOG_LV::LOG_WARN, class_id, funct_id, crt_epoch.str_mjdsod("can not process the rec : " + crt_rec + " and the sat : " + crt_sat + " "));
					continue;
				}
			}

			codeOmc = equ_result.get_codeomc();
			//reset recclk omc, check range
			double recclk;
			double sigma;

			int k = Check_Range(codeOmc, recclk, sigma);
			int idx = _lsq->_x_solve.getParam(crt_rec, par_type::CLK, "");
			if (idx >= 0 && !while_valid && (fabs(recclk / CLIGHT) > 1e-6 ||
				(k >= 1 && fabs(_lsq->_x_solve[idx].value() / CLIGHT) < 1e-15))) {
				// reset recclk
				recclk += _lsq->_x_solve[idx].value();
				_lsq->_x_solve[idx].value(recclk);
				_bias_model->update_obj_clk(crt_rec, crt_epoch, recclk / CLIGHT);

				codeOmc.clear();
				while_valid = true;
			}
			else {
				while_valid = false;
			}
		} while (while_valid);
		return true;
	}

	bool t_gpcelsqIF::_solveEpoch(vector<t_gsatdata>& crt_obs)
	{
		_glog->logDebug("t_gpcelsqIF", "_solveEpoch", "Begin solve epoch clock");
		for (const auto& equ : _map_all_equ) {
			_lsq->add_equation(equ.second, _crt_time, _write_equ);
		}

		_get_obs_crt_num();
		_check_ref_clk(_lsq);

		if (!_epo_solved && _obs_crt_num > _sat_list.size() + _rec_list.size()) {
			t_glsq lsqepo(*_lsq);
			try {
				lsqepo.solve_x();
				_get_clk_crt(&lsqepo, false);
				_epo_solved = true;
				_crt_sigma = lsqepo.sigma0();
			}
			catch (...) {
				_glog->comment(t_glog::LOG_LV::LOG_ERROR, "t_gpcelsqIF", "_solveEpoch", _crt_time.str_mjdsod("solve_x throw error"));
				for (const auto& sat : _sat_list) _clk_crt[sat] = 0;
				return false;
			}
		}

		double crt_mjd = _crt_time.mjd();
		if (crt_mjd > _crt_mjd)
		{
			//delete last day
			delete _satclkfile;
			_satclkfile = nullptr;
			//creat current day
			string satclk_file = "clk_" + _crt_time.str_yyyydoy() + "_epo";
			_satclkfile = new t_giof(satclk_file);
			_satclkfile->tsys(t_gtime::GPS);
			_satclkfile->append(false);

			_crt_mjd = crt_mjd;
		}
		if (_satclkfile) _print_clk_solved();

		_glog->logDebug("t_gpcelsqIF", "_solveEpoch", "Finish solve epoch clock");
		return true;
	}

	bool t_gpcelsqIF::_ref_clk_valid(t_glsq* lsq)
	{
		bool valid = false;
		int idx = -1;
		if (_ref_clk.length() == 3) {
			idx = lsq->_x_solve.getParam("", par_type::CLK_SAT, _ref_clk);
		}
		else if (_ref_clk.length() == 4) {
			idx = lsq->_x_solve.getParam(_ref_clk, par_type::CLK, "");
		}

		if (idx >= 0) valid = lsq->par_alive(idx);
		if (!valid) _glog->logInfo("t_gpcelsqIF", "_ref_clk_valid", "reference clock not valid: " + _ref_clk);
		return valid;
	}

	bool t_gpcelsqIF::_check_ref_clk(t_glsq* lsq)
	{
		if (_ref_clk.empty()) { lsq->lsq_clk_constraint(); return true; }
		if (_ref_clk_valid(lsq)) return true;

		// if use receiver clock as reference, do not change to a satellite
		if (_ref_clk.length() != 3) return false;
		// change the reference clock
		const string& sat_ref1 = _sat_obs_max(lsq);
		if (!sat_ref1.empty()) {
			int idx0 = -1;
			if (_ref_clk.length() == 3) idx0 = lsq->_x_solve.getParam("", par_type::CLK_SAT, _ref_clk);
			else idx0 = lsq->_x_solve.getParam(_ref_clk, par_type::CLK, "");
			int idx1 = lsq->_x_solve.getParam("", par_type::CLK_SAT, sat_ref1);
			if (idx0 > 0 && idx1 > 0) {
				double sig0 = lsq->_x_solve.getParSig(idx0);
				double sig1 = lsq->_x_solve.getParSig(idx1);
				lsq->_x_solve.setParSig(idx0, sig1);
				lsq->_x_solve.setParSig(idx1, sig0);
				cout << "reference clock not valid: " << _ref_clk << ", change to " << sat_ref1 << endl;
				_glog->logInfo("t_gpcelsqIF", "_check_ref_clk", "reference clock not valid: " + _ref_clk + ", change to " + sat_ref1);
				_ref_clk = sat_ref1;
				_ref_clk_crt = sat_ref1;
				return true;
			}
		}
		return false;
	}

	void t_gpcelsqIF::_get_clk_crt(t_glsq* lsq, bool update_std)
	{
		const ColumnVector& dx = lsq->dx();
		for (const auto& sat : _sat_list) {
			int idx = lsq->_x_solve.getParam("", par_type::CLK_SAT, sat, _crt_time, _crt_time);
			if (idx < 0) continue;
			double dclk = dx(idx + 1);
			if (dclk != 0) {
				_clk_crt[sat] = lsq->_x_solve.getParValue(idx) + dclk;
				if (update_std) _clk_std_crt[sat] = lsq->stdx(idx + 1);
				else _clk_std_crt[sat] = 0.05;
			}
		}
	}

	void t_gpcelsqIF::_get_obs_crt_num()
	{
		int value = 0;
		for (auto& iter : _map_all_equ) {
			auto& equ = iter.second;
			value += equ.num_equ();
		}
		_obs_crt_num = value;
	}

	string t_gpcelsqIF::_sat_obs_max(t_glsq* lsq)
	{
		map<string, size_t> sat_count;
		set<GSYS> gsys;
		for (auto& iter : _map_all_equ) {
			auto& equ = iter.second;
			for (int i = 0; i < equ.num_equ(); ++i) {
				const string& sat = iter.second.get_satname(i);
				sat_count[sat]++;
				gsys.insert(t_gsys::sat2gsys(sat));
			}
		}
		// choose the used GSYS
		// site as reference: choose a G/E/C/R satellite instead
		// satellite as reference: choose a same system satellite instead
		GSYS gs = GSYS::GPS;
		vector<GSYS> gs_order = { GSYS::GPS, GSYS::GAL, GSYS::BDS, GSYS::GLO, GSYS::QZS };
		if (_ref_clk.length() == 3 && gsys.find(t_gsys::sat2gsys(_ref_clk)) != gsys.end()) {
			gs = t_gsys::sat2gsys(_ref_clk);
		}
		else {
			for (const auto& iter : gs_order) {
				if (gsys.find(iter) != gsys.end()) { gs = iter; break; }
			}
		}
		size_t max = 0;
		string sat_max;
		for (const auto& iter : sat_count) {
			if (t_gsys::sat2gsys(iter.first) != gs) continue;
			int idx = lsq->_x_solve.getParam("", par_type::CLK_SAT, iter.first, _crt_time, _crt_time);
			if (idx < 0) continue;
			if (iter.second > max) {
				sat_max = iter.first;
				max = iter.second;
			}
		}
		return sat_max;
	}

	bool t_gpcelsqIF::_initOneEpoch()
	{
		if (_lsq->mode() == LSQMODE::EPO) {
			_gall_nav->clean_invalid();
			_gall_nav->clean_duplicit();

			// update log file every 6h
			if (_crt_time.sod() % (6 * 3600) == 0) {
				int jpos = _crt_time.sod() / (6 * 3600);
				stringstream ss;
				ss << "great_pcelsq.app_log_" << _crt_time.str_yyyydoy() << "_" << setw(2) << setfill('0') << jpos * 6 << setfill(' ');
				const string& logfile = ss.str();
				int verb = _glog->verb();
				if (_clk_log) {
					delete _clk_log;
					_clk_log = nullptr;
				}
				_clk_log = new t_glog(logfile);
				_clk_log->cache_size(99);
				_clk_log->tsys(t_gtime::GPS);
				_clk_log->time_stamp(true);
				_clk_log->verb(verb);
				_glog = _clk_log;
			}
		}
		
		_map_all_equ.erase(_map_all_equ.begin(), _map_all_equ.end());
		_obs_crt_num = 0;

		_ref_clk_crt = _ref_clk;
		_epo_solved = false;
		for (const auto& sat : _sat_list) {
			_clk_crt[sat] = 0;
			_clk_std_crt[sat] = 0;
		}

		auto end_epo_time = chrono::high_resolution_clock::now();
		_prepare_time = chrono::duration_cast<chrono::milliseconds>(end_epo_time - _beg_epo_time).count() / 1000.0;
		_beg_epo_time = chrono::high_resolution_clock::now();

		return true;
	}

	void t_gpcelsqIF::_print_clk_solved()
	{
		ostringstream os;
		t_gtime creat_t(t_gtime::UTC);
		creat_t.add_secs(28800);

		// write sp3 clk file header
		if (!_satclkfile->is_open()) {
			os << "     2.00           C                                       RINEX VERSION / TYPE" << endl;
			os << left << setw(20) << "GREAT" << setw(20) << "GREAT" << setw(20) << creat_t.str("%Y%m%d %H%M%S") << "PGM / RUN BY / DATE" << endl;
			os << right << setw(9) << fixed << setprecision(2) << _obs_intv << setw(51) << " " << "INTERVAL" << endl;
			// write # / TYPES OF DATA
			os << "    " << right << setw(2) << "1" << "    " << "AS" << setw(48) << "" << "# / TYPES OF DATA" << endl;
			// write # OF CLK REF
			os << "    " << right << setw(2) << "1" << setw(54) << "" << "# OF CLK REF" << endl;
			// write ANALYSIS CLK REF
			os << setw(4) << _ref_clk << setw(56) << " " << "ANALYSIS CLK REF" << endl;
			// write END OF HEADER
			os << setw(60) << " " << "END OF HEADER" << endl;
		}

		// Print sp3 clk file header
		if (_satclkfile) {
			_satclkfile->write(os.str().c_str(), os.str().size());
			_satclkfile->flush();
			os.str("");
		}
		else {
			_glog->comment(2, "t_gpcelsqIF::_write_clk_header", "have no output file!");
		}

		// write sp3 clk file data
		for (const auto& sat : _sat_list) {
			if (_clk_crt.find(sat) == _clk_crt.end() || double_eq(_clk_crt[sat], 0)) continue;
			os << "AS" << " " << setw(3) << sat << "  " <<
				setw(26) << _crt_time.str("%Y %C %L %K %O %P") << " " << setw(2)
				<< "1" << "   " << setw(19) << fixed << scientific << setprecision(12) << _clk_crt[sat] / CLIGHT << endl;
		}

		// Print sp3 clk file data
		if (_satclkfile) {
			_satclkfile->write(os.str().c_str(), os.str().size());
			_satclkfile->flush();
			os.str("");
		}
		else {
			_glog->comment(2, "t_gpcelsqIF::_printout", "have no output file!");
		}
	}
}