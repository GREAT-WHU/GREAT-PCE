/**
 * @file         gupdatepar.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gupdatepar.h"

namespace great
{

	t_gupdateparinfo::t_gupdateparinfo()
	{
	}
	t_gupdateparinfo::~t_gupdateparinfo()
	{
	}

	bool t_gupdateparinfo::exist(int id)
	{
		return _remove_id.count(id)!=0 ? true : false;
	}

	void t_gupdateparinfo::add(int id)
	{
		if (_remove_id.count(id) != 0) 
		{
			throw std::logic_error("Remove id have repeated!");
		}

		_remove_id.insert(id);

	}

	void t_gupdateparinfo::add(t_gpar newpar)
	{
		_new_parlist.push_back(newpar);
	}

	void t_gupdateparinfo::add(int id, t_gpar newpar)
	{
		this->add(id);
		_new_parlist.push_back(newpar);
	}

	void t_gupdateparinfo::add(vector<pair<int, t_gpar> > update_par, t_glsqEquationMatrix state_equ)
	{
		for (auto iter : update_par) {
			//this->add(iter.first, iter.second);
			this->add(iter.first);
			_equ_parlist.push_back(iter.second);
		}
		_state_equ.add_equ(state_equ);
	}

	void t_gupdateparinfo::get(vector<int>& remove_id)
	{
		remove_id = vector<int>(_remove_id.begin(),_remove_id.end());
	}

	void t_gupdateparinfo::get(vector<t_gpar>& newparlist)
	{
		newparlist = _new_parlist;
	}

	void t_gupdateparinfo::get(vector<t_gpar>& update_newpar, t_glsqEquationMatrix& equ)
	{
		update_newpar = _equ_parlist;
		equ = _state_equ;
	}

	void t_gupdateparinfo::get(vector<int>& remove_id, vector<t_gpar>& newparlist, vector<t_gpar>& equ_parlist, t_glsqEquationMatrix & equ)
	{
		get(remove_id);
		get(newparlist);
		get(equ_parlist, equ);
	}



	t_gupdatepar::t_gupdatepar()
	{
	}

	t_gupdatepar::~t_gupdatepar()
	{
	}

	t_gupdatepar::t_gupdatepar(const t_gupdatepar & Other):
		_state_mode(Other._state_mode),
		_cycleslip(Other._cycleslip)
	{
		
	}

	void t_gupdatepar::set_interval(double interval)
	{
		_intv = interval;
	}

	void t_gupdatepar::set_sig_ion(double sigion)
	{
		_sig_ion = sigion;
	}

	void t_gupdatepar::set_cycleslip(shared_ptr<t_gcycleslip> cycleslip)
	{
		_cycleslip = cycleslip;
	}
	
	void t_gupdatepar::set_amb_update_way(bool way)
	{
	    _lite_update_amb = way;
	}

	void t_gupdatepar::set_par_state_mode(par_type type, int order, double dt, double noise)
	{
		_state_mode[type] = t_statemode(order, dt, noise);
	}

	void t_gupdatepar::update_amb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)
	{
		this->_update_amb_pars(epoch, allpars, obsdata, update_info);
	}

	void t_gupdatepar::update_isb_pars(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)
	{
		map<string, set<string>> site_sys;
		map<string, set<string>> sys_site;
		map<pair<string, string>, int> nobs_site_sys;
		bool est_ifb = _sys_bias == SYSBIASMODEL::AUTO_RWK || _sys_bias == SYSBIASMODEL::AUTO_WHIT ||
			_sys_bias == SYSBIASMODEL::AUTO_PWC || _sys_bias == SYSBIASMODEL::AUTO_CON;
		map<string, int> glo_fid;
		set<int> fid_glo;

		for (const auto& obs : obsdata) {
			const auto& sat = obs.sat();	
			// consider ISB between BDS2 and BDS3
			if (_use_bds2_isb && sat.substr(0, 1) == "C" && sat < "C17") {
				site_sys[obs.site()].insert("C2");
				sys_site["C2"].insert(obs.site());
				nobs_site_sys[make_pair<string, string>(obs.site(), "C2")]++;
			}
			else if (sat.substr(0, 1) == "R") {
				string gsys = "R";
				if (est_ifb)  gsys += int2str(obs.channel(), 2);
				site_sys[obs.site()].insert(gsys);
				sys_site[gsys].insert(obs.site());
				nobs_site_sys[{obs.site(), gsys}]++;
				if (glo_fid.find(obs.sat()) == glo_fid.end()) glo_fid[obs.sat()] = obs.channel();
				fid_glo.insert(obs.channel());
			}
			else {
				site_sys[obs.site()].insert(sat.substr(0, 1));
				sys_site[sat.substr(0, 1)].insert(obs.site());
				nobs_site_sys[make_pair<string, string>(obs.site(), sat.substr(0, 1))]++;
			}
		}
		vector<string> ref_order = { "G", "E", "C", "C2", "J" };
		set<par_type> isb_list = { par_type::BDS_ISB, par_type::BD2_ISB, par_type::GAL_ISB, par_type::GLO_ISB, par_type::QZS_ISB };
		// find CLK_SAT parameter to determine wether need a reference ISB for one system
		set<string> sys_est_clk;
		for (const auto& par : allpars.getAllPar()) {
			if (par.parType == par_type::CLK_SAT) {
				const string& sat = par.prn;
				if (_use_bds2_isb && sat.substr(0, 1) == "C" && sat < "C17")
					sys_est_clk.insert("C2");
				else if (sat.substr(0, 1) == "R" && est_ifb && glo_fid.find(sat) != glo_fid.end()) 
					sys_est_clk.insert("R" + int2str(glo_fid[sat], 2));
				else 
					sys_est_clk.insert(sat.substr(0, 1));
			}
		}
		set<string> new_ref_sys;
		for (const string& gsys : sys_est_clk) {
			set<string> sites = sys_site[gsys];
			if (sites.empty()) continue;
			if (sites.find(_sys_ref_site[gsys]) != sites.end()) continue;
			if (sites.find(_ref_clk) != sites.end()) {
				_sys_ref_site[gsys] = _ref_clk;
				new_ref_sys.insert(gsys);
			}
			else {
				string ref_site;
				int max_sys = 1;
				int max_nobs = 0;
				for (const auto& site : sites) {
					if (site_sys[site].size() >= max_sys && nobs_site_sys[{site, gsys}] > max_nobs) {
						ref_site = site;
						max_sys = site_sys[site].size();
						max_nobs = nobs_site_sys[{site, gsys}];
					}
				}
				_sys_ref_site[gsys] = ref_site;
				new_ref_sys.insert(gsys);
			}
		}

		const int idx_max = allpars.maxIndex();
		// int num = 0;
		for (const auto& iter : site_sys) {
			const auto& site = iter.first;
			const auto& sys_list = iter.second;
			map<string, int> isb_par_idx;
			for (const auto& par : isb_list) {
				if (est_ifb && par == par_type::GLO_ISB) continue;
				int idx_isb = allpars.getParam(site, par, "");
				if (idx_isb >= 0) isb_par_idx[allpars[idx_isb].str_type()] = idx_isb;
			}
			if (est_ifb) {
				for (const auto& fid : fid_glo) {
					int idx_ifb = allpars.getParam(site, par_type::GLO_IFB, "", FIRST_TIME, LAST_TIME, fid);
					if (idx_ifb >= 0) isb_par_idx[allpars[idx_ifb].str_type()] = idx_ifb;
				}
			}
			// system num < 2, remove all ISB pars
			if (sys_list.size() < 2) {
				for (const auto& par : isb_par_idx) {
					if (allpars[par.second].end > epoch - _intv) allpars[par.second].end = epoch - _intv;
					update_info.add(par.second + 1);
				}
				continue;
			}
			string ref_sys;
			for (const auto& sys : ref_order) {
				if (sys_list.find(sys) != sys_list.end()) { ref_sys = sys;	break; }
			}
			// GLONASS only, choose one channel as reference
			if (ref_sys.empty() && est_ifb) {
				for (const auto& sys : sys_list) {
					if (sys.substr(0, 1) == "R") { ref_sys = sys; break; }
				}
			}
			// if the referece system changes, all ISB should be updated
			bool ref_change = false;
			if (_site_ref_sys.find(site) == _site_ref_sys.end()) {
				_site_ref_sys[site] = ref_sys;
			}
			else if (_site_ref_sys[site] != ref_sys) {
				for (const auto& par : isb_par_idx) {
					if (allpars[par.second].end > epoch - _intv) allpars[par.second].end = epoch - _intv;
					update_info.add(par.second + 1);
				}
				_site_ref_sys[site] = ref_sys;
				ref_change = true;
			}
			for (const auto& sys : sys_list) {
				if (sys == ref_sys) continue;
				string isb_type;
				par_type isb_par = par_type::NO_DEF;
				int channel = DEF_CHANNEL;
				if (sys == "C") { 
					isb_type = "BDS_ISB"; isb_par = par_type::BDS_ISB; 
				}
				else if (sys == "C2") {
					isb_type = "BD2_ISB"; isb_par = par_type::BD2_ISB;
				}
				else if (sys == "E") {
					isb_type = "GAL_ISB"; isb_par = par_type::GAL_ISB;
				}
				else if (sys == "J") {
					isb_type = "QZS_ISB"; isb_par = par_type::QZS_ISB;
				}
				else if (est_ifb && sys.substr(0, 1) == "R") {
					channel = str2int(sys.substr(1));
					isb_type = "GLO_IFB_" + sys.substr(1);
					isb_par = par_type::GLO_IFB;
				}
				else if (sys == "R") {
					isb_type = "GLO_ISB"; isb_par = par_type::GLO_ISB;
				}
				if (isb_type.empty() || isb_par == par_type::NO_DEF) continue;
				// create new ISB parameters: 
				// 1. the reference site of one sys is changed, all ISB of this sys have to be reset
				// 2. the reference sys of one site is changed, all ISB of this site have to be reset
				// 3. no ISB parameter for corresponding observations, create the ISB for this site and sys
				if (new_ref_sys.find(sys) != new_ref_sys.end() || ref_change || isb_par_idx.find(isb_type) == isb_par_idx.end()) {
					if (isb_par_idx.find(isb_type) != isb_par_idx.end() && !ref_change) {
						const int idx_isb = isb_par_idx[isb_type];
						if (allpars[idx_isb].end > epoch - _intv) allpars[idx_isb].end = epoch - _intv;
						update_info.add(idx_isb + 1);
					}
					t_gpar par_ISB;
					par_ISB = t_gpar(site, isb_par, idx_max + 1 + _new_par, "");
					par_ISB.channel = channel;
					par_ISB.value(0.0);
					if (_sys_ref_site[sys] == site) par_ISB.apriori(0.001);
					else par_ISB.apriori(3000.0);
					if (_sys_bias == SYSBIASMODEL::AUTO_CON || _sys_bias == SYSBIASMODEL::ISB_CON) {
						par_ISB.setTime(epoch, _end_time);
					}
					else if (_sys_bias == SYSBIASMODEL::AUTO_RWK || _sys_bias == SYSBIASMODEL::ISB_RWK) {
						par_ISB.setTime(epoch, epoch);
						set_par_state_mode(isb_par, 1, 24, 3000.0);
					}
					else if (_sys_bias == SYSBIASMODEL::AUTO_PWC || _sys_bias == SYSBIASMODEL::ISB_PWC) {
						t_gtime pwc_time = epoch + 3600 * _isb_intv - _intv;
						if (pwc_time > _end_time) pwc_time = _end_time;
						par_ISB.setTime(epoch, pwc_time);
						set_par_state_mode(isb_par, 1, _isb_intv, 3000.0);
					}
					else {
						par_ISB.setTime(epoch, epoch);
					}
					_new_par++;
					update_info.add(par_ISB);
				}
				// update existing ISB par
				else {	
					if (_sys_bias == SYSBIASMODEL::AUTO_CON || _sys_bias == SYSBIASMODEL::ISB_CON) continue;
					const int idx_isb = isb_par_idx[isb_type];
					if (allpars[idx_isb].end > epoch - _intv) {
						if (_sys_bias == SYSBIASMODEL::AUTO_PWC || _sys_bias == SYSBIASMODEL::ISB_PWC) continue;
						allpars[idx_isb].end = epoch - _intv;
					}
					if (_sys_bias == SYSBIASMODEL::AUTO_WHIT || _sys_bias == SYSBIASMODEL::ISB_WHIT)
						_update_process_par(epoch, idx_isb, allpars, update_info);
					else
						_update_state_par(epoch, idx_isb, allpars, update_info);
				}
			}
		}
		// remove old ISB parameter
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (t_gpar::is_sysbias(allpars[i].parType) && !update_info.exist(i + 1))
			{
				const string& site = allpars[i].site;
				const string& gsys = t_gpar::sysbias2cgsys(allpars[i].parType, allpars[i].channel);
				if (allpars[i].end < epoch) {
					update_info.add(i + 1);
				}
				else if (site_sys[site].find(gsys) == site_sys[site].end() && allpars[i].parType != par_type::GLO_IFB) {
					allpars[i].end = epoch - _intv;
					update_info.add(i + 1);
				}
			}
		}
	}

	void t_gupdatepar::update_isb_pars_new(const t_gtime& epoch, t_gallpar& allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo& update_info)
	{
		// ----------------- for test! ------------------
		if (_bd2_isb.empty()) {
			string site;
			double bd2_isb;
			fstream fs("isb_bds2");
			
			string line;
			while (getline(fs, line)) {
				stringstream ss(line);
				ss >> site >> bd2_isb;
				_bd2_isb[site] = bd2_isb;
			}
			if (fs) fs.close();
		}

		map<string, set<string>> site_sys;
		map<string, set<string>> sys_site;
		map<pair<string, string>, int> nobs_site_sys;
		for (const auto& obs : obsdata) {
			const auto& sat = obs.sat();
			// consider ISB between BDS2 and BDS3
			if (sat.substr(0, 1) == "C" && sat < "C17") {
				site_sys[obs.site()].insert("C2");
				sys_site["C2"].insert(obs.site());
				nobs_site_sys[make_pair<string, string>(obs.site(), "C2")]++;
			}
			else {
				site_sys[obs.site()].insert(sat.substr(0, 1));
				sys_site[sat.substr(0, 1)].insert(obs.site());
				nobs_site_sys[make_pair<string, string>(obs.site(), sat.substr(0, 1))]++;
			}
		}
		vector<string> ref_order = { "G", "E", "C", "C2", "J", "R" };
		set<par_type> isb_list = { par_type::BDS_ISB, par_type::BD2_ISB, par_type::GAL_ISB, par_type::GLO_ISB, par_type::QZS_ISB };
		// find CLK_SAT parameter to determine wether need a reference ISB for one system
		set<string> sys_est_clk;
		for (const auto& par : allpars.getAllPar()) {
			if (par.parType == par_type::CLK_SAT) {
				const string& sat = par.prn;
				if (sat.substr(0, 1) == "C" && sat < "C17") sys_est_clk.insert("C2");
				else sys_est_clk.insert(sat.substr(0, 1));
			}
		}
		set<string> new_ref_sys;
		for (const string& gsys : sys_est_clk) {
			set<string> sites = sys_site[gsys];
			if (sites.empty()) continue;
			if (sites.find(_sys_ref_site[gsys]) != sites.end()) continue;
			if (sites.find(_ref_clk) != sites.end()) {
				_sys_ref_site[gsys] = _ref_clk;
				new_ref_sys.insert(gsys);
			}
			else {
				string ref_site;
				int max_sys = 1;
				int max_nobs = 1;
				for (const auto& site : sites) {
					if (site_sys[site].size() >= max_sys && nobs_site_sys[{site, gsys}] > max_nobs) {
						ref_site = site;
						max_sys = site_sys[site].size();
						max_nobs = nobs_site_sys[{site, gsys}];
					}
				}
				_sys_ref_site[gsys] = ref_site;
				new_ref_sys.insert(gsys);
			}
		}

		const int idx_max = allpars.maxIndex();
		// int num = 0;
		for (const auto& iter : site_sys) {
			const auto& site = iter.first;
			const auto& sys_list = iter.second;
			map<par_type, int> isb_par_idx;
			for (const auto& par : isb_list) {
				int idx_isb = allpars.getParam(site, par, "");
				if (idx_isb > 0) isb_par_idx[par] = idx_isb;
			}
			// system num < 2, remove all ISB pars
			if (sys_list.size() < 2) {
				for (const auto& par : isb_par_idx) {
					if (allpars[par.second].end > epoch - _intv) allpars[par.second].end = epoch - _intv;
					update_info.add(par.second + 1);
				}
				continue;
			}
			// BDS_ISB denote BDS3 or BDS2, if the meaning of BDS_ISB change, update BDS_ISB
			bool bds_change = false;
			string ref_bds = "C";
			if (site_sys.find("C") == site_sys.end() && site_sys.find("C2") != site_sys.end()) ref_bds = "C2";
			if (_site_bds_ref.find(site) == _site_bds_ref.end()) {
				_site_bds_ref[site] = ref_bds;
			}
			else if (_site_bds_ref[site] != ref_bds && isb_par_idx.find(par_type::BD2_ISB) == isb_par_idx.end()) {
				_site_bds_ref[site] = ref_bds; bds_change = true;
			}
			// if the referece system changes, all ISB should be updated
			string ref_sys = "G";
			for (const auto& sys : ref_order) {
				if (sys_list.find(sys) != sys_list.end()) { ref_sys = sys;	break; }
			}
			bool ref_change = false;
			if (_site_ref_sys.find(site) == _site_ref_sys.end()) {
				_site_ref_sys[site] = ref_sys;
			}
			else if (_site_ref_sys[site] != ref_sys) {
				if (_site_ref_sys[site].substr(0, 1) != "C" || ref_sys.substr(0, 1) != "C" || bds_change) {
					for (const auto& par : isb_par_idx) {
						if (allpars[par.second].end > epoch - _intv) allpars[par.second].end = epoch - _intv;
						update_info.add(par.second + 1);
					}
					ref_change = true;
				}
				_site_ref_sys[site] = ref_sys;
			}

			set<string> sys_list_tmp;
			for (const string& sys : sys_list) sys_list_tmp.insert(sys.substr(0, 1));
			for (const auto& sys : sys_list_tmp) {
				if (sys == ref_sys) continue;
				par_type isb_type = par_type::NO_DEF;
				if (sys == "C") isb_type = par_type::BDS_ISB;
				else if (sys == "E") isb_type = par_type::GAL_ISB;
				else if (sys == "J") isb_type = par_type::QZS_ISB;
				else if (sys == "R") isb_type = par_type::GLO_ISB;
				if (isb_type == par_type::NO_DEF) continue;
				// init the ISB between BDS2 and BDS3
				if (sys == "C" && sys_list.find("C2") != sys_list.end() && sys_list.find("C") != sys_list.end() && 
					isb_par_idx.find(par_type::BD2_ISB) == isb_par_idx.end()) {
					t_gpar par_ISB = t_gpar(site, par_type::BD2_ISB, idx_max + 1 + _new_par, "");
					// -------- test ! -------
					if (_bd2_isb.find(site) != _bd2_isb.end()) {
						par_ISB.value(_bd2_isb[site]);
						par_ISB.apriori(0.001);
					}
					else {
						par_ISB.value(0.0);
						if (_sys_ref_site["C2"] == site) par_ISB.apriori(0.001);
						else par_ISB.apriori(3000.0);
					}
					par_ISB.setTime(epoch, _end_time);
					_new_par++;
					update_info.add(par_ISB);
				}
				// create new ISB parameters: 
				// 1. the reference site of one sys is changed, all ISB of this sys have to be reset
				// 2. the reference sys of one site is changed, all ISB of this site have to be reset
				// 3. no ISB parameter for corresponding observations, create the ISB for this site and sys
				if (new_ref_sys.find(sys) != new_ref_sys.end() || ref_change || 
					isb_par_idx.find(isb_type) == isb_par_idx.end() || (sys == "C" && bds_change)) {
					if (isb_par_idx.find(isb_type) != isb_par_idx.end() && !ref_change) {
						const int idx_isb = isb_par_idx[isb_type];
						if (allpars[idx_isb].end > epoch - _intv) allpars[idx_isb].end = epoch - _intv;
						update_info.add(idx_isb + 1);
					}
					t_gpar par_ISB = t_gpar(site, isb_type, idx_max + 1 + _new_par, "");
					par_ISB.value(0.0);
					if (_sys_ref_site[sys] == site) par_ISB.apriori(0.001);
					else par_ISB.apriori(3000.0);
					if (_sys_bias == SYSBIASMODEL::AUTO_CON || _sys_bias == SYSBIASMODEL::ISB_CON) {
						par_ISB.setTime(epoch, _end_time);
					}
					else if (_sys_bias == SYSBIASMODEL::AUTO_RWK || _sys_bias == SYSBIASMODEL::ISB_RWK) {
						par_ISB.setTime(epoch, epoch);
						set_par_state_mode(isb_type, 1, 24, 3000.0);
					}
					else if (_sys_bias == SYSBIASMODEL::AUTO_PWC || _sys_bias == SYSBIASMODEL::ISB_PWC) {
						t_gtime pwc_time = epoch + 3600 * _isb_intv - _intv;
						if (pwc_time > _end_time) pwc_time = _end_time;
						par_ISB.setTime(epoch, pwc_time);
						set_par_state_mode(isb_type, 1, _isb_intv, 3000.0);
					}
					else {
						par_ISB.setTime(epoch, epoch);
					}
					_new_par++;
					update_info.add(par_ISB);
				}
				// update existing ISB par
				else {
					if (_sys_bias == SYSBIASMODEL::AUTO_CON || _sys_bias == SYSBIASMODEL::ISB_CON) continue;
					const int idx_isb = isb_par_idx[isb_type];
					if (allpars[idx_isb].end > epoch - _intv) {
						if (_sys_bias == SYSBIASMODEL::AUTO_PWC || _sys_bias == SYSBIASMODEL::ISB_PWC) continue;
						allpars[idx_isb].end = epoch - _intv;
					}
					if (_sys_bias == SYSBIASMODEL::AUTO_WHIT || _sys_bias == SYSBIASMODEL::ISB_WHIT)
						_update_process_par(epoch, idx_isb, allpars, update_info);
					else
						_update_state_par(epoch, idx_isb, allpars, update_info);
				}
			}
		}
		// remove old ISB parameter
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			const string& site = allpars[i].site;
			if (t_gpar::is_sysbias(allpars[i].parType) && !update_info.exist(i + 1))
			{
				if (allpars[i].end < epoch) {
					update_info.add(i + 1);
				}
				else if (site_sys[site].find(t_gpar::sysbias2cgsys(allpars[i].parType)) == site_sys[site].end()) {
					// the ISB between BDS2 and BDS3 is considered to be constant
					if (allpars[i].parType != par_type::BD2_ISB) {
						allpars[i].end = epoch - _intv;
						update_info.add(i + 1);
					}
				}
			}
		}
	}

	t_gupdateparinfo t_gupdatepar::get_all_update_parameters(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata)
	{
		t_gupdateparinfo update_info;
		_new_par = 0;
		// update ISB
		update_isb_pars(epoch, allpars, obsdata, update_info);

		// update AMB
		_update_amb_pars(epoch,allpars,obsdata, update_info);

		// update Other pars
		_update_process_pars(epoch, allpars, update_info);

		return update_info;
	}

	void t_gupdatepar::_update_process_pars(const t_gtime& epoch, t_gallpar& allpars,t_gupdateparinfo & update_info)
	{

		// update process pars
		vector<t_gpar> parnew;
		for (unsigned int ipar = 0; ipar < allpars.parNumber(); ipar++)
		{
			if (t_gpar::is_amb(allpars[ipar].parType) || t_gpar::is_sysbias(allpars[ipar].parType)) continue;
			// according to time,lremove,update_info
			if (allpars[ipar].end < epoch && allpars[ipar].lremove && !update_info.exist(ipar+1))
			{	
				// judge state mode
				if (_state_mode.find(allpars[ipar].parType) == _state_mode.end())
				{
					_update_process_par(epoch, ipar, allpars, update_info);
				}
				else {
					_update_state_par(epoch, ipar, allpars, update_info);
				}
			}
		}
	}

	void t_gupdatepar::_update_state_par(const t_gtime& epoch, int update_id, t_gallpar& allpars, t_gupdateparinfo& update_info)
	{
		// record the new add par location
		vector<pair<int, t_gpar> > parnewlist;
		t_glsqEquationMatrix state_equ;
		t_glsqEquationMatrix state_equ_tmp;
		vector<t_gpar> equ_parlist;

		update_info.get(equ_parlist, state_equ_tmp);
		int newpar_location = allpars.parNumber() + 1 + equ_parlist.size();
		par_type par_type = allpars[update_id].parType;

		for (int i = 0; i < _state_mode[par_type].order; i++)
		{
			t_gpar par_temp = allpars[update_id + i];
			double intv = par_temp.end - par_temp.beg;
			t_gtime pwc_time = epoch + intv;
			if (pwc_time > _end_time) pwc_time = _end_time;
			par_temp.setTime(epoch, pwc_time);
			parnewlist.push_back(make_pair(update_id + i + 1, par_temp));
		}

		// create the state_equation
		for (int Row = 1; Row <= _state_mode[par_type].order; Row++)
		{
			vector<pair<int, double> > B;
			double P;
			// for B
			B.push_back(make_pair(update_id + Row, 1.0));
			for (int Col = Row; Col <= _state_mode[par_type].order; Col++) 
			{
				B.push_back(make_pair(newpar_location + Col - 1, -1.0));
			}

			// for P
			P = _state_mode[par_type].P(Row, Row);

			t_gobscombtype type;
			state_equ.add_equ(B, P, 0.0, "", "", type, false);
		}

		update_info.add(parnewlist, state_equ);

	}

	void t_gupdatepar::_update_process_par(const t_gtime & epoch, int update_id, t_gallpar & allpars, t_gupdateparinfo& update_info)
	{
		t_gpar par_tmp = allpars[update_id];
		double intv = par_tmp.end - par_tmp.beg;
		par_tmp.setTime(epoch, epoch + intv);
		update_info.add(update_id + 1, par_tmp);
	}

	void t_gupdatepar::_udpate_gps_rec_ifb_pars(const t_gtime & epoch, t_gallpar & allpars, const vector<t_gsatdata>& obsdata, t_gupdateparinfo & update_info)
	{
		// for debug
		double my_mjd = epoch.mjd();
		double my_sod = epoch.sod();


		std::set<string> mapRec;
		for (auto& data : obsdata)
		{
			auto  rec = data.site();
			if (mapRec.count(rec) > 0) continue;
			mapRec.insert(rec);
			// add new ifb par
			int idx_gps_ifb = allpars.getParam(rec, par_type::GPS_REC_IFB_C3, "");
			if (idx_gps_ifb < 0)
			{
				int idx_max = allpars[allpars.parNumber() - 1].index;
				t_gpar par_gps_rec_ifb = t_gpar(rec, par_type::GPS_REC_IFB_C3, idx_max + 1, "");
				par_gps_rec_ifb.value(0.0);
				par_gps_rec_ifb.setTime(epoch, epoch);
				par_gps_rec_ifb.apriori(9000);
				update_info.add(par_gps_rec_ifb);
			}
		}

		// remove old ifb par
		for (unsigned int i = 0; i < allpars.parNumber(); i++)
		{
			if (allpars[i].parType == par_type::GPS_REC_IFB_C3 && mapRec.count(allpars[i].site) == 0)
			{
				update_info.add(i + 1);
			}
		}
	}



}
