/**
*
* @verbatim
History
-1.0 ZhengHJ  2019-09-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		  gallrecover.cpp
* @brief	  Storage all recover info
*
* @author     ZhengHJ, Wuhan University
* @version	  1.0.0
* @date		  2019-09-25
*
*/

#include "gallrecover.h"

namespace great
{
	t_gallrecover::t_gallrecover() :
		t_gdata()
	{
		_type = t_gdata::ALLRECOVER;
	}


	t_gallrecover::~t_gallrecover()
	{
		for(t_grecover_data* data : _recover_data)
		{
			delete data;
			data = nullptr;
		}
		_time_equmap.clear();
		_time_parmap.clear();
		_site_sat_time_equmap.clear();
		_type_parmap.clear();
	}

	void t_gallrecover::add_recover_equation(const t_grecover_equation& recover_equ)
	{
		t_grecover_equation* equ_data = new t_grecover_equation(recover_equ);

		_add_common_data(equ_data);

		_recover_head.site_list.insert(recover_equ.site_name);                              
		_recover_head.sat_list.insert(recover_equ.sat_name);


		// default set && add
		if (_time_equmap.find(recover_equ.time) == _time_equmap.end()) {
			_time_equmap[recover_equ.time] = vector<t_grecover_equation*>();
		}
		// add equ map
		_time_equmap[recover_equ.time].push_back(equ_data);

		if (_site_sat_time_equmap.find(recover_equ.site_name) == _site_sat_time_equmap.end()) {
			_site_sat_time_equmap[recover_equ.site_name] = t_map_sat_equ();
		}
		t_map_sat_equ& temp_sat_equ = _site_sat_time_equmap[recover_equ.site_name];
		if (temp_sat_equ.find(recover_equ.sat_name) == temp_sat_equ.end()) {
			temp_sat_equ[recover_equ.sat_name] = t_map_time_equ();
		}
		t_map_time_equ& temp_time_equ = temp_sat_equ[recover_equ.sat_name];
		if (temp_time_equ.find(recover_equ.time) == temp_time_equ.end()) {
			temp_time_equ[recover_equ.time] = vector<t_grecover_equation*>();
		}
		temp_time_equ[recover_equ.time].push_back(equ_data);

	}

	void t_gallrecover::add_recover_par(const t_grecover_par & recover_par)
	{
		t_grecover_par* par_data = new t_grecover_par(recover_par);

		_add_common_data(par_data);

		if (_time_parmap.find(recover_par.get_recover_time()) == _time_parmap.end())
		{
			_time_parmap[recover_par.get_recover_time()] = vector<t_grecover_par*>();
			_time_parmap_pce[recover_par.get_recover_time()] = vector<t_grecover_par>();
		}
		_time_parmap[recover_par.get_recover_time()].push_back(par_data);
		_time_parmap_pce[recover_par.get_recover_time()].push_back(*par_data);
		
		if (_type_parmap.find(recover_par.par.parType) == _type_parmap.end())
		{
			_type_parmap[recover_par.par.parType]= vector<t_grecover_par*>();
		}
		_type_parmap[recover_par.par.parType].push_back(par_data);
	}

	void t_gallrecover::get_clkdata(t_gallprec& clkdata, t_gallprec::clk_type type)
	{
		set<t_gallprec::clk_type> type_list;
		if (type != t_gallprec::UNDEF) {
			type_list.insert(type);
		}
		else {
			type_list.insert(t_gallprec::AR);
			type_list.insert(t_gallprec::AS);
		}
		
		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) {
				string obj = "";
				if (par->par.parType == par_type::CLK && type_list.find(t_gallprec::AR) != type_list.end() ){
					obj = par->par.site;
				}
				else if (par->par.parType == par_type::CLK_SAT && type_list.find(t_gallprec::AS) != type_list.end() ) {
					obj = par->par.prn;
				}
				else {
					continue;
				}

				if (par->correct_value == 0.0) {
					continue;
				}

				double clk[3] = { 0.0,0.0,0.0 };
				double var[3] = { 0.0,0.0,0.0 };

				clk[0] = (par->par.value() + par->correct_value) / CLIGHT;
				clkdata.addclk(obj, epoch, clk, var);
			}
		}
	}

	void t_gallrecover::get_stadata(t_gallsta& stadata)
	{
		map<string, double> dx;
		map<string, double> dy;
		map<string, double> dz;

		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) {
				string obj = "";
				if (par->par.parType == par_type::CRD_X) {
					obj = par->par.site;
					dx[obj] = par->correct_value;
				}
				else if (par->par.parType == par_type::CRD_Y) {
					obj = par->par.site;
					dy[obj] = par->correct_value;
				}
				else if (par->par.parType == par_type::CRD_Z) {
					obj = par->par.site;
					dz[obj] = par->correct_value;
				}
				else {
					continue;
				}

			}
		}

		for (auto iter = dx.begin(); iter != dx.end(); iter++)
		{
			stadata.addsta(iter->first, iter->second, dy[iter->first], dz[iter->first]);
		}

	}

	void t_gallrecover::get_iondata(vector<t_tuple_ion>& ion_data)
	{
		for (auto iter = _time_parmap.begin(); iter != _time_parmap.end(); iter++)
		{
			t_gtime epoch = iter->first;
			for (auto par : iter->second) 
			{
				string type = "";
				if (par->par.parType == par_type::SION)
				{
					type = "SION";
				}
				else if (par->par.parType == par_type::VION)
				{
					type = "VION";
				}
				else
				{
					continue;
				}

				if (double_eq(par->correct_value, 0.0)) continue;

				double value = 0.0;
				double sigma = 0.5;
				
				value = (par->par.value() + par->correct_value);
				auto tmp = make_tuple(type, par->par.site, par->par.prn,
					value, sigma, par->par.beg, par->par.end);
				ion_data.push_back(tmp);
			}
		}
	}

	vector<t_grecover_equation> t_gallrecover::get_first_phase_recover_equation(string site, string sat, string freq)
	{
		string trim_freq = (trim(freq));

		map<t_gobscombtype, int> phase_order = 
		{
			{t_gobscombtype("LC12"),13},
			{t_gobscombtype("LC13"),12},
			{t_gobscombtype("LC14"),11},
			{t_gobscombtype("LC15"),10},
			{t_gobscombtype("L1"),9},
			{t_gobscombtype("L2"),8},
			{t_gobscombtype("L3"),7},
			{t_gobscombtype("L4"),6},
			{t_gobscombtype("L5"),5},
		};

		vector<t_grecover_equation> first_phase_recover_equation;

		if (phase_order.find(t_gobscombtype(trim_freq)) != phase_order.end())
		{
			phase_order[t_gobscombtype(trim_freq)] = 14;
		}
		else
		{
			throw runtime_error("check your xml file, the edit res freq should be LC1X or LCX (X means 1-5)");
		}

		for(auto time_equ : _site_sat_time_equmap[site][sat]) 
		{		
			t_grecover_equation first_equ(time_equ.first,site,sat);

			bool find = false;
			int  max_order = 0;
			for(auto equ : time_equ.second)
			{
				if (!equ->obstype.is_phase()) 
				{
					continue;
				}

				string trim_equ_freq = trim(("L" + gfreqseq2str(equ->obstype.getFreq_1()) + gfreqseq2str(equ->obstype.getFreq_2())));
				if (trim_freq != trim_equ_freq)
				{
					continue;
				}

				if (phase_order[equ->obstype] > max_order)
				{
					max_order = phase_order[equ->obstype];
					first_equ = *equ;
					find = true;
					continue;
				}
			}

			if(find)  first_phase_recover_equation.push_back(first_equ);

		}
		return first_phase_recover_equation;
	}

	bool t_gallrecover::get_first_phase_recover_equation(string site, string sat, vector<t_grecover_equation>& equ, string freq)
	{
		string trim_freq = (trim(freq));

		map<t_gobscombtype, int> phase_order =
		{
			{t_gobscombtype("LC12"),13},
			{t_gobscombtype("LC13"),12},
			{t_gobscombtype("LC14"),11},
			{t_gobscombtype("LC15"),10},
			{t_gobscombtype("L1"),9},
			{t_gobscombtype("L2"),8},
			{t_gobscombtype("L3"),7},
			{t_gobscombtype("L4"),6},
			{t_gobscombtype("L5"),5},
		};

		vector<t_grecover_equation> first_phase_recover_equation;

		bool freq_find = false;
		int  max_order = 0;
		if (phase_order.find(t_gobscombtype(trim_freq)) != phase_order.end())
		{
			phase_order[t_gobscombtype(trim_freq)] = 14;
			max_order = 14;
			freq_find = true;
		}
		else
		{
			throw runtime_error("check your xml file, the edit res freq should be LC1X or LCX (X means 1-5)");
		}

		for (auto time_equ : _site_sat_time_equmap[site][sat])
		{
			t_grecover_equation first_equ(time_equ.first, site, sat);

			bool find = false;
			for (auto equ : time_equ.second)
			{
				if (!equ->obstype.is_phase())
				{
					continue;
				}

				if (phase_order[equ->obstype] == max_order)
				{
					if (find == false) first_equ = *equ;
					find = true;
					continue;
				}

				if (phase_order[equ->obstype] > max_order)
				{
					max_order = phase_order[equ->obstype];
					first_equ = *equ;
					find = true;
					continue;
				}
			}

			if (find)  first_phase_recover_equation.push_back(first_equ);

		}

		if (first_phase_recover_equation.empty()) return false;
		equ = first_phase_recover_equation;

		return true;
	}

	vector<t_grecover_par> t_gallrecover::get_recover_par(par_type parType)
	{
		vector<t_grecover_par> tmp;
		for (auto type_par : _type_parmap[parType])
		{
			t_grecover_par par(type_par->par, type_par->correct_value);
			tmp.push_back(par);
		}
		return tmp;
	}


	t_gtime t_gallrecover::get_beg_time() const
	{
		return _recover_head.get_beg_time();
	}

	t_gtime t_gallrecover::get_end_time() const
	{
		return _recover_head.get_end_time();
	}

	t_gtime t_gallrecover::get_equ_beg_time() const
	{
		return _time_equmap.begin()->first;
	}

	double t_gallrecover::get_interval() const
	{
		return _recover_head.interval;
	}

	double t_gallrecover::get_sigma0() const
	{
		return _recover_head.sigma0;
	}

	void t_gallrecover::set_interval(double intv)
	{
		_recover_head.interval = intv;
	}

	void t_gallrecover::set_sigma0(double sigma0)
	{
		_recover_head.sigma0 = sigma0;
	}

	const vector<t_grecover_data*>& t_gallrecover::get_all_recover_data() const
	{
		return _recover_data;
	}

	double t_gallrecover::get_recover_data_value(par_type parType)
	{
		vector<t_grecover_data*>::iterator it;
		for (it = _recover_data.begin(); it != _recover_data.end(); it++)
		{
			double correct_value;
			double par_value;
			if (parType == dynamic_cast<t_grecover_par*>(*it)->par.parType)
			{
				par_value = dynamic_cast<t_grecover_par*>(*it)->par.value();
				correct_value = dynamic_cast<t_grecover_par*>(*it)->correct_value;
				if (double_eq(par_value, 0.0) || double_eq(correct_value, 0.0)) return 0;
				return par_value + correct_value;
			}
		}
		return 0;
	}

	const t_map_time_par& t_gallrecover::get_map_time_par() const
	{
		return _time_parmap;
	}

	const t_map_time_equ& t_gallrecover::get_map_time_equ() const
	{
		return _time_equmap;
	}
	const t_map_site_equ& t_gallrecover::get_map_site_equ() const
	{
		return _site_sat_time_equmap;
	}

	set<string> t_gallrecover::get_sat_list() const
	{
		return _recover_head.sat_list;
	}

	set<string> t_gallrecover::get_site_list() const
	{
		return _recover_head.site_list;
	}

	set<t_gtime> t_gallrecover::get_time_list() const
	{
		return _recover_head.time_list;
	}

	vector<t_gtime> t_gallrecover::get_all_obs_time_list() const
	{
		t_gtime beg = _time_equmap.begin()->first;
		t_gtime end = _time_equmap.rbegin()->first;

		int nepoch = floor((end - beg) / _recover_head.interval) + 1;
		vector<t_gtime> ans;
		for (int i = 0; i < nepoch; i++) {
			ans.push_back(beg + _recover_head.interval * i);
		}
		return ans;
	}

	t_gallpar t_gallrecover::get_all_pars()
	{
		t_gallpar tmp;

		t_gtime beg  = get_beg_time();
		t_gtime end  = get_end_time();
		double  intv = get_interval();
		while (beg <= end)
		{
			auto par_end = _time_parmap[beg].end();
			auto par_crt = _time_parmap[beg].begin();
			for (par_crt; par_crt != par_end; par_crt++)
			{
				// 加上改正数
				auto recover_par = (*par_crt)->par;
				recover_par.value(recover_par.value() + (*par_crt)->correct_value);
				tmp.addParam(recover_par);
			}
		}
		
		return  tmp;
	}

	map<t_gtime, int> t_gallrecover::get_sat_number_per_epo() const
	{
		map<t_gtime, int> nsat;
		for (auto it : _time_equmap)
		{
			nsat[it.first] = it.second.size() / 2;
		}
		return nsat;
	}

	void t_gallrecover::_add_common_data(t_grecover_data * data)
	{
		_recover_data.push_back(data);
		_recover_head.time_list.insert(data->get_recover_time());
	}

}