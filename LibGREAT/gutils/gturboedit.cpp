/**
*
* @verbatim
	History
	 -1.0 GREAT	    2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  gturboedit.cpp
* @brief		  CRS coordinate transform to ACR coordinate
*
* @author         GREAT, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/

#include "gutils/gturboedit.h"
#include "gset/gsetrec.h"
#include "gset/gsetinp.h"
#include "gutils/gfileconv.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <set>
#include "gio/gfile.h"
#include "gcoders/ambflag.h"

#ifdef WIN32
	#include <io.h>
#endif // WIN32

using namespace std;
using namespace gnut;

namespace great
{
	great::t_gturboedit::t_gturboedit()
	{
	}

	great::t_gturboedit::t_gturboedit(t_gsetbase * gset, t_glog * glog, int index) : t_gcycleslip(gset, glog)
	{
		_gset = gset;
		_glog = glog;

		_index = index;
		t_gtime beg_time = dynamic_cast<t_gsetgen*>(gset)->beg();

		set<string> rec_list = dynamic_cast<t_gsetgen*>(gset)->rec_all(); 

		bool lite_turboedit = dynamic_cast<t_gsetturboedit*>(gset)->liteMode();
		if (_slip_model == SLIPMODEL::DEF_DETECT_MODEL) lite_turboedit = true;

		if (!lite_turboedit)
		{
			_read_logfile(rec_list, beg_time, index);
		}
	}

	great::t_gturboedit::~t_gturboedit()
	{
	}

	set<string> great::t_gturboedit::get_sitelist_of_logfile() const
	{
		set<string> site_list;
		for (auto iter = _amb_info_file_exist.begin(); iter != _amb_info_file_exist.end(); iter++)
		{
			if (iter->second)
			{
				site_list.insert(iter->first);
			}
		}
		return site_list;
	}

	void great::t_gturboedit::merge_logfile_exist(const map<string, bool>& logfiles)
	{
		for (const auto& file : logfiles)
		{
			if (_amb_info_file_exist.count(file.first) == 0)
			{
				_amb_info_file_exist[file.first] = false;
			}
			else
			{
				_amb_info_file_exist[file.first] = (file.second && _amb_info_file_exist[file.first]);
			}
		}
	}

	void great::t_gturboedit::_read_logfie(const set<string>& rec, const t_gtime& epoch, int index)
	{
		multimap<IFMT, string> inp_xml = dynamic_cast<t_gsetinp*>(_gset)->inputs_all();
		map<string, string> inp_log_xml;
		for (auto& inp : inp_xml)
		{
			if (inp.first == IFMT::AMBFLAG12_INP && index == 2)
			{
				string _file_path = inp.second.substr(inp.second.find_first_of("/") + 2);
				string _file_name = file_name(inp.second);
				int pos = -1;
				if (index == 2) pos = _file_name.find_last_not_of("log");
				else continue;

				if (pos < 12) continue;

				string rec_low = _file_name.substr(pos - 12, 4);
				inp_log_xml.insert(make_pair(rec_low, _file_path));
			}
		}

		fstream amb_info_file;
		for (auto rec_iter = rec.begin(); rec_iter != rec.end(); rec_iter++)
		{
			stringstream log_suffix;

			if (index == 2)
			{
				log_suffix << setw(3) << setfill('0') << epoch.doy() << "0."
					<< setw(2) << setfill('0') << epoch.yr() << "o.log";
			}

			string rec_name_low;
			transform((*rec_iter).begin(), (*rec_iter).end(), back_inserter(rec_name_low), ::tolower);

			// first find the log file in xml input setting
			// if no log file in the path of xml file, find it in log_tb path
			string rec_amb_info_file_name = inp_log_xml[rec_name_low];
			if (ACCESS(rec_amb_info_file_name.c_str(), 0) != 0)
			{
				rec_amb_info_file_name = "log_tb/" + rec_name_low + log_suffix.str();
			}

			amb_info_file.open(rec_amb_info_file_name, ios::in);
			if (!amb_info_file.is_open())
			{
				cout << rec_amb_info_file_name << ":can't open!" << endl;
				_amb_info_file_exist[*rec_iter] = false;
				continue;
			}

			if (_slip_model == SLIPMODEL::DEF_DETECT_MODEL)
			{
				_amb_info_file_exist[*rec_iter] = false;
				continue;
			}

			// read
			string line_txt;
			getline(amb_info_file, line_txt);
			line_txt = line_txt.substr(27);
			stringstream input(line_txt);
			double sod, intv;
			int temp, jd;
			input >> jd >> sod >> temp >> intv;

			t_gtime beg_time(t_gtime::GPS);
			beg_time.from_mjd(jd, int(sod), sod - int(sod));
			line_txt.clear();
			while (!amb_info_file.eof())
			{
				getline(amb_info_file, line_txt);
				if (line_txt.size() <= 1)
				{
					line_txt.clear();
					continue;
				}
				if (line_txt[0] == '%' && line_txt[1] == 'M')
				{
					stringstream oss(line_txt);
					string str1, str2, str3, str4, str5, str6;
					int max_amb;
					oss >> str1 >> str2 >> str3 >> str4 >> str5 >> str6 >> max_amb;
					_active_amb[*rec_iter] = max_amb;
					line_txt.clear();
					continue;
				}

				if (line_txt[0] < 'A' || line_txt[0] > 'Z')
				{
					line_txt.clear();
					continue;
				}
				stringstream os(line_txt);
				string identify;
				string sat_name;
				int beg_idx, end_idx, amb_flag;
				os >> identify >> sat_name >> beg_idx >> end_idx >> amb_flag;

				t_gtime beg, end;
				beg = beg_time + (beg_idx - 1) * intv;
				end = beg_time + (end_idx - 1) * intv;

				if (amb_flag == 1 || amb_flag == 2 || identify == "AMB")
				{
					_cycle_flag[*rec_iter][sat_name].push_back(make_pair(beg, end));
				}
				else
				{
					_cycle_flag_unused[*rec_iter][sat_name].push_back(make_pair(beg, end));
				}

				line_txt.clear();
			}

			amb_info_file.close();

			_amb_info_file_exist[*rec_iter] = true;
		}
	}

	void great::t_gturboedit::_read_logfile(const set<string>& rec, const t_gtime& epoch, int index)
	{
		if (!_gambflag) _gambflag = make_shared<t_gallambflag>(t_gdata::AMBFLAG);

		// first, get all log file write in xml file
		multimap<IFMT, string> inp_xml = dynamic_cast<t_gsetinp*>(_gset)->inputs_all();
		map<string, string> inp_log_xml;
		for (auto& inp : inp_xml) {
			if (inp.first == IFMT::AMBFLAG12_INP && index == 2) {
				string _file_path = inp.second.substr(inp.second.find_first_of("/") + 2);
				string _file_name = file_name(inp.second);
				int pos;
				if (index == 2) pos = _file_name.find_last_not_of("log");
				string rec_low = _file_name.substr(pos - 12, 4);
				inp_log_xml.insert(make_pair(rec_low, _file_path));
			}
		}
		t_gdata* gdata = nullptr;
		for (auto site_iter = rec.begin(); site_iter != rec.end(); ++site_iter) {
			stringstream log_suffix;
			if (index == 2) { log_suffix << setw(3) << setfill('0') << epoch.doy() << "0." << setw(2) << setfill('0') << epoch.yr() << "o.log"; }
			else if (index == 3) { log_suffix << setw(3) << setfill('0') << epoch.doy() << "0." << setw(2) << setfill('0') << epoch.yr() << "o.log13"; }
			else if (index == 4) { log_suffix << setw(3) << setfill('0') << epoch.doy() << "0." << setw(2) << setfill('0') << epoch.yr() << "o.log14"; }
			else if (index == 5) { log_suffix << setw(3) << setfill('0') << epoch.doy() << "0." << setw(2) << setfill('0') << epoch.yr() << "o.log15"; }

			string site_low;
			transform((*site_iter).begin(), (*site_iter).end(), back_inserter(site_low), ::tolower);

			// first find the log file in xml input setting
			// if no log file in the path of xml file, find it in log_tb path
			string site_file_name = inp_log_xml[site_low];
			if (ACCESS(site_file_name.c_str(), 0) != 0) 
			{
				//cout << "can not find file : " + site_file_name << endl;
				site_file_name = "log_tb/" + site_low + log_suffix.str();
			}

			if (ACCESS(site_file_name.c_str(), 0) != 0)
			{
				_amb_info_file_exist[*site_iter] = false;
				continue;
			}

			if (_slip_model == SLIPMODEL::DEF_DETECT_MODEL)
			{
				_amb_info_file_exist[*site_iter] = false;
				continue;
			}

			gdata = _gambflag.get();
			t_gcoder* gcoder = new t_ambflag(_gset, "", 4096);
			t_gio* gio = new t_gfile;
			string path("file://" + site_file_name);
			gio->glog(_log);
			gio->path(path);

			// Put the file into gcoder
			gcoder->clear();
			gcoder->path(path);
			gcoder->glog(_log);

			// Put the data container into gcoder
			gcoder->add_data("ID0", gdata);
			gio->coder(gcoder);
			gio->run_read();

			// Delete 
			delete gio;
			delete gcoder;

			_amb_info_file_exist[*site_iter] = true;
		}
	}

	bool great::t_gturboedit::use_of_obs(const string& site, const string& prn, const t_gtime& time) {
		if (!_gambflag) return false;
		string lower_site(site);
		transform(site.begin(), site.end(), lower_site.begin(), ::tolower);
		if (_gambflag->getAllAmbFlag().find(lower_site) == _gambflag->getAllAmbFlag().end())
			return false;
		auto ambflag = _gambflag->getOneAmbFlag(lower_site);
		int pos = 0;
		return ambflag.isValid(prn, time, pos);
	}

	int great::t_gturboedit::num_of_amb_arc(const string& site, const string& prn, const t_gtime& time)
	{
		if (!_gambflag) return -1;
		string lower_site(site);
		transform(site.begin(), site.end(), lower_site.begin(), ::tolower);
		if (_gambflag->getAllAmbFlag().find(lower_site) == _gambflag->getAllAmbFlag().end())
			return -1;
		auto ambflag = _gambflag->getOneAmbFlag(lower_site);
		return ambflag.get_amb_pos(prn, time, _apply_carrier_range) + 1;
	}

	bool great::t_gturboedit::cycle_slip(const t_gsatdata& obsdata, const t_gtime& time)
	{
		// change some codes
		string sat = obsdata.sat();
		string rec = obsdata.site();

		if (_amb_info_file_exist[rec])
		{
			if (_amb_flag.find(rec) != _amb_flag.end() && _amb_flag[rec].find(sat) != _amb_flag[rec].end() && _amb_flag[rec][sat] == 0)
			{
				throw runtime_error("Can't judge cycle slip unless already in one amb arc");
			}
			int flag_now = num_of_amb_arc(rec, sat, time);
			if (flag_now && flag_now != _amb_flag[rec][sat])
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		return false;
	}

	bool great::t_gturboedit::cycle_slip123(t_gsatdata& obsdata, const t_gtime& time)
	{
		// change some codes
		GSYS gs = obsdata.gsys();
		string sat = obsdata.sat();
		string rec = obsdata.site();

		if (_amb_info_file_exist[rec])
		{
			if (_amb_flag[rec][sat] == 0)
			{
				//throw runtime_error("Can't judge cycle slip unless already in one amb arc");
			}
			int flag_now = num_of_amb_arc(rec, sat, time);
			if (flag_now && flag_now != _amb_flag[rec][sat])
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		return false;
	}

	int great::t_gturboedit::get_active_amb(string site)
	{
		return _active_amb[site];
	}

	void great::t_gturboedit::set_active_amb(string site, int active_num)
	{
		if (_slip_model != SLIPMODEL::TURBO_EDIT) _slip_model = SLIPMODEL::TURBO_EDIT;
		_amb_info_file_exist[site] = true;
		_active_amb[site] = active_num;
	}

	void great::t_gturboedit::set_amb_flag(const string& rec, const string& sat, int flag)
	{
		_amb_flag[rec][sat] = flag;
		_new_amb[rec][sat] = true;
	}

	int great::t_gturboedit::get_amb_flag(const string& rec, const string& sat)
	{
		return _amb_flag[rec][sat];
	}

	bool great::t_gturboedit::new_amb(const string& rec, const string& sat)
	{
		if (_new_amb.find(rec) == _new_amb.end())
		{
			return false;
		}
		if (_new_amb[rec].find(sat) == _new_amb[rec].end())
		{
			return false;
		}
		return _new_amb[rec][sat];
	}

	void great::t_gturboedit::set_new_amb(const string& rec, const string& sat, bool isNew)
	{
		_new_amb[rec][sat] = isNew;
	}

	t_gtime great::t_gturboedit::get_crt_amb_end(const string& site, const string& sat)
	{
		if (_amb_flag[site][sat] <= 0)
		{
			string tmp = "sat : " + sat + " rec : " + site + " freq " + int2str(_index);
			throw runtime_error("Crt amb isn't exist!!! " + tmp);
		}

		string lower_site(site);
		transform(site.begin(), site.end(), lower_site.begin(), ::tolower);
		auto ambflag = _gambflag->getOneAmbFlag(lower_site);
		auto one_ambflag = ambflag.getSatAmbFlag(sat).at(_amb_flag[site][sat] - 1);
		return ambflag.epoch2time(one_ambflag->end_epo);
	}

	void great::t_gturboedit::add_ambflag(string site, string sat, string description, t_gtime beg, t_gtime end)
	{
		if (description == "AMB")
		{
			_cycle_flag[site][sat].push_back(make_pair(beg, end));
		}
		else
		{
			_cycle_flag_unused[site][sat].push_back(make_pair(beg, end));
		}

	}
}

