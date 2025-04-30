/**
*
* @verbatim
History
-1.0 Hongjie Zheng  2019-09-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		recover.cpp
* @brief	The base class used to decode and encode recover file information.
*
* @author   Hongjie Zheng, Wuhan University
* @version	1.0.0
* @date		2019-09-25
*
*/

#include <gcoders/recover.h>
#include <gall/gallrecover.h>
#include <gdata/grecoverdata.h>

namespace great
{
	const string t_resfile::END_OF_HEADER = "##END OF HEADER";
	const string t_resfile::TIME_HEADER = "##Time&Interval :";
	const string t_resfile::SIGMA_HEADER = "##Sigma:";
	const string t_resfile::SAT_HEADER = "##SAT:";
	const string t_resfile::SITE_HEADER = "##STA:";
	const string t_resfile::OBS_DATA = "RES:=";
	const string t_resfile::PAR_DATA = "PAR:=";

	t_resfile::t_resfile(t_gsetbase * s, string version, int sz):
		t_gcoder(s,version,sz),
		_recover_data(nullptr)
	{
	}

	t_resfile::~t_resfile()
	{
	}

	int t_resfile::decode_head(char * buff, int sz, vector<string>& errmsg)
	{
		_mutex.lock();

		if (!_recover_data) {
			if (_log) {
				_log->comment(0, "t_resfile::decode_head", " ERROR: Have no storeage allrecover for decoding resfile");
			}
			_mutex.unlock();
			return -1;
		}

		int tmpsize = 0;
		int linesize = 0;
		string line;
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

		while ((linesize = t_gcoder::_getline(line,tmpsize) )>= 0) {
			tmpsize += linesize;
			int idx = 0;
			if ( line.find(t_resfile::END_OF_HEADER) !=-1 ) {
				t_gcoder::_consume(tmpsize);
				_mutex.unlock();
				return -1;
			}
			else if (line.find(t_resfile::TIME_HEADER) != -1) {
				line = line.substr(t_resfile::TIME_HEADER.length());
				string temp_ymd, temp_hms;
				double interval;
				stringstream templine(line);
				templine >> temp_ymd >> temp_hms >> interval;
				_recover_data->set_interval(interval);
			}
			else if (line.find(t_resfile::SIGMA_HEADER) != -1) {
				line = line.substr(t_resfile::SIGMA_HEADER.length());
				double sigma;
				stringstream templine(line);
				templine >> sigma;
				_recover_data->set_sigma0(sigma);
			}
			else {
				continue;
			}
		}

		t_gcoder::_consume(tmpsize);
		_mutex.unlock();
		return tmpsize;
	}

	int t_resfile::decode_data(char * buff, int sz, int & cnt, vector<string>& errmsg)
	{
		_mutex.lock();
		if (!_recover_data) {
			if (_log) {
				_log->comment(0, "t_resfile::decode_head", " ERROR: Have no storeage allrecover for decoding resfile");
			}
			_mutex.unlock();
			return -1;
		}

		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

		int linesize = 0;
		int tmpsize = 0;
		string line;

		while ((linesize = t_gcoder::_getline(line, tmpsize)) >= 0)
		{
			tmpsize += linesize;
			int idx = 0;
			if ((idx = line.find(t_resfile::OBS_DATA)) != -1) {

				_recover_data->add_recover_equation(strline2recover_equation(line.substr(idx + t_resfile::OBS_DATA.length())));
			}
			else if ((idx = line.find(t_resfile::PAR_DATA)) != -1) {

				_recover_data->add_recover_par(strline2recover_par(line.substr(idx + t_resfile::OBS_DATA.length())));
			}
			else {
				continue;
			}
		}

		tmpsize = t_gcoder::_consume(tmpsize);
		_mutex.unlock();
		return tmpsize;

	}

	int t_resfile::encode_head(char * buff, int sz, vector<string>& errmsg)
	{
		_mutex.lock();


		if (!_recover_data) {
			if (_log) {
				_log->comment(0, "t_resfile::encode_head", " ERROR: Have no data for encode in resfile encoder");
			}
			_mutex.unlock();
			return -1;
		}

		if (_ss_position == 0)
		{
			_ss << t_resfile::TIME_HEADER
				<< setiosflags(ios::right)
				<< setw(30) << _recover_data->get_beg_time().str() 
				<< setw(15) << _recover_data->get_interval() 
				<< endl;

			_ss << t_resfile::SIGMA_HEADER
				<< setw(15) << fixed << setprecision(3) << _recover_data->get_sigma0() << endl;

			_ss << t_resfile::SITE_HEADER;
			int count = 0;
			for (string site : _recover_data->get_site_list()) {
				_ss << setw(5) << setiosflags(ios::right) << site;
				if (++count == 10) {
					_ss << endl << t_resfile::SITE_HEADER;
					count = 0;
				}
			}
			_ss << endl;

			_ss << t_resfile::SAT_HEADER;
			count = 0;
			for (string sat : _recover_data->get_sat_list()) {
				_ss << setw(4) << setiosflags(ios::right) << sat;
				if (++count == 10) {
					_ss << endl << t_resfile::SAT_HEADER;
					count = 0;
				}
			}
			_ss << endl;

			_ss << t_resfile::END_OF_HEADER << endl;

		}
		int size = _fill_buffer(buff, sz);
		_mutex.unlock();
		return size;
	}

	int t_resfile::encode_data(char * buff, int sz, int & cnt, vector<string>& errmsg)
	{
		_mutex.lock();

		if (_ss_position == 0)
		{
			for(t_grecover_data* recover_record : _recover_data->get_all_recover_data()) 
			{
				_ss << recover_record->convert2strline();
			}
		}
		int size = _fill_buffer(buff, sz);
		_mutex.unlock();
		return size;
	}

	void t_resfile::_add_data(string id, t_gdata * data)
	{
		if (data->id_type() == t_gdata::ALLRECOVER) {
			_recover_data = dynamic_cast<t_gallrecover*>(data);
		}
	}


	 t_grecover_par strline2recover_par(string line)
	{
		stringstream templine(line);
		string str_partype;
		string time_beg1, time_beg2, time_end1, time_end2;
		double inital_value, correct_value;
		templine >> str_partype >> time_beg1 >> time_beg2 >> time_end1 >> time_end2 >> inital_value >> correct_value;


		t_gpar par = str2gpar(str_partype);
		t_gtime beg_time(GPS), end_time(GPS);
		beg_time.from_str("%Y-%m-%d %H:%M:%S", time_beg1 + " " + time_beg2);
		end_time.from_str("%Y-%m-%d %H:%M:%S", time_end1 + " " + time_end2);
		par.beg = beg_time;
		par.end = end_time;
		par.value(inital_value);
		t_grecover_par recover_par(par, correct_value);
		return recover_par;

	}

	t_grecover_equation strline2recover_equation(string line)
	{
		stringstream templine(line);

		string time1, time2;
		string site, sat, str_obstype;
		double weight, residual;
		int is_newamb;
		templine >> time1 >> time2 >> is_newamb >>  site >> sat >> str_obstype >> weight >> residual;
		t_gtime time;
		time.from_str("%Y-%m-%d %H:%M:%S", time1 + " " + time2);

		t_grecover_equation recover_equ(time, site, sat);
		recover_equ.set_recover_equation(t_gobscombtype(str_obstype), make_pair(weight, residual),is_newamb);

		return recover_equ;
	}

}