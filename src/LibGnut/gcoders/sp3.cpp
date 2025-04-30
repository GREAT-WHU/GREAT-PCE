/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>

#include "gall/gallprec.h"
#include "gcoders/sp3.h"
#include "gutils/gtriple.h"
#include "gutils/gtypeconv.h"
 
using namespace std;

namespace gnut {

/* ----------
 * constructor
 */
t_sp3::t_sp3( t_gsetbase* s, string version, int sz )
  : t_gcoder( s, version, sz )
{
  gtrace("t_sp3::constructor");

  _start.tsys(t_gtime::GPS);
  _lastepo.tsys(t_gtime::GPS);
  _nrecord = -1;
  _nrecmax = -1;
}


/* ----------
 * SP3 header
 */
int t_sp3::decode_head(char* buff, int sz, vector<string>& errmsg)
{
	gtrace("t_sp3::decode_head");

	_mutex.lock();

	if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

	string tmp;
	int consume = 0;
	int tmpsize = 0;

	while ((tmpsize = t_gcoder::_getline(tmp)) >= 0) {

		consume += tmpsize;

		// first information
		if (tmp.substr(0, 1) == "#")
		{
			// first line
			if (tmp.substr(1, 1) == "a" ||  // 60-columns
				tmp.substr(1, 1) == "c"  || tmp.substr(1, 1) == "d"    // 80-columns
				) {

				_version = tmp.substr(1, 1);
				_start.from_str("%Y %m %d  %H %M %S", tmp.substr(3, 30));
				_nepochs = str2int(tmp.substr(32, 7));
				_orbrefs = tmp.substr(46, 5);
				_orbtype = tmp.substr(52, 3);
				if (tmp.size() > 60)
				{
					_agency = tmp.substr(56, 4);
				}
				else
				{
					_agency = "";
				}
				// second line
			}
			else if (tmp.substr(1, 1) == "#") {
				_orbintv = (long)str2dbl(tmp.substr(24, 14)); // [sec]
				if (_log) {
					ostringstream ltmp;
					ltmp << "start time = " << _start.str(" %Y-%m-%d %H:%M:%S")
						<< "  refs = " << _orbrefs
						<< "  type = " << _orbtype
						<< "  intv = " << _orbintv
						<< "  nepo = " << _nepochs;
					_log->comment(2, "sp3", ltmp.str());
				}
				
			}
			else {
				if (_log) _log->comment(1, "sp3", " unknown record: " + tmp);
			}

		}
		else if (tmp.substr(0, 2) == "+ ") {
			ostringstream ltmp("reading satellites:");
			for (int i = 9; i < 60; i = i + 3) {
				string ssat;
				try {
					ssat = tmp.substr(i, 3);
				}
				catch (out_of_range) {
					break;
				}
				if (ssat.substr(0, 2) == "  ")
				{
					ssat = "G0" + ssat.substr(2, 1);
				}
				else if (ssat.substr(0, 1) == " ")
				{
					ssat = "G" + ssat.substr(1, 2);
				}
				_prn.push_back(ssat);
				ltmp << " " << tmp.substr(i, 3);
			}
			if (_log) _log->comment(2, "sp3", ltmp.str());

		}
		else if (tmp.substr(0, 2) == "++") {
			ostringstream ltmp("reading accuracies:");
			for (int i = 9; i < 60; i = i + 3) {
				int acc = 0;
				try {
					acc = str2int(tmp.substr(i, 3));
				}
				catch (out_of_range) {
					break;
				}
				_acc.push_back(acc);
				ltmp << " " << tmp.substr(i, 3);
			}
			if (_log) _log->comment(2, "sp3", ltmp.str());

		}
		else if (tmp.substr(0, 2) == "%c") {
				_timesys.push_back(tmp.substr(3, 2));
				_timesys.push_back(tmp.substr(6, 2));
				_timesys.push_back(tmp.substr(9, 3));
				_timesys.push_back(tmp.substr(13, 3));
				_timesys.push_back(tmp.substr(17, 4));
				_timesys.push_back(tmp.substr(22, 4));
				_timesys.push_back(tmp.substr(27, 4));
				_timesys.push_back(tmp.substr(32, 4));
				_timesys.push_back(tmp.substr(37, 5));
				_timesys.push_back(tmp.substr(43, 5));
				_timesys.push_back(tmp.substr(49, 5));
				_timesys.push_back(tmp.substr(55, 5));
				if (_log) _log->comment(2, "sp3", "reading satellite systems");

		}
		else if (tmp.substr(0, 2) == "%f") {
			_accbase.push_back(str2int(tmp.substr(3, 9)));
			_accbase.push_back(str2int(tmp.substr(14, 9)));
			if (_log) _log->comment(2, "sp3", "reading PV base");

		}
		else if (tmp.substr(0, 2) == "%i") {
			if (_log) _log->comment(2, "sp3", "additional info");

		}
		else if (tmp.substr(0, 2) == "/*") {
			if (_log) _log->comment(2, "sp3", "comments");

		}
		else if (tmp.substr(0, 2) == "* ") {  // END OF HEADER !!!
			if (_log) _log->comment(2, "sp3", "END OF HEADER");
			_mutex.unlock(); return -1;

		}
		else {
			if (_log) _log->comment(2, "sp3", "END OF HEADER");
			else       cerr << "warning: unknown header message :" << tmp << endl;
		}

		t_gcoder::_consume(tmpsize);
	}

	_mutex.unlock(); return consume;
}

   
/* ----------
 * SP3 body
 */
int t_sp3::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
	gtrace("t_sp3::decode_data");

	_mutex.lock();

	if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

	string tmp;
	//  int consume = 0;
	//  int recsize = 0;
	int tmpsize = 0;
	//  bool epoch_defined = false;

	//  while( ( tmpsize = t_gcoder::_getline( tmp, recsize ) ) >= 0 ){
	while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0) {

#ifdef DEBUG
		cout << " 0: " << tmp;
#endif

		// EPOCH record
		if (tmp.substr(0, 1) == "*" ||
			tmp.substr(0, 3) == "EOF"
			) {

			if (tmp.substr(0, 3) == "EOF") {
				if (_log) _log->comment(3, "sp3", "EOF found");
				t_gcoder::_consume(tmpsize);
				_mutex.unlock(); return tmpsize;
			}

			if (_nrecord > 0 && _nrecord != _nrecmax) {
				ostringstream ltmp;
				ltmp << "warning: not equal number of satellites " << _nrecord << " " << _nrecmax << "!";
				if (_log) _log->comment(0, "sp3", ltmp.str());
				else       cerr << ltmp.str();
			}

			char dummy;
			int yr, mn, dd, hr, mi;
			double sc;
			stringstream ss(tmp);
			ss >> dummy >> yr >> mn >> dd >> hr >> mi >> sc;

			if (ss.fail()) {
				printf("Warning: incorrect SP3 epoch record:%s \n\n", ss.str().c_str());
				t_gcoder::_consume(tmpsize);
				_mutex.unlock(); return -1;
			}

			int sod = hr * 3600 + mi * 60 + (int)sc;
			_lastepo.from_ymd(yr, mn, dd, sod);
			ostringstream ltmp;
			ltmp << "reading EPOCH [" << _nrecord << "] - " << _lastepo.str(" %Y-%m-%d %H:%M:%S");
			if (_log) _log->comment(3, "sp3", ltmp.str());
			_nrecord = -1;
		}

		// POSITION reccord
		if (tmp.substr(0, 1) == "P") { // and epoch_defined ){

			t_gtriple  xyz(0.0, 0.0, 0.0);
			t_gtriple dxyz(0.0, 0.0, 0.0);
			double t = 0.0, dt = 0.9;

			char sat[3 + 1]; sat[3] = '\0';
			char flg;
			double pos[4] = { 0.0, 0.0, 0.0, 0.0 };
			//      double var[4] = {0.0, 0.0, 0.0, 0.0};
			stringstream ss(tmp);

			ss >> noskipws >> flg >> sat[0] >> sat[1] >> sat[2]
				>> skipws >> pos[0] >> pos[1] >> pos[2] >> pos[3];
			string prn;
			prn = t_gsys::eval_sat(string(sat));
				

			for (int i = 0; i < 3; i++)
				if (pos[i] == 0.0)
					xyz[i] = UNDEFVAL_POS; // gephprec
				else  xyz[i] = pos[i] * 1000;  // km -> m

			if (pos[3] > 999999)
				t = UNDEFVAL_CLK;        // gephprec
			else  t = pos[3] / 1000000;      // us -> s

#ifdef DEBUG
			cout << " reading:" << _lastepo.str(" %Y-%m-%d %H:%M:%S ") << flg << " " << prn << " " << sat
				<< " " << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3] * 1000000 << endl;
#endif

			// fill single data record
			if (!_filter_gnss(prn)) {
				if (_log) _log->comment(4, "sp3", "skip " + prn);
			}
			else {
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end()) {
					if (it->second->id_type() == t_gdata::ALLPREC)
					{
						((t_gallprec*)it->second)->addpos(prn, _lastepo, xyz, t, dxyz, dt);	
						((t_gallprec*)it->second)->add_interval(prn,_orbintv);
						((t_gallprec*)it->second)->add_agency(_agency);

					}

					it++;
				}
			}
			
			if (ss.fail()) {
				if (_log) _log->comment(1, "sp3", "warning: incorrect SP3 data record: " + ss.str());
				else       cerr << "warning: incorrect SP3 data record: " << ss.str() << endl;
				t_gcoder::_consume(tmpsize);
				_mutex.unlock(); return -1;
			}
			cnt++;
		}

		// VELOCITY reccord
		if (tmp.substr(0, 1) == "V") { // and epoch_defined ){

			char sat[3 + 1]; sat[3] = '\0';
			char flg;
			double vel[4] = { 0.0, 0.0, 0.0, 0.0 };
			double var[4] = { 0.0, 0.0, 0.0, 0.0 };

			stringstream ss(tmp);
			ss >> noskipws >> flg >> sat[0] >> sat[1] >> sat[2]
				>> skipws >> vel[0] >> vel[1] >> vel[2] >> vel[3];

			string prn;
			prn = t_gsys::eval_sat(string(sat));

			// fill single data record
			map<string, t_gdata*>::iterator it = _data.begin();
			while (it != _data.end()) {
				if (!_filter_gnss(prn))
				{
					if (_log) _log->comment(4, "sp3", "skip " + prn);
					it++;
				}
				else
				{
					if (it->second->id_type() == t_gdata::ALLPREC)
					{
						((t_gallprec*)it->second)->addvel(prn, _lastepo, vel, var);

					}
					it++;
					//     }
				}
				
			}

			if (ss.fail()) {
				if (_log) _log->comment(1, "sp3", "warning: incorrect SP3 data record: " + ss.str());
				else       cerr << "warning: incorrect SP3 data record: " << ss.str() << endl;
				t_gcoder::_consume(tmpsize);
				_mutex.unlock(); return -1;
			}
		}

		_nrecord++;

		if (_nrecord > _nrecmax) _nrecmax = _nrecord;

		t_gcoder::_consume(tmpsize);
	}
	_mutex.unlock(); return tmpsize;
}

int t_sp3::encode_head(char* buff, int sz, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	// init ss position
	_mutex.lock();
	int gpsweek, sow, year, month, day, hour, minute, sec;

	double dsec;
	string data_used = "ORBIT";
	string frame = "IGS14";
	string ctype = "FIT";
	string agency = "SGG";
	string time_type = "GPST";
	vector<string> acc;
	int nsats;
	try
	{
		if (_ss_position == 0)
		{

			// initial val from data
			auto data_iter = _data.begin();
			for (data_iter; data_iter != _data.end(); data_iter++)
			{
				if (data_iter->second->id_type() != t_gdata::ALLPREC)
					continue;

				_start = ((t_gallprec*)data_iter->second)->get_beg();
				_lastepo = ((t_gallprec*)data_iter->second)->get_end();

				gpsweek = _start.gwk();
				sow = _start.sow();
				_start.ymd(year, month, day);
				_start.hms(hour, minute, sec);
				dsec = sec * 1.0;

				_sat = ((t_gallprec*)data_iter->second)->get_sats();

				_prn = ((t_gallprec*)data_iter->second)->get_sat3();
				_data_type = ((t_gallprec*)data_iter->second)->get_data_type();
				_sattype = ((t_gallprec*)data_iter->second)->get_sat_type();
				_orbintv = ((t_gallprec*)data_iter->second)->intv(_prn[0]);
				_nepochs = int((_lastepo - _start) / _orbintv) + 1;
				if (double_eq(_orbintv, 0.0))
				{
					_orbintv = ((t_gallprec*)data_iter->second)->intv();
				}
				agency = ((t_gallprec*)data_iter->second)->get_agency();
				_maxsats = _prn.size();
				nsats = _sat.size();

				for (int ii = 0; ii < _maxsats; ii++)
				{
					acc.push_back("  0");
				}
			}
			_ss << "#d" << _data_type
				<< setw(4) << right << year
				<< setw(3) << right << month
				<< setw(3) << right << day
				<< setw(3) << right << hour
				<< setw(3) << right << minute
				<< setw(12) << fixed << setprecision(8) << showpoint << right << dsec
				<< " " << setw(7) << right << _nepochs
				<< " " << setw(5) << right << data_used
				<< " " << setw(5) << right << frame
				<< " " << setw(3) << right << ctype;

			_ss << " " << setw(4) << right << time_type << endl;
			_ss << "## " << setw(4) << right << gpsweek
				<< setw(16) << fixed << setprecision(8) << showpoint << right << double(sow)
				<< setw(15) << fixed << setprecision(8) << showpoint << right << double(_orbintv)
				<< " " << setw(5) << right << _start.mjd()
				<< " " << setw(15) << fixed << setprecision(13) << showpoint << _start.sod() / 86400 << endl;

			_ss << "+" << "  " << setw(3) << right << nsats << "   ";

			for (int i = 0; i < _prn.size(); i++)
			{

				if (i % 17 == 0 && i != 0)
				{
					_ss << endl;
					_ss << "+" << "        ";
				}
				_ss << setw(3) << right << _prn[i];
			}


			for (int i = 0; i < acc.size(); i++)
			{
				if (i % 17 == 0)
				{
					_ss << endl;
					_ss << "++" << "       ";
				}
				_ss << setw(3) << right << acc[i];
			}

			_ss << endl;
			_ss << "%c" << " " << "M" << " " << " cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc" << endl;

			_ss << "/* generated by GREAT" << endl;


		}

		int size = _fill_buffer(buff, sz);

		_mutex.unlock();  return size;
	}

	catch (...)
	{
		if (_log)
		{
			_log->comment(0, "ERROR : t_sp3::encode_head throw exception");
		}
		else
		{
			cout << "ERROR : t_sp3::encode_head throw exception" << endl;
		}
		_mutex.unlock();
		return -1;
	}
}

int t_sp3::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_mutex.lock();
	try
	{
		if (_ss_position == 0)
		{
			auto data_iter = _data.begin();
			for (data_iter; data_iter != _data.end(); data_iter++)
			{
				if (data_iter->second->id_type() != t_gdata::ALLPREC)
					continue;

				t_gtime ref;
				
				int nsize = _sat.size();
				t_gtime epoch = _start;
				int year, month, day, hour, minute,sec;
				
				while (epoch <= _lastepo)
				{
					epoch.ymd(year, month, day);
					epoch.hms(hour, minute, sec);
					
					double xyz[3] = { 0.0,0.0,0.0 }, vel[3] = { 0.0,0.0,0.0 };
					double clk = 0.0;
					int obs_num = 0;

					for (int i = 0; i < nsize; i++)
					{
						
						if(!((t_gallprec*)data_iter->second)->get_pos_vel(_prn[i], epoch, xyz, vel, obs_num)) continue;
						if (double_eq(xyz[0]*xyz[1]*xyz[2],0.0))
						{
							continue;
						}
						if (i == 0)
						{
							_ss << "*" << setw(6) << right << year << setw(3) << right << month
								<< setw(3) << right << day << setw(3) << right << hour
								<< setw(3) << right << minute
								<< setw(12) << fixed << setprecision(8) << showpoint << right << double(sec) << endl;
						}
						_ss << "P" << _prn[i]
							<< setw(14) << fixed << setprecision(6) << showpoint << right << xyz[0]
							<< setw(14) << fixed << setprecision(6) << showpoint << right << xyz[1]
							<< setw(14) << fixed << setprecision(6) << showpoint << right << xyz[2]
							<< setw(14) << fixed << setprecision(6) << showpoint << right << clk
							<< setw(4) << right << obs_num << setw(4) << right << 0 << endl;
						if (_data_type == "V")
						{
							_ss << "V" << _prn[i]
								<< setw(14) << fixed << setprecision(6) << showpoint << right << vel[0]
								<< setw(14) << fixed << setprecision(6) << showpoint << right << vel[1]
								<< setw(14) << fixed << setprecision(6) << showpoint << right << vel[2]
								<< setw(14) << fixed << setprecision(6) << showpoint << right << clk
								<< setw(4) << right << 0 << endl;
						}
					}
					epoch = epoch + _orbintv;
				}
					_ss << "EOF" << endl;
					
			}
			
			//_get_delta_pos_vel		
		}
		int size = _fill_buffer(buff, sz);
		_mutex.unlock();
		return size;
	}
	catch (...)
	{
		if (_log)
		{
			_log->comment(0, "ERROR : t_sp3::encode_data throw exception");
		}
		else
		{
			cout << "ERROR : t_sp3::encode_data throw exception" << endl;
		}
		return -1;
	}

}


} // namespace
