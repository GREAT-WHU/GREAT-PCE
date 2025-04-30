/**
 * @file         poleut1.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        The base class used to decode poleut1 file information.
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include <string>
#include <algorithm>

#include"gcoders/poleut1.h"
#include "gdata/gpoleut1.h"
#include "gall/gallobj.h"

using namespace std;
namespace great {

	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_poleut1::t_poleut1(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz), _begtime(0), _endtime(0), _interval(0), _parnum(0)
	{
		gtrace("t_dvpteph405::constructor");
	}

	/**
	* @brief decode header of poleut1 file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_poleut1::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) {
			_mutex.unlock(); 
			return 0; 
		};

 
		string tmp;
		int consume = 0;
		int tmpsize = 0;
		string timetype;
		string time;
		string param;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp)) >= 0) {

				consume += tmpsize;

				// first information
				if (tmp.substr(0, 2) == "% ")
				{
					istringstream istr(tmp.substr(2));
					if (tmp.substr(2, 1) == "U")
					{
						istr >> timetype >> timetype >> timetype >> _timetype;
					}
					else if (tmp.substr(2, 1) == "S")
					{
						istr >> time >> time >> _begtime >> _endtime >> _interval;
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::ALLPOLEUT1)
                            {
								istr >> time >> time >> _begtime >> _endtime >> _interval;
                                //((t_gpoleut1*)it->second)->setBegEndTime(_begtime, _endtime);
								dynamic_cast<t_gpoleut1*>(it->second)->setBegEndTime(_begtime, _endtime);
                            }
							it++;
						}
					}
					else if (tmp.substr(2, 1) == "N")
					{
						istr >> param >> param >> param >> param >> _parnum;
						string unitstr;
						while (istr >> unitstr)
						{
							replace(unitstr.begin(), unitstr.end(), 'D', 'E');
							_parunit.push_back(atof(unitstr.c_str()));
						}
					}
				}
				else if (tmp.substr(0, 2) == "+p")
				{
					istringstream istr2(tmp.substr(2));
				}
			else if (tmp.substr(0, 2) == "%%")
				{
					if (tmp.substr(3, 1) == "M")
					{
						istringstream istr3(tmp.substr(6));
						string name;
						int i = 0;
						while (i < _parnum)
						{
							istr3 >> name;
							_parname.push_back(name);
							i++;
						}
					}
				}
				else {
					if (_log) _log->comment(2, "poleut1", "END OF HEADER");
					else       cerr << "warning: unknown header message :" << tmp << endl;
					_mutex.unlock(); return -1;
				}

				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock(); return consume;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_poleut1::decode_head throw exception");
			}
			else
			{
				cout << "ERROR : t_poleut1::decode_head throw exception" << endl;
			}
			return -1;
		}
	}

	/**
	* @brief decode data body of poleut1 file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_poleut1::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0){ 
			_mutex.unlock(); 
			return 0; 
		};

		string tmp;
		int tmpsize = 0;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				if (tmp.substr(0, 1) == " ")
				{
					double t;
					istringstream istr(tmp.substr(0, 65));
					vector<double> par;
					istr >> t;
					double pa;
					while (istr >> pa)
					{
						par.push_back(pa);
					}
					int mjd = floor(t);
					int sec = (int)((t - mjd) * 86400 / 3600.0) * 3600;
					t_gtime T;
					T.from_mjd(mjd, sec, 0.0);
					map<string, double> data;
					for (int i = 0; i < _parname.size(); i++)
					{
						data[_parname[i]] = par[i] * _parunit[i];
					}
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (T.mjd() > _endtime)
							break;
						if (it->second->id_type() == t_gdata::ALLPOLEUT1)
							//((t_gpoleut1*)it->second)->setEopData(T, data, _timetype, _interval);
							dynamic_cast<t_gpoleut1*>(it->second)->setEopData(T, data, _timetype, _interval);
						it++;
					}
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock(); return tmpsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_poleut1::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_poleut1::decode_data throw exception" << endl;
			}
			return -1;
		}
	}

	int t_poleut1::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
		try
		{
			if (_ss_position == 0)
			{
				for (auto itt = _data.begin(); itt != _data.end(); itt++)
				{

					if (itt->second->id_type() == t_gdata::ALLPOLEUT1)
					{
						map<t_gtime, map<string, double> >* mapEOP = ((t_gpoleut1*)itt->second)->getPoleUt1DataMap();
						auto iter_beg = mapEOP->begin();
						auto iter_end = (--mapEOP->end());
						_ss << "%% eopupd              yuanyongqiang       25-Jun-17" << endl;
						_ss << "%% Bulletin A file : finals2000A.data" << endl;
						_ss << "%% Input EOP  file :" << endl;
						_ss << "%%" << endl;
						_ss << "+pole&ut1" << endl;
						_ss << "% UT1 type = UT1R" << endl;
						_ss << "% Start&End%Interval =" << setw(9) << right << iter_beg->first.mjd() - 1 << setw(8) << right << iter_end->first.mjd() - 1 << setw(7) << right << fixed << setprecision(2) << 1.0 << endl;
						_ss << "% Num. of Vars&Units =     5  0.1D+01  0.1D+01  0.1D+01  0.1D-02  0.1D-02" << endl;
						_ss << "% Format = (f9.2,1x,2f10.6,f15.7,2f10.3,3(1x,a1))" << endl;
						_ss << "%% MJD        XPOLE     YPOLE      UT1-TAI        DPSI     DEPSI    PRED_ID" << endl;
						
						for (auto iterAll = mapEOP->begin(); iterAll != mapEOP->end(); iterAll++)
						{

							_ss << setw(9) << right << fixed << setprecision(2) << iterAll->first.mjd()-1.0;
							_ss << setw(11) << right << fixed << setprecision(6) << iterAll->second["XPOLE"];
							_ss << setw(10) << right << fixed << setprecision(6) << iterAll->second["YPOLE"];
							_ss << setw(15) << right << fixed << setprecision(7) << iterAll->second["UT1"];
							_ss << setw(10) << right << fixed << setprecision(3) << 0.0;
							_ss << setw(10) << right << fixed << setprecision(3) << 0.0;
							_ss << " I I I" << endl;
						}
					}
				}
			}


		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_result::encode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_result::encode_data throw exception" << endl;
			}
			return -1;
		}

		int size = _fill_buffer(buff, sz);

		_mutex.unlock();  return size;
	}

	
}
