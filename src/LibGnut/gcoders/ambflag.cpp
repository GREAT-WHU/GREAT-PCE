/**
*
* @verbatim
	History
	 -1.0 GREAT	    2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  ambflag.cpp
* @brief
*
* @author         GREAT, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/

#include"gcoders/ambflag.h"

using namespace std;
namespace great {

	t_ambflag::t_ambflag(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_ambflag::constructor");
	}

	/** @brief destructor. */
	t_ambflag::~t_ambflag()
	{
		gtrace("t_ambflag::destructor");
	}

	int t_ambflag::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		};
#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		string tmp;
		int tmpsize = 0;
		int consume = 0;
		string str;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				consume += tmpsize;
				istringstream istr(tmp);
				if (tmp.find("End of header") == string::npos)
				{
					if (tmp.find("Start time and interval :") != string::npos)
					{
						istr >> str >> str >> str >> str >> str >> _ambflag_head.beg_mjd >> _ambflag_head.beg_sod >> _ambflag_head.duration >> _ambflag_head.intv;
					}
					else if (tmp.find("Max ambc in one epoch   :") != string::npos)
					{
						istr >> str >> str >> str >> str >> str >> str >> _ambflag_head.max_amb_1epo;
					}
					else if (tmp.find("Old remvoed observations:") != string::npos)
					{
						istr >> str >> str >> str >> _ambflag_head.old_rm_obs;
					}
					else if (tmp.find("New removed observations:") != string::npos)
					{
						istr >> str >> str >> str >> _ambflag_head.new_rm_obs;
					}
					else if (tmp.find("Existed    ambiguities  :") != string::npos)
					{
						istr >> str >> str >> str >> _ambflag_head.exist_amb;
					}
					else if (tmp.find("New    ambiguities      :") != string::npos)
					{
						istr >> str >> str >> str >> _ambflag_head.new_amb;
					}
					else if (tmp.find("Available observations  :") != string::npos)
					{
						istr >> str >> str >> str >> _ambflag_head.avaiable_obs;
					}
					else
					{
						if(_log)
						{
							_log->comment(2, "t_ambflag::decode_head", "WARNING: unknown ambflag-head message :" + tmp);
						}
						else
						{
							cout << "WARNING: t_ambflag::decode_head ,unknown ambflag-head message :" << tmp << endl;
						}
						t_gcoder::_consume(tmpsize);
						_mutex.unlock();
						return -1;
					}
					t_gcoder::_consume(tmpsize);
				}
				else
				{
					//add head data
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()){
						if (it->second->id_type() == t_gdata::AMBFLAG || it->second->id_type() == t_gdata::AMBFLAG13)
						{
							string site_name;
							site_name = this->_fname.substr(_fname.rfind(".log") -12).substr(0, 4);
							((t_gallambflag*)it->second)->addAmbFlagHead(site_name, _ambflag_head);
						}
						it++;
					}
					t_gcoder::_consume(tmpsize);
					_mutex.unlock();
					return -1;
				}
			}
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_ambflag::decode_head", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambflag::decode_head , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
		_mutex.unlock();
		return consume;
	}

	int t_ambflag::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();

		if (t_gcoder::_add2buffer(buff, sz) == 0)
		{
			_mutex.unlock();
			return 0;
		};
#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		string tmp;
		int tmpsize = 0;
		ambflag_data data_tmp;
		string sat_name;
		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				istringstream istr(tmp);
				if (tmp.substr(0, 3) == "AMB" || tmp.substr(0, 3) == "IAM" || tmp.substr(0, 3) == "DEL" || tmp.substr(0, 3) == "BAD")
				{
					istr >> data_tmp.identify 
						 >> sat_name 
						 >> data_tmp.beg_epo 
						 >> data_tmp.end_epo
						 >> data_tmp.iflag 
						 >> data_tmp.C1 
						 >> data_tmp.C2 
						 >> data_tmp.reason;

					//fill data loop
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end())
					{
						if (it->second->id_type() == t_gdata::AMBFLAG || it->second->id_type() == t_gdata::AMBFLAG13)
						{
							string rec_name = this->_fname.substr(_fname.rfind(".log") - 12).substr(0, 4);
							((t_gallambflag*)it->second)->addAmbFlagData(rec_name, sat_name, data_tmp);
						}
						it++;
					}
				}
				else
				{
					if (_log)
					{
						_log->comment(2, "t_ambflag::decode_data", "WARNING: unknown ambflag-data message :" + tmp);
					}
					else
					{
						cout << "WARNING: unknown ambflag-data message :" << tmp << endl;
					}
					_mutex.unlock();
					return -1;
				}
				t_gcoder::_consume(tmpsize);
			}
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_ambflag::decode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambflag::decode_data , unknown mistake" << endl;
			}
			throw(-1);
		}
		_mutex.unlock();
		return tmpsize;
	}

	int t_ambflag::encode_head(char* buff, int sz, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				//get upd head info
				shared_ptr<ambflag_hd> head;
				map<string, t_gdata*>::iterator it = _data.begin();
				while (it != _data.end())
				{
					if (it->second->id_type() == t_gdata::AMBFLAG || it->second->id_type() == t_gdata::AMBFLAG13)
					{
						string site_name;
						site_name = this->_fname.substr(_fname.rfind(".log") - 12).substr(0, 4);
						head = dynamic_cast<t_gallambflag*>(it->second)->getOneAmbFlag(site_name).getAmbFlagHead();
					}
					it++;
				}

				_max_epo = floor(head->duration / head->intv);

				// fill head data
				_ss << fixed << "%Start time and interval :" << setw(6) << head->beg_mjd << setw(10) << setprecision(3) << head->beg_sod
					<< setw(9) << head->duration << setw(7) << setprecision(2) << head->intv << endl;
				_ss << "%Max ambc in one epoch   :" << setw(6) << head->max_amb_1epo << endl;
				_ss << "%Old remvoed observations:" << setw(12) << head->old_rm_obs << endl;
				_ss << "%New removed observations:" << setw(12) << head->new_rm_obs << endl;
				_ss << "%Existed    ambiguities  :" << setw(12) << head->exist_amb << endl;
				_ss << "%New    ambiguities      :" << setw(12) << head->new_amb << endl;
				_ss << "%Available observations  :" << setw(12) << head->avaiable_obs << endl;
				_ss << "%End of header" << endl;
			}
			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_ambflag::encode_head", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambflag::encode_head , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
	}

	int t_ambflag::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		try
		{
			if (_ss_position == 0)
			{
				//get data from _data
				t_map_sat_ambflag data_tmp;
				auto it = _data.begin();
				for (it = _data.begin(); it != _data.end(); ++it)
				{
					if (it->second->id_type() == t_gdata::AMBFLAG || it->second->id_type() == t_gdata::AMBFLAG13)
					{
						string rec_name = this->_fname.substr(_fname.rfind(".log") - 12).substr(0, 4);
						data_tmp = dynamic_cast<t_gallambflag*>(it->second)->getOneAmbFlag(rec_name).getAmbFlagData();
					}
				}

				for (auto itsat = data_tmp.begin(); itsat != data_tmp.end(); ++itsat)
				{
					for (auto itvet = itsat->second.begin(); itvet != itsat->second.end(); itvet++)
					{
						for (int i = 0; i < _max_epo; i++)
						{
							if ((*itvet)->beg_epo == i + 1)
							{
								_ss << fixed << setw(3) << (*itvet)->identify << setw(4) << itsat->first << setw(7) << (*itvet)->beg_epo
									<< setw(7) << (*itvet)->end_epo << setw(4) << (*itvet)->iflag << setw(20)
									<< setprecision(3) << (*itvet)->C1 << setw(20) << setprecision(3) << (*itvet)->C2
									<< "  " << (*itvet)->reason << endl;
							}
						}
					}
				}
			}
			int size = _fill_buffer(buff, sz);
			_mutex.unlock();
			return size;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(2, "t_ambflag::encode_data", "ERROR: unknown mistake");
			}
			else
			{
				cout << "ERROR: t_ambflag::encode_data , unknown mistake" << endl;
			}
			return -1;
			throw(-1);
		}
		
	}
}//namespace