/**
 *
 * @verbatim
	 History
	  -1.0	Bo Wong			2019-01-04	creat the file.
	  -1.1  Jie Li			2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		leapsecond.cpp
 * @brief		The base class used to decode leapsecond file information.
 *
 * @author      Bo Wong, Wuhan University
 * @version		1.0.0
 * @date		2019-01-04
 */

#include"gcoders/leapsecond.h"

using namespace std;
namespace great {

	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_leapsecond::t_leapsecond(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{
		gtrace("t_leapsecond::constructor");
	}

	/** @brief destructor. */
	t_leapsecond::~t_leapsecond()
	{

	}

	/**
	* @brief decode header of leap_second file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_leapsecond::decode_head(char* buff, int sz, vector<string>& errmsg)
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
		_mutex.unlock();
		return -1;
	}

	/**
	* @brief decode data body of leap_second file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size for data body decoding
	*/
	int t_leapsecond::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
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

		try
		{
			while ((tmpsize = t_gcoder::_getline(tmp, 0)) >= 0)
			{
				if (tmp.substr(0, 1) == "+")
				{
					t_gcoder::_consume(tmpsize);
					continue;
				}
				else if (tmp.substr(0, 1) == " ")
				{
					istringstream istr(tmp);
					istr >> _mjd >> _leap;
					//fill data loop
					map<string, t_gdata*>::iterator it = _data.begin();
					while (it != _data.end()) {
						if (it->second->id_type() == t_gdata::LEAPSECOND)
							((t_gleapsecond*)it->second)->add_data(_mjd, _leap);
						it++;
					}
				}
				else if (tmp.substr(0, 1) == "-")
				{
					_mutex.unlock();
					return tmpsize;
				}
				t_gcoder::_consume(tmpsize);
			}
			_mutex.unlock();
			return tmpsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_leapsecond::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_leapsecond::decode_data throw exception" << endl;
			}
			_mutex.unlock();
			return -1;
		}
	}

}//namespace