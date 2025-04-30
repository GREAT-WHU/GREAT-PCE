/**
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang	2019-01-05	creat the file.
	  -1.1  Jie Li			2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		oceantide.cpp
 * @brief		The base class used to decode oceantide file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-01-05
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <algorithm>

#include "newmat/newmat.h"
#include"gcoders/oceantide.h"
#include"gdata/gotdata.h"

using namespace std;

namespace great {

	/**
	* @brief constructor.
	* @param[in]  s        setbase control
	* @param[in]  version  version of the gcoder
	* @param[in]  sz       size of the buffer
	*/
	t_oceantide::t_oceantide(t_gsetbase* s, string version, int sz)
		:t_gcoder(s, version, sz)
	{

	}

	/**
	* @brief decode header of oceantide file
	* @param[in]  buff        buffer of the data
	* @param[in]  sz          buffer size of the data
	* @param[in]  errmsg      error message of the data decoding
	* @return consume size of header decoding
	*/
	int t_oceantide::decode_head(char* buff, int sz, vector<string>& errmsg)
	{
		if (t_gcoder::_add2buffer(buff, sz) == 0) { return 0; };
		return -1;
	}
	/**
		* @brief decode data body of oceantide file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size for data body decoding
		*/
	int t_oceantide::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_mutex.lock();
		if (t_gcoder::_add2buffer(buff, sz) == 0) { _mutex.unlock(); return 0; };

#ifdef DEBUG
		cout << " BUFFER : \n" << _buffer << "\n size = " << sz << " END OF BUFFER \n\n"; cout.flush();
#endif
		//t_oceantide oct[IMAXCONSTITS];
		string line;
		string _line;
		int consume = 0;
		int tmpsize = 0;
		int sitsize = 0;
		int ii, jj, num;
		bool complete = false;
		int nconstits;
		int ndegree;
		Matrix dindex;
		double coeff[5];
		t_gocean_tide t_oceantd;
		dindex.ReSize(1, 6);


		tmpsize = t_gcoder::_getline(line);
		sitsize += tmpsize;
		consume += tmpsize;
		complete = false;

		//find the starting place
		/*while (line.find("+coefficient") == string::npos ) {
			t_gcoder::_consume(tmpsize);
			tmpsize = t_gcoder::_getline(line, sitsize);
			sitsize += tmpsize;
			consume += tmpsize;
			complete = false;

			}*/

		nconstits = 0;
		ndegree = 0;
		try
		{
			//find the ending place
			while ((tmpsize = t_gcoder::_getline(line)) >= 0)
			{
				sitsize += tmpsize;
				consume += tmpsize;

				if (line.find("-coefficient") != string::npos)
				{
					_tab = 1;
				}

				if (line.substr(0, 1) == " "&&_tab == 0)
				{
					nconstits += 1;
					if (nconstits > MAX_CONSTITS)
					{
						if (_log) _log->comment(2, "oceantide", "To many constits.");
					}
					istringstream istr(line);

					istr >> n >> m;

					//the index of Doodson
					for (int i = 1; i <= 6; i++)
					{
						istr >> dindex(1, i);
					}

					//the index of spheric harmoics
					for (int i = 1; i <= 2; i++)
					{
						istr >> _line;
						replace(_line.begin(), _line.end(), 'd', 'e');
						coeff[i] = atof(_line.c_str());
					}
					string line1 = line.substr(73, 19);
					replace(line1.begin(), line1.end(), 'd', 'e');
					string line2 = line.substr(93, 19);
					replace(line2.begin(), line2.end(), 'd', 'e');
					istringstream istr1(line1);
					istr1 >> coeff[3];
					istringstream istr2(line2);
					istr2 >> coeff[4];

					for (ii = 0; ii < 6; ii++)
					{
						index[ii] = dindex(1, ii + 1);
					}
					//coefficients of the ocean tide(C_PROG, C_RETG, S_PROG, S_RETG)
					for (jj = 0; jj < 4; jj++)
					{
						coeff[jj] = coeff[jj + 1];
					}
					if (n > ndegree)
						ndegree = n;
					if (m > ndegree)
						ndegree = m;
					if (ndegree > MAX_OCEAN_DEGREE)
					{
						if (_log) _log->comment(2, "oceantide", " degree/order larger than maximum");
					}
					else
					{
						complete = true;
					}

					if (complete = true)
					{
						//t_oceantd.setoct(_m, _n, _dindex, _coeff);
						map<string, t_gdata*>::iterator it = _data.begin();
						while (it != _data.end())
						{
							if (it->second->id_type() == t_gdata::OCEANTIDE)
							{
								((t_gocean_tide*)it->second)->set_oct(m, n, index, coeff);
								//the number of parameters to be eatimated for one satellite
								//((t_gotdata*)it->second)->set(satname, tab_data);
							}

							it++;
						}
						complete = false;
					}

				}

				if (line.substr(0, 1) == " "&&_tab == 1)
				{
					if (line.substr(0, 1) == " ")
					{
						string line3;
						istringstream istr(line);
						istr >> num;
						istr >> line3;
						replace(line3.begin(), line3.end(), 'D', 'e');
						_load = atof(line3.c_str());
						if (!istr.fail())
						{
							//
							//t_oceantd.setload(_load);
							map<string, t_gdata*>::iterator it = _data.begin();
							while (it != _data.end())
							{
								if (it->second->id_type() == t_gdata::OCEANTIDE)
								{
									((t_gocean_tide*)it->second)->set_load(_load);
								}

								it++;
							}
						}
					}

				}


				t_gcoder::_consume(tmpsize);
			}


			/*while (line.find("+load coefficient") == string::npos)
			{
			t_gcoder::_consume(tmpsize);
			tmpsize = t_gcoder::_getline(line, sitsize);
			sitsize += tmpsize;
			consume += tmpsize;
			if (line == "")
			{
			return tmpsize;
			}
			}*/

			_mutex.unlock(); return sitsize;
		}
		catch (...)
		{
			if (_log)
			{
				_log->comment(0, "ERROR : t_oceantide::decode_data throw exception");
			}
			else
			{
				cout << "ERROR : t_oceantide::decode_data throw exception" << endl;
			}
			return -1;
		}

	}

}

