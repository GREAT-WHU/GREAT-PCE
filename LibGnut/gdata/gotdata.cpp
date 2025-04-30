/**
 * @file         gotdata.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        The class for storaging ocean_tide data.
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#include"gotdata.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gocean_tide::t_gocean_tide() :t_gdata()
	{
		id_type(t_gdata::OCEANTIDE);
	}


	/** @brief whether the ocean tide data is empty.
	* @return  bool
	*	@retval   true   ocean tide data is empty
	*   @retval   false  ocean tide data is existent
	*/
	bool t_gocean_tide::is_empty()
	{
		if (_oct.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	/**
	* @brief get all oceantide coefficient.
	* @return   all oceantide coefficient
	*/
	vector<t_octrec>& t_gocean_tide::oct()
	{ 
		return _oct; 
	}

	/**
	* @brief get all load coefficient.
	* @return   all load coefficient
	*/
	vector<double>& t_gocean_tide::load() 
	{ 
		return _load; 
	}

	/**
	* @brief add one oceanide data record to _oct.
	* @param[in]   m		index m
	* @param[in]   n		index n
	* @param[in]   index	Doodson's index
	* @param[in]   coeff	coefficient(c_prog,c_retg,s_prog,s_cretg)
	*/
	void t_gocean_tide::set_oct(const int& m, const int& n, int index[], double coeff[])
	{
    #ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
    #endif 
		_gmutex.lock();
		int i, j;
		t_octrec oct;
		oct.m = m;
		oct.n = n;
		for ( i = 0; i < 6; i++)
		{
			oct.index[i] = index[i];
		}
		for ( j = 0; j < 4; j++)
		{
			oct.coeff[j] = coeff[j];
		}
		_oct.push_back(oct);
		_gmutex.unlock();
		return;
	}

	/**
	* @brief add a load coefficient value to vector.
	* @param[in]   load		load coefficient value
	*/
	void t_gocean_tide::set_load(const double& load)
	{
    #ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
    #endif 
		_gmutex.lock();
		_load.push_back(load);
		_gmutex.unlock();
		return;
	}

	/** @brief constructor. */
	t_octrec::t_octrec()
	{

	}

}

