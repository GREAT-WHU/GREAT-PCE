/**
 *
 * @verbatim
	 History
	  -1.0	Bo Wong		2019-01-04	creat the file.
	  -1.1  Jie Li		2019-04-01	add code in doxygen style.
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		gleapsecond.cpp
 * @brief		The base class used to storage the leapsecond data.
 *
 * @author      Bo Wong, Wuhan University
 * @version		1.0.0
 * @date		2019-01-04
 */


#include"gdata/gleapsecond.h"

using namespace std;
namespace great {

	/** @brief constructor. */
	t_gleapsecond::t_gleapsecond(){
		gtrace("t_gleapsecond::constructor");
		id_type(t_gdata::LEAPSECOND);
	}

	/** @brief destructor. */
	t_gleapsecond::~t_gleapsecond(){}

	/**
	* @brief add leapsecond data to map.
	* @param[in]   mjd   mjd of the leapsecon data
	* @param[in]   leap  leapsecon data
	*/
	void t_gleapsecond::add_data(int mjd, int leap){
		_leapseconds[mjd] = leap;
	}

	/** @brief whether the leapsecond data is empty.
	* @return  bool
	*	@retval   true     leapsecond data is empty
	*   @retval   false    leapsecond data is existent
	*/
	bool t_gleapsecond::is_empty()
	{
		if (_leapseconds.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}//namespace