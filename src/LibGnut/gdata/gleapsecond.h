/**
 *
 * @verbatim
	 History
	  -1.0  Bo Wong  2019-01-04 Creat the file.
	  -1.1  Jie Li   2019-03-30 Code comments are added in Doxygen format.
   @endverbatim
 * Copyright (c) 2018, Wuhan University.All rights reserved.
 *
 * @file		gleapsecond.h
 * @brief		The class for storaging leapsecond data.
 *
 * @author      Bo Wong, Wuhan University
 * @version		1.0.0
 * @date		2019-03-30
 *
 */

#ifndef GLEAPSECOND_H
#define GLEAPSECOND_H

#include "gcoders/leapsecond.h"
#include "gdata/gdata.h"

#include "gexport/ExportLibGnut.h"

using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for save leapsecond data
	*
	* This class save leapsecond data and provide the interface to get data.
	*/
	class LibGnut_LIBRARY_EXPORT t_gleapsecond : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gleapsecond();

		/** @brief default destructor. */
		virtual ~t_gleapsecond();

		/**
		* @brief add leapsecond data to map.
		* @param[in]   mjd   mjd of the leapsecon data
		* @param[in]   leap  leapsecon data
		*/
		void add_data(int mjd, int leap);

		/** @brief whether the leapsecond data is empty.
		* @return  bool
		*	@retval   true     leapsecond data is empty
		*   @retval   false    leapsecond data is existent
		*/
		bool is_empty();

		/**
		* @brief get hole leapsecond data
		* @return leapsecond data include time and seconds
		*/
		map<int, int> &get_leap() { return _leapseconds; };

	protected:
		map<int, int> _leapseconds;		///< leapsecond data(mjd,value)

	private:

	};

}
#endif // !GALLPLANETEPH_H
