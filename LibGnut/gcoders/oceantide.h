/**
 * The oceantide file is used for calculating the effect of ocean tide 
	to gravitational coefficient
 * Following is an example for the ocentide file.
 *
 * @verbatim
		+coefficient			! begin symbol of the coefficient
		=================================================================================================================
		#  n   m   ---Doodson's Index---  ------C_PROG------  ------C_RETG------  ------S_PROG------  ------S_RETG------
		   2   0   0   0   0   0   1   0  0.967395410000d+00  0.000000000000d+00  0.000000000000d+00  0.000000000000d+00
		   3   0   0   0   0   0   1   0  0.160936920000d-01  0.000000000000d+00  0.000000000000d+00  0.000000000000d+00
		   4   0   0   0   0   0   1   0 -0.163259180000d+00  0.000000000000d+00  0.000000000000d+00  0.000000000000d+00
		   5   0   0   0   0   0   1   0  0.320794360000d+00  0.000000000000d+00  0.000000000000d+00  0.000000000000d+00
					...							...							...					...
		-coefficient			! end symbol of the coefficient 
		+load coefficient		! begin symbol of load coefficient
		========================
		order		value
		========================
		1   .00000000000000D+00
		2  -.30750000000000D+00
				...
		-load coefficient		! end symbol of load coefficient
   @endverbatim
 *
 * @verbatim
	 History
	  -1.0	Hongmin Zhang  2019-01-06 creat the file.
	  -1.1  Jie Li		   2019-03-30 adding doxygen style code remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		oceantide.h
 * @brief		The base class used to decode ocean_tide file information.
 *
 * @author      Hongmin Zhang, Wuhan University
 * @version		1.0.0
 * @date		2019-03-30
 */

#ifndef OCEANTIDE_H
#define OCEANTIDE_H

#include "gexport/ExportLibGnut.h"
#include "gcoders/gcoder.h"
#include "gdata/gdata.h"

#define MAX_CONSTITS     2000
#define MAX_OCEAN_DEGREE 35

using namespace gnut;

namespace great
{
	/**
	*@brief	   Class for decoding the oceantide data
	*
	* The gcoder t_oceantide corresponding to the gdata otdata
	*/
	class LibGnut_LIBRARY_EXPORT t_oceantide : public t_gcoder
	{
	public:
		/**
		 * @brief default constructor.
		 *
		 * @param[in]  s        setbase control
		 * @param[in]  version  version of the gcoder
		 * @param[in]  sz       size of the buffer
		 */
		t_oceantide(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_oceantide() {};

		/**
		* @brief decode data body of oceantide file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return 
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

		/**
		* @brief decode header of oceantide file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return 
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg);

	protected:
		int     n;            ///< the index of spheric harmoics
		int     m;            ///< the index of spheric harmoics
		int     index[6];     ///< the index of Doodson
		double  coeff[4];     ///< coefficients of the ocean tide(C_PROG, C_RETG, S_PROG, S_RETG)
		double  _load;		  ///< load coefficient
		int     _tab = 0;	  ///< data type symbol for decoding

	}; //namespace
}
#endif