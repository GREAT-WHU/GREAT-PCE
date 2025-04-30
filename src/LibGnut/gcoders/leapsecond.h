/**
 * The leap_second file store the date in which happened leapsecond
 * Following is an example for the leap_second file.
 *
 * @verbatim
		+leap sec		! begin of the data
		41317  10                 ! 31 DEC 1971 LEAP SEC INCREMENT (at the end of the day)	! date happened leapping, second
		41499  11                 ! 30 JUN 1972 LEAP SEC INCREMENT
		41683  12                 ! 31 DEC 1972 LEAP SEC INCREMENT
		42048  13                 ! 31 DEC 1973 LEAP SEC INCREMENT
		42413  14                 ! 31 DEC 1974 LEAP SEC INCREMENT
		42778  15                 ! 31 DEC 1975 LEAP SEC INCREMENT
		43114  16                 ! 31 DEC 1976 LEAP SEC INCREMENT
		43509  17                 ! 31 DEC 1977 LEAP SEC INCREMENT
		43874  18                 ! 31 DEC 1978 LEAP SEC INCREMENT
		44239  19                 ! 31 DEC 1979 LEAP SEC INCREMENT
		44786  20                 ! 30 JUN 1981 LEAP SEC INCREMENT
		45151  21                 ! 30 JUN 1982 LEAP SEC INCREMENT
		45516  22                 ! 30 JUN 1983 LEAP SEC INCREMENT
		46247  23                 ! 30 JUN 1985 LEAP SEC INCREMENT
		47161  24                 ! 31 DEC 1987 LEAP SEC INCREMENT
		47892  25                 ! 31 DEC 1989 LEAP SEC INCREMENT
		48257  26                 ! 31 DEC 1990 LEAP SEC INCREMENT
		48804  27                 ! 30 JUN 1992 LEAP SEC INCREMENT
		49169  28                 ! 30 JUN 1993 LEAP SEC INCREMENT
		49534  29                 ! 30 JUN 1994 LEAP SEC INCREMENT
		50083  30                 ! 31 DEC 1995 LEAP SEC INCREMENT
		50630  31                 ! 30 JUN 1997 LEAP SEC INCREMENT
		51179  32                 ! 31 DEC 1998 LEAP SEC INCREMENT
		53736  33                 ! 31 DEC 2005 LEAP SEC INCREMENT
		54832  34                 ! 31 DEC 2008 LEAP SEC INCREMENT
		56109  35                 ! 30 JUN 2012 LEAP SEC INCREMENT
		57205  36                 ! 30 JUN 2015 LEAP SEC INCREMENT
		57752  37                 ! 31 DEC 2016 LEAP SEC INCREMENT (predicted)
		-leap sec		! end of the data
   @endverbatim
 *
 * @verbatim
	 History
	  -1.0	Bo Wong	2019-01-04 creat the file.
	  -1.1  Jie Li	2019-03-30 adding doxygen style code remarks
   @endverbatim
 * Copyright (c) 2018, Wuhan University. All rights reserved.
 *
 * @file		leapsecond.h
 * @brief		The base class used to decode leapsecond file information.
 *
 * @author      Bo Wong, Wuhan University
 * @version		1.0.0
 * @date		2019-03-30
 */

#ifndef LEAPSECOND_H
#define LEAPSECOND_H

#include "gexport/ExportLibGnut.h"
#include "gcoders/gcoder.h"
#include "gdata/gleapsecond.h"

using namespace gnut;

namespace great
{

	/**
	*@brief	   Class for decoding the leapsecond data
	*
	* The document contains date and value of leap seond
	* The gcoder t_leapsecond corresponding to the gdata leapsecond
	*/
	class LibGnut_LIBRARY_EXPORT t_leapsecond : public t_gcoder {

	public:
		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_leapsecond(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_leapsecond();

		/**
		* @brief decode header of leap_second file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of header decoding
			@retval <0  finish reading
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief decode data body of leap_second file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return
			@retval >=0 consume size of body decoding
			@retval <0  finish reading
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

	protected:
		int _mjd;    ///< the mjd of leap seconds happen
		int _leap;   ///< leap seconds /s
	};
} //namespace

#endif