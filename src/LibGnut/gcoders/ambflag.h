/**
*
* @verbatim
	History
	 -1.0 GREAT	    2019-01-04 creat the file.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			  ambflag.h
* @brief
*
* @author         GREAT, Wuhan University
* @version		  1.0.0
* @date			  2019-01-04
*
*/

#ifndef AMBFLAG_H
#define AMBFLAG_H

#include "gexport/ExportLibGnut.h"
#include "gcoders/gcoder.h"
#include "gall/gallambflag.h"

using namespace std;
using namespace gnut;

namespace great
{

	/**
	*@brief	   Class for decode/encode ambflag file
	*/
	class LibGnut_LIBRARY_EXPORT t_ambflag : public t_gcoder
	{

	public:

		/**
		* @brief constructor.
		* @param[in]  s        setbase control
		* @param[in]  version  version of the gcoder
		* @param[in]  sz       size of the buffer
		*/
		t_ambflag(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		/** @brief default destructor. */
		virtual ~t_ambflag();

		/**
		* @brief decode header of ambflag file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief decode data body of ambflag file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size for data body decoding
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

		/**
		* @brief encode header of ambflag file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of header encoding
		*/
		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief encode data body of ambflag file
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size for data body encoding
		*/
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

	protected:
		ambflag_hd    _ambflag_head;  ///< map for storaging ambflag data header
		ambflag_data  _ambflag_data;  ///< map for storaging ambflag data body
		int           _max_epo;       ///< max epoch begin
	private:

	};

} // namespace

#endif
