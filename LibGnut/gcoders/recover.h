/**
*
* @verbatim
History
-1.0 Hongjie Zheng  2019-09-25  creat the file.
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		recover.h
* @brief	The base class used to decode and encode recover file information.
*
* @author   Hongjie Zheng, Wuhan University
* @version	1.0.0
* @date		2019-09-25
*
*/

#ifndef RECOVER_H
#define RECOVER_H


#include <gcoders/gcoder.h>
#include <gdata/grecoverdata.h>
#include <gall/gallrecover.h>
using namespace gnut;
namespace great
{

	class LibGnut_LIBRARY_EXPORT t_resfile :public t_gcoder
	{
	public:
		const static string END_OF_HEADER;
		const static string TIME_HEADER;
		const static string SIGMA_HEADER;
		const static string SAT_HEADER;
		const static string SITE_HEADER;
		const static string OBS_DATA;
		const static string PAR_DATA;

	public:

		t_resfile(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

		virtual ~t_resfile();

		/**
		* @brief decode header of resfile 
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief decode data of  resfile
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return consume size of header decoding
		*/
		virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

		/**
		* @brief encode header of  resfile
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of header encoding
		*/
		virtual  int encode_head(char* buff, int sz, vector<string>& errmsg)override;

		/**
		* @brief encode data of  resfile
		* @param[in]  buff        buffer of the data
		* @param[in]  sz          buffer size of the data
		* @param[in]  errmsg      error message of the data decoding
		* @return  size of data body encoding
		*/
		virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)override;

	protected:

		void _add_data(string id, t_gdata* data) override;
		
		t_gallrecover* _recover_data;
	};


	LibGnut_LIBRARY_EXPORT t_grecover_par strline2recover_par(string line);
	LibGnut_LIBRARY_EXPORT t_grecover_equation strline2recover_equation(string line);

}


#endif // !RECOVER_H
