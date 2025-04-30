/**
*
* @verbatim
	History

	@endverbatim
*
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

*
* @file        gsetout.h
* @brief       implements output setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/


#ifndef GSETOUT_H
#define GSETOUT_H

#define XMLKEY_OUT "outputs"        ///< The defination of outputs module in xml file

#define DEFAULT_FILE_VER ""         ///< default version for file format
#define DEFAULT_FILE_UPD     0      ///< default auto update [min] for file saving
#define DEFAULT_FILE_LEN     0      ///< default auto length [min] for file saving
#define DEFAULT_FILE_SMP     0      ///< default auto sample [sec] for file saving
#define DEFAULT_FILE_OFF     0      ///< default file offset [min] for file saving
#define DEFAULT_FILE_SYS  "UTC"     ///< default file system       for file saving

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut
{

	///< the order is important here!
	enum OFMT
	{
		XXX_OUT,
		LOG_OUT,
		PPP_OUT,
		FLT_OUT,

		RECCLK_OUT,      ///< rec clk out
		SATCLK_OUT,      ///< sat clk out

		CLK_OUT, 
		EPODIR_OUT
	};

	/// The class for settings of output
	class LibGnut_LIBRARY_EXPORT t_gsetout : public virtual t_gsetbase
	{
	public:
		/// constructor
		t_gsetout();
		/// destructor
		~t_gsetout();

		/**
		 * @brief change from string to OFMT
		 * @param[in] s file format
		 * @return OFMT : file format
		 */
		static OFMT   str2ofmt(const string& s);

		/**
		 * @brief change from OFMT to string
		 * @param[in] f file format
		 * @return string : file format
		 */
		static string ofmt2str(const OFMT&   f);

		/// settings check
		void check();         
		/// settings help
		void help();                               

		// attributes
		/**
		 * @brief  get verbosity attribute
		 * @return int : verbosity attribute
		 */
		int verb();   

		/**
		 * @brief  get append request
		 * @return bool : append request
		 */
		bool append(); 

		int sp3_obslimit();

		/// update glog (mask,verb)
		virtual void _upd_glog();                  

		// elements
		/**
		 * @brief  get format output size
		 * @param[in] fmt file format
		 * @return int : format output size
		 */
		int output_size(const string& fmt);    

		/**
		 * @brief  get string outputs
		 * @param[in] fmt file format
		 * @return string : string outputs
		 */
		string outputs(const string& fmt);


		/**
		 * @brief  get formats
		 * @return set<string> : all the outputs
		 */
		set<string> oformats();                     

		/**
		 * @brief  get string output version
		 * @param[in] fmt file format
		 * @return string : string output version
		 */
		string version(const string& fmt);   
		
		/**
		 * @brief  get length [min] for output file content
		 * @param[in] fmt file format
		 * @return int : length [min] for output file content
		 */
		int   out_length(const string& fmt);     
		
		/**
		 * @brief  get sample [sec] for output file data
		 * @param[in] fmt file format
		 * @return float : sample [sec] for output file data
		 */
		float out_sample(const string& fmt);  

	protected:

		/**
		 * @brief  get string output file 
		 * @param[in] fmt file format
		 * @return string : string output file 
		 */
		string _outputs(const string& fmt);

		/**
		 * @brief  get all string output file
		 * @return set<string> : all string output file
		 */
		set<string> _oformats();

		set<OFMT> _OFMT_supported;                  ///< vector of supported OFMTs (app-specific)
		bool  _append;                              ///< append mode
		int   _verb;                                ///< output verbosity   
		int   _upd;                                 ///< update [min] for output file update
		int   _len;                                 ///< length [min] for output file content
		float _smp;                                 ///< sample [sec] for output file data

	private:

	};

} // namespace


#endif
