/**
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  @endverbatim
*
* @file		gsetgen.h
* @brief	implements common general settings
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#ifndef GSETGEN_H
#define GSETGEN_H

#include "gexport/ExportLibGnut.h"
#include <set>
#include <string>
#include <iostream>
#include <vector>

#include "gio/glog.h"
#include "gutils/gsys.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"

#define XMLKEY_GEN   "gen" ///< The defination of gen node
#define DEF_RECEIVER "   " ///< Default receiver : all !
#define DEF_SAMPLING 30    ///< Default sampling : 30s !

using namespace std;

namespace gnut
{

	class LibGnut_LIBRARY_EXPORT t_gobj;

	class LibGnut_LIBRARY_EXPORT t_gsetgen : public virtual t_gsetbase
	{
	public:

		/**
		 * @brief default constructor, distinguish GNSS/nonGNSS app
		 * @param[in] gnss : gnss process or not, default is true
		 */
		t_gsetgen(bool gnss = true);                 

		/**@brief destructor */
		~t_gsetgen();

		/**@brief settings check */
		void check() override;

		/**@brief settings help */
		void help() override;

		/**
		 * @brief get the beg time of process
		 * @return t_gtime : the beg time of process
		 */
		t_gtime beg(bool conv = true);

		/**
		 * @brief get the end time of process
		 * @return t_gtime : the end time of process
		 */
		t_gtime end(bool conv = true);

		/**
		 * @brief get the sampling time of process
		 * @return double : the sampling time of process
		 */
		double sampling();

		/**
		 * @brief get the default sampling time of process
		 * @return double : the default sampling time of process
		 *  @retval DEF_SAMPLING default sampling 
		 */
		double sampling_default() { return DEF_SAMPLING; }       

		/**
		 * @brief get the decimals scale of sampling time in process
		 * @return int : the decimals scale of sampling time in process
		 */
		int    sampling_scalefc() { return (int)pow(10, _dec); } 

		/**
		 * @brief get the decimals for sampling interval (for high-rate) in process
		 * @return int : decimals for sampling interval (for high-rate) in process
		 */
		int    sampling_decimal() { return _dec; }

		/**
		 * @brief get the List of system names
		 * @return set<string> : List of system names
		 */
		virtual set<string> sys();

		set<string> recs();
		set<string> rec_all();

		virtual vector<string> list_base();
		virtual vector<string> list_rover();

		// Return crd of site
		virtual vector<double> crd(string site); 

		// add for clk est.
		virtual string refsat();  // ref sat clk
		virtual string refsite(); // ref site clk
		virtual double sig_refclk(); // sig_refclk;

		/**
		 * @brief get name of estimator
		 * @return string : estimator name
		 */
		virtual string estimator();

		/**
		* @brief add for remove unused satellites
		* @return set<string> : satellites which will be removed
		*/
		virtual set<string> sat_rm();

	protected:

		bool   _gnss; ///< gnss or not
		string _sys;  ///< sys name
		int    _dec;  ///< sampling 

	private:
	};

} // namespace

#endif
