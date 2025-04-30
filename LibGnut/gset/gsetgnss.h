/**
*
* @verbatim
	 (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  @endverbatim
*
* @file		gsetgnss.h
* @brief	implements data extraction setting
* @author   Jan Dousa
* @version	1.0.0
* @date		2012-10-23
*
*/

#ifndef GSETGNSS_H
#define GSETGNSS_H

#define XMLKEY_GNSS "gnss" ///< The defination of gnss node 

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gsys.h"
#include "gutils/gobs.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut
{
	/// The class for gnss module in xml file
	class LibGnut_LIBRARY_EXPORT t_gsetgnss : public virtual t_gsetbase
	{
	public:
		/// constructor
		t_gsetgnss();
		/// destructor
		~t_gsetgnss();

		/// settings check
		void check() override;
		/// settings help
		void help() override;

		/**
		 * @brief get sats of all systems 
		 * @return set<string> : sats of all systems
		 */
		set<string>  sat();       

		/**
		 * @brief get sats of single system def=true (give default values at least)
		 * @param[in] gsys system
		 * @param[in] def  if empty, return default setting or not. default is true
		 * @return set<string> : sats of single system def=true (give default values at least)
		 */
		set<string>  sat(GSYS gsys, bool def = true);  

		/**
		 * @brief get obs(eg.C1C ..) of single system def=true (give default values at least)
	     * @param[in] gsys system
		 * @param[in] def  if empty, return default setting or not. default is true
		 * @return set<string> : obs of single system def=true (give default values at least)
		 */
		set<string>  obs(GSYS gsys, bool def = true);  

		/**
		 * @brief get nav of single system def=true (give default values at least)
		 * @param[in] gsys system
		 * @param[in] def  if empty, return default setting or not. default is true
		 * @return set<string> : nav of single system def=true (give default values at least)
		 */
		set<string>  nav(GSYS gsys, bool def = true); 

		/**
		 * @brief get extending gobs list with complete singals
		 * @param[in] gsys system
		 * @return set<string> : extending gobs list with complete singals
		 */
		set<string> gobs(GSYS gsys); 

		/**
		 * @brief get nav of system:gsys
		 * @param[in] gsys system
		 * @return vector<GNAVTYPE> : nav of system:gsys
		 */
		vector<GNAVTYPE> gnav(GSYS gsys);

		/**
		 * @brief get obs type of system:gsys
		 * @param[in] gsys system
		 * @return vector<GOBSTYPE> : obs type of system:gsys
		 */
		vector<GOBSTYPE> type(GSYS gsys);

		/**
		 * @brief get obs band of system:gsys
		 * @param[in] gsys system
		 * @return vector<GOBSBAND> : obs band of system:gsys
		 */
		vector<GOBSBAND> band(GSYS gsys);


		/**
		 * @brief get the freq used in proc
		 * @param[in] gsys system
		 * @return vector<FREQ_SEQ> : GNSS freq Sequence
		 */
		vector<FREQ_SEQ> freqs(GSYS gsys);

		/**
		 * @brief get the band order
		 * @param[in] gsys system
		 * @return map<FREQ_SEQ, GOBSBAND> : band index
		 */
		map<FREQ_SEQ, GOBSBAND> band_index(GSYS gsys);

		/**
		* @brief get the freq order
		* @param[in] gsys system
		* @return map<FREQ_SEQ, GOBSBAND> : freq index
		*/
		map<GOBSBAND, FREQ_SEQ> freq_index(GSYS gsys);

		/**
		 * @brief get obs ATTR of system:gsys
		 * @param[in] gsys system
		 * @return vector<GOBSATTR> : obs ATTR of system:gsys
		 */
		vector<GOBSATTR> attr(GSYS gsys);

		set<GSYS> SYS_SUPPORTED();

		/**
		 * @brief get sigma value of L obs for gsys
		 * @param[in] gsys system
		 * @return double : sigma value of L obs for gsys
		 */
		double           sigma_L(GSYS gsys);

		/**
		 * @brief get sigma value of C obs for gsys
		 * @param[in] gsys system
		 * @return double : sigma value of C obs for gsys
		 */
		double           sigma_C(GSYS gsys);

		/**
		 * @brief get sigma value of D obs for gsys
		 * @param[in] gsys system
		 * @return double : sigma value of D obs for gsys
		 */
		double           sigma_D(GSYS gsys);

		/**
		 * @brief get max res of L obs for gsys
		 * @param[in] gsys system
		 * @return double : max res of L obs for gsys
		 */
		double           maxres_L(GSYS gsys);

		/**
		 * @brief get max res of C obs for gsys
		 * @param[in] gsys system
		 * @return double : max res of C obs for gsys
		 */
		double           maxres_C(GSYS gsys);

		/**
		 * @brief get max res of D obs for gsys
		 * @param[in] gsys system
		 * @return double : max res of D obs for gsys
		 */
		double           maxres_D(GSYS gsys);

	protected:
		/**
		 * @brief get obs type of system:gsys form XML
		 * @param[in] gsys system
		 * @return vector<GOBSTYPE> : obs type of system:gsys form XML
		 */
		vector<GOBSTYPE> _type(GSYS gsys);

		/**
		 * @brief get obs band of system:gsys form XML
		 * @param[in] gsys system
		 * @return vector<GOBSBAND> : obs band of system:gsys form XML
		 */
		vector<GOBSBAND> _band(GSYS gsys);

		/**
		 * @brief get obs band of system:gsys form XML
		 * @param[in] gsys system
		 * @return vector<GOBSBAND> : obs band of system:gsys form XML
		 */
		vector<FREQ_SEQ> _sysfreq(GSYS gsys);
		
		/**
		 * @brief get the freq used in proc
		 * @param[in] gsys system
		 * @return vector<FREQ_SEQ> : GNSS freq Sequence
		 */
		vector<FREQ_SEQ> _freqs(GSYS gsys);

		/**
		 * @brief get obs attributes of system:gsys form XML
		 * @param[in] gsys system
		 * @return vector<GOBSATTR> : obs attributes of system:gsys form XML
		 */
		vector<GOBSATTR> _attr(GSYS gsys);

		/**
		 * @brief change GSYS to string::gsys
		 * @param[in] gsys system
		 * @return string : gsys
		 */
		string           _gsys(GSYS gsys);

		/**
		 * @brief get the double value of sys.attribute(sigma_L)
		 * @param[in] gsys system
		 * @return string : the double value of sys.attribute(sigma_L)
		 */
		double           _sigma_L(GSYS gsys);

		/**
		 * @brief get the double value of sys.attribute(sigma_C)
		 * @param[in] gsys system
		 * @return string : the double value of sys.attribute(sigma_C)
		 */
		double           _sigma_C(GSYS gsys);

		/**
		 * @brief get the double value of sys.attribute(sigma_D)
		 * @param[in] gsys system
		 * @return string : the double value of sys.attribute(sigma_D)
		 */
		double           _sigma_D(GSYS gsys);

		/**
		 * @brief get the double value of sys.attribute(maxres_L)
		 * @param[in] gsys system
		 * @return string : the double value of sys.attribute(maxres_L)
		 */
		double           _maxres_L(GSYS gsys);

		/**
		 * @brief get the double value of sys.attribute(maxres_C)
		 * @param[in] gsys system
		 * @return string : the double value of sys.attribute(maxres_C)
		 */
		double           _maxres_C(GSYS gsys);

		/**
		 * @brief get the double value of sys.attribute(maxres_D)
		 * @param[in] gsys system
		 * @return string : the double value of sys.attribute(maxres_D)
		 */
		double           _maxres_D(GSYS gsys);

	protected:
		map<GSYS, vector<string> > _band_str;           ///< default set
		map<GSYS, vector<string> > _type_str;           ///< default set
		map<GSYS, vector<string> > _attr_str;           ///< default set

		map<GSYS, t_gpair>        _sigma_def;          ///< default set
		map<GSYS, t_gpair>        _maxres_def;         ///< default set
		map < GSYS, double > _sigma_def_doppler;       ///< default set
		map < GSYS, double > _maxres_def_doppler;      ///< default set

	private:
	};

} // namespace

#endif
