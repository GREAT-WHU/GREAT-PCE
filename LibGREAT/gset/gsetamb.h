/**
 * @file         gsetamb.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GSETAMB_H
#define GSETAMB_H

#define XMLKEY_AMBIGUITY "ambiguity"

#include "gexport/ExportLibGREAT.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace gnut;
namespace great
{
	enum class AMB_TYPE {
		UD, // undifferenced
		SD, // Single-differenced, usually between satellites
		DD, // double-differenced
		UNDEF
	};

	enum class WL_MODE {
		OBS,  // get WL amb from RINEXO
		WM,   // get WL amb from MW obs
		AMB   // get WL amb from ambinp
	};

	/** @brief enum of FIX mode. */
	enum class FIX_MODE {
		NO,
		ROUND,
		SEARCH,
		HOLD
	};

	/** @brief enum of UPD mode. */
	enum class UPD_MODE {
		UPD,  // wl upd + nl upd
	};

	/**
	*@brief	   Class for set ambiguity fixed xml
	*/
	class LibGREAT_LIBRARY_EXPORT t_gsetamb : public virtual t_gsetbase
	{
	public:

		/** @brief default constructor. */
		t_gsetamb();

		/** @brief default destructor. */
		virtual ~t_gsetamb();

		/**
		* @brief settings check.
		* @return void
		*/
		void check();

		/**
		* @brief settings help.
		* @return void
		*/
		void help();

		/**
		* @brief  get ambiguity fixing mode .
		* @return ambiguity fixing mode
		*/
		FIX_MODE fix_mode();

		/**
		* @brief  get upd mode .
		* @return upd mode
		*/
		UPD_MODE upd_mode();

		/**
		* @brief  get lambda ratio .
		* @return lambda ratio
		*/
		double lambda_ratio();

		/**
		* @brief  get bootstrapping rate.
		* @return bootstrapping rate
		*/
		double bootstrapping();

		/**
		* @brief  get minimum common time of two observation arc.
		* @return minimum common time of two observation arc
		*/
		double min_common_time();

		/**
		* @brief   get ambiguity decision.
		* @return  ambiguity decision.
		*/
		map<string, double> get_amb_decision(string str);

		/**
		* @brief   whether take partial ambiguity fixed mode.
		* @return  true or false.
		*/
		bool part_ambfix();

		/**
		* @brief   value's size which take partial ambiguity fixed mode.
		* @return  true or false.
		*/
		int part_ambfix_num();

		/**
		* @brief change string to FIX_MODE.
		* @param[in]  str   fix mode string
		* @return	  FIX_MODE
        */
		FIX_MODE str2fixmode(string str);

		/**
		* @brief change string to UPD_MODE.
		* @param[in]  str   upd mode string
		* @return	  UPD_MODE
		*/
		UPD_MODE str2upd_mode(string str);

		bool isSetRefSat();

	protected:
		map<string, map<string, double>> _default_decision = {
			{"EWL", {{"maxdev", 0.07}, {"maxsig", 0.10}, {"alpha", 1000}}},
			{  "WL", {{"maxdev", 0.25}, {"maxsig", 0.10}, {"alpha", 1000}}},
			{   "NL", {{"maxdev", 0.25}, {"maxsig", 0.10}, {"alpha", 1000}}}
		};
	private:

	};
}
#endif
