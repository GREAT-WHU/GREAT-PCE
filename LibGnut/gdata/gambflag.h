/**
*
* @verbatim
History
-1.0 BoWong  2019-02-25  creat the file.
-1.1 BoWong  2019-04-08  Adding Doxygen Style Code Remarks
@endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		gambflag.h
* @brief	Storage the XXXXddd0.yyo.log/log13 files' data(only one site)
*				XXXX  ---- SITE name
*				 ddd  ---- Doy of the file
*				  yy  ---- year
*.
* @author   BoWong, Wuhan University
* @version	1.0.0
* @date		2019-04-08
*
*/

#ifndef GAMBFLAG_H
#define GAMBFLAG_H

#include <stdio.h>
#include <string>

#include "gdata/gobj.h"
#include "gdata/gdata.h"
#include "gexport/ExportLibGnut.h"

using namespace std;
using namespace gnut;

namespace great
{

	/**
	*@brief	   Class for one XXXXddd0.yyo.log/log13 file head data storaging
	*/
	class ambflag_hd
	{
	public:

		/** @brief default constructor. */
		ambflag_hd() {};

		/** @brief default destructor. */
		virtual ~ambflag_hd() {};

		int     beg_mjd = 0;        ///< start time mjd
		double  beg_sod = 0.0;      ///< start time sod
		int     duration = 0;        ///< duration
		double  intv = 30.0;         ///< sampling interval
		int     max_amb_1epo = 0;    ///< max ambc in one epoch
		int     old_rm_obs = 0;      ///< old remvoed observations
		int     new_rm_obs = 0;      ///< new removed observations
		int     exist_amb = 0;       ///< existed ambiguities
		int     new_amb = 0;         ///< new ambiguities
		int     avaiable_obs = 0;    ///< available observations

		string      rec_name;
		set<string> sys_name;
		set<string> sat_name;
		map<string, string> freq_name; // get freq by band
		map<string, string> band_name; // get band by freq
	};


	/**
	*@brief	   Class for one XXXXddd0.yyo.log/log13 file one satellite data storaging
	*/
	class ambflag_data {
	public:

		/** @brief default constructor. */
		ambflag_data() {};

		/** @brief default destructor. */
		virtual ~ambflag_data() {};

		string identify;        ///< identification including AMB,DEL and BAD
		int    beg_epo = 0;     ///< begin epoch
		int    end_epo = 0;     ///< end epoch
		int    nmb_arc = 0;     ///< number of arc
		string iflag;           ///< flag for compute
		double C1 = 0.0;        ///< coefficient C1
		double C2 = 0.0;        ///< coefficient C2
		string reason;          ///< cause of ambiguity 
		string log_type;        ///< log type : L12 L13 L14 L15 or more
		//add by xiong for IRC
		int yr = 0;
		int mon = 0;
		int day = 0;
		int hr = 0;
		int min = 0;
		int sec = 0;
		double C3 = 0.0;
		double C4 = 0.0;
		double C5 = 0.0;
	};

	/** vector container for storaging ambflag data ptr satellite/one station/all epoch-arc in a log file */
	typedef vector<shared_ptr<ambflag_data>> t_vec_ambflag;

	/** map container using satellite name as a index for storaging t_vec_ambflag data of one station */
	typedef map<string, t_vec_ambflag >      t_map_sat_ambflag;

	/**
	*@brief	   Class for one XXXXddd0.yyo.log/log13 file data storaging
	*/
	class LibGnut_LIBRARY_EXPORT t_gambflag : public t_gdata
	{

	public:

		/** @brief default constructor. */
		t_gambflag();

		/** @brief default destructor. */
		virtual ~t_gambflag();

		/**
		* @brief set the file's belonging, it's station name.
		* @param[in] site site name
		*/
		void setSite(string site) { _site = site; }

		/**
		* @brief add head data of ambflag file.
		* @param[in]  head		ambflag file head data.
		* @return	  void
		*/
		void setAmbFlagHead(ambflag_hd head);

		/**
		* @brief reset the iflag of a data of one station/one satellite/one epoch-arc.
		* @param[in]  prn		satellite name.
		* @param[in]  flag		iflag in class ambflag_data.
		* @param[in]  pos		the serial number of epoch-arc where current epoch in.
		* @return	  void
		*/
		void reset_iflag(const string& prn, const string& flag, const int& pos);

		/**
		* @brief get ambflag head data of one file(one station datum's head).
		* @return ambflag head data
		*/
		shared_ptr<ambflag_hd>& getAmbFlagHead();

		/**
		* @brief get t_map_sat_ambflag data of one station/one satellite.
		* @return	 t_map_sat_ambflag data of one station/one satellite
		*/
		t_map_sat_ambflag& getAmbFlagData();

		t_vec_ambflag& getSatAmbFlag(const string& sat);
		/**
		* @brief add ambflag data of one station/one satellite/one epoch-arc.
		* @param[in]  prn		satellite name.
		* @param[in]  data		ambflag data of one station/one satellite/one epoch-arc.
		* @return	  void
		*/
		void addAmbFlagData(string prn, const ambflag_data& data);

		/**
		* @brief judge if an epoch-arc data is usable.
		* @param[in]  prn		satellite name.
		* @param[in]  time		current epoch.
		* @param[in]  pos		the serial number of epoch-arc where current epoch in.
		* @param[out] pos
		* @return	  true is available/ flase not available
		*/
		bool isValid(const string& prn, const t_gtime& time, int& pos);

		int get_amb_pos(const string& prn, const t_gtime& time, bool apply_carrier_range);

		t_gtime epoch2time(const int &epo);

		set<string> getAllSatSet();

		t_gtime getMaxEndEpoch();

		void clearAllAmbFlagData();

	protected:
		string                          _site;              ///< site name
		shared_ptr<ambflag_hd>          _ambflag_head;      ///< ptr of ambflag head which storaging log file head of one station
		t_map_sat_ambflag               _all_sats_ambflag;  ///< map container using satellite name as a index for storaging 
															///<  t_vec_ambflag data of one station
	private:
	};

} // namespace

#endif