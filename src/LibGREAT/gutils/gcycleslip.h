/**
 * @file         gcycleslip.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        cycleslip
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GCYCLESLIP_H
#define GCYCLESLIP_H

#include "gset/gsetbase.h"
#include "gdata/gsatdata.h"
#include "gio/glog.h"
#include "gexport/ExportLibGREAT.h"

using namespace gnut;
namespace great
{
	/**
	* @brief        class for cycleslip
	*/
	class LibGREAT_LIBRARY_EXPORT t_gcycleslip
	{
	public:
		/** @brief default constructor. */
		t_gcycleslip();
		t_gcycleslip(t_gsetbase* set, t_glog* log);
		virtual ~t_gcycleslip();
		/** @brief set amb flag. */
		virtual void	set_amb_flag(const string& rec, const string& sat, int flag) = 0;
		/** @brief get amb flag. */
		virtual int		get_amb_flag(const string& rec, const string& sat) = 0;
		/** @brief get active amb. */
		virtual int		get_active_amb(string site) = 0;
		/** @brief set active amb. */
		virtual void    set_active_amb(string site, int active_num) = 0;
		virtual map<string, int  >		get_active_amb() = 0;
		/** @brief number of amb arc. */
		virtual int		num_of_amb_arc(const string& site, const string& prn, const t_gtime& time) = 0 ;
		virtual bool	use_of_obs(const string& site, const string& prn, const t_gtime& time) = 0;
		virtual bool	cycle_slip(const t_gsatdata& obsdata,const t_gtime& time) = 0;
		bool cycle_slip123(t_gsatdata& obsdata, const t_gtime& time) { return false; }
		virtual void    add_ambflag(string site, string sat, string description, t_gtime beg, t_gtime end) = 0;
		
	protected:
		t_glog*     _log;
		t_gsetbase* _set;
		SLIPMODEL   _slip_model;
	};
}
#endif // !G_CYCLESLIP_H
