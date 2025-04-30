
/**
* @verbatim
	History
	2012-05-11  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gstat3d.h
* @brief       Purpose: statistical function (3D)
* @author      JD
* @version     1.0.0
* @date        2012-05-11
*
*/

#ifndef GSTAT3D_H
#define GSTAT3D_H 

#include "gutils/gstat.h"
#include "gutils/gtriple.h"
#include "gexport/ExportLibGnut.h"

#include <vector>

using namespace std;

namespace gnut {
	/** @brief class for t_gstat3d. */
	class t_gstat3d : public t_gstat
	{
	public:
		/** @brief default constructor. */
		t_gstat3d(vector<t_gtriple>& data, double cint = CONF_INT);
		~t_gstat3d();

		virtual int calc_stat(double sig = 0.0);

		t_gtriple  get_rms3d() { return _rms3d; }
		t_gtriple  get_var3d() { return _var3d; }
		t_gtriple  get_sdev3d() { return _sdev3d; }
		t_gtriple  get_mean3d() { return _mean3d; }
		t_gtriple  get_median3d() { return _medi3d; }

	protected:
		void         _add_data(vector<t_gtriple>& data);

		virtual void _clear();

		// calculate statistics with outlier detection and set validity status
		virtual int  _statistics(double sig = 0.0);

		// internal purposes only (no validity status changed)
		virtual int  _calc_median();
		virtual int  _calc_mean();
		virtual int  _calc_sdev();
		virtual int  _calc_rms();

		virtual int  _chk_resid3d(double sig = 0.0);

		vector<t_gtriple> _data3d;
		vector<t_gtriple> _outl3d;

		t_gtriple         _rms3d;
		t_gtriple         _var3d;
		t_gtriple         _sdev3d;
		t_gtriple         _mean3d;
		t_gtriple         _medi3d;
	};

} // namespace

#endif
