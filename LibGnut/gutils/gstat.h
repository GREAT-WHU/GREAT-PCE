
/**
* @verbatim
	History
	2012-05-11  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gstat.h
* @brief       Purpose: statistical function (1D)
* @author      JD
* @version     1.0.0
* @date        2012-05-11
*
*/

#ifndef GSTAT_H
#define GSTAT_H 

#include "gdata/gdata.h"
#include "gutils/gpair.h"
#include "gexport/ExportLibGnut.h"

#include <map>
#include <vector>

#define CONF_INT 3    // confident interval factor

using namespace std;

namespace gnut {
	/** @brief class for t_gdata. */
	class t_gstat : public t_gdata
	{
	public:

		typedef map<t_gpair, int>  t_hist;
		/** @brief default constructor. */
		t_gstat(double cint = CONF_INT);                           // just construct
		t_gstat(vector<double>& data, double cint = CONF_INT);     // construct + calculate statistics
		virtual ~t_gstat();

		void add_data(vector<double>& data);                       // clear + calculate statistics

		virtual bool valid() { return _valid; }                     // get status

		virtual int calc_stat(double sig = 0.0);                   // calculate statistics

		virtual int calc_quartiles(double& low, double& upp);
		virtual int calc_iqrlimits(double& low, double& upp);

		virtual double get_min() { return _min; }
		virtual double get_max() { return _max; }
		virtual double get_var() { return _var; }
		virtual double get_rms() { return _rms; }
		virtual double get_sdev() { return _sdev; }
		virtual double get_mean() { return _mean; }
		virtual double get_median() { return _medi; }

		virtual int    get_size() { return _data.size(); }
		virtual int    get_outl() { return _outl.size(); }

		t_hist histogram(vector<double>& data, set<double>& bound);

	protected:

		void         _add_data(vector<double>& data);

		virtual void _clear();

		// calculate statistics with outlier detection and set validity status
		virtual int  _statistics(double sig = 0.0);

		// internal purposes only (no validity status changed)
		virtual int  _calc_minmax();
		virtual int  _calc_median();
		virtual int  _calc_mean();
		virtual int  _calc_sdev();
		virtual int  _calc_rms();

		virtual int  _chk_resid(double sig = 0.0); // check residuals (optional use of external sigma)

		bool _valid;

		vector<double> _data;
		vector<double> _outl;

		double _cint;
		double _sdev;
		double _var;
		double _rms;
		double _mean;
		double _medi;
		double _min;
		double _max;
	};

} // namespace

#endif
