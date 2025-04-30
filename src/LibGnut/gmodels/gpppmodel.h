
/**
*
* @verbatim
    History
    2012-09-20  PV: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gpppmodel.h
* @brief      Purpose: various PPP models
*.
* @author     PV
* @version    1.0.0
* @date       2012-09-20
*
*/

#ifndef GPPPMODEL_H
#define GPPPMODEL_H


#include <string>
#include <map>
#include <cmath>

#include "gmodels/gsppmodel.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gephplan.h"
#include "gall/gallobj.h"

using namespace std;

namespace gnut
{
	/** @brief class for t_gpppmodel derive from t_gsppmodel. */
	class LibGnut_LIBRARY_EXPORT t_gpppmodel : public t_gsppmodel
	{
	public:
		/** @brief constructor 1. */
		t_gpppmodel(string site, t_glog* glog, t_gsetbase* settings);
		t_gpppmodel();
		virtual ~t_gpppmodel();

		virtual double windUp(t_gsatdata& satdata, const ColumnVector&);
		virtual double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata&, t_gobs& gobs, bool com);
		virtual double tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple ell, t_gsatdata& satdata);

		/** @brief set object. */
		virtual void setOBJ(t_gallobj* obj);

		// attitude modeling - public interface
		int  attitude_old(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k); // from RTKlib (to remove)
		int  attitude(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		int  attitude(t_gsatdata& satdata, double yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		int  attitude(string antype, string prn, ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, ColumnVector& i, ColumnVector& j, ColumnVector& k);

	protected:

		// From RTKlib - needs to be removed
		int  _yaw(t_gsatdata& satdata, string antype, ColumnVector& xs, ColumnVector& ys, ColumnVector& zs);

		// attitude niminal modeling
		void _ysm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _ysm(string prn, double bata, double mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _onm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _onm(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _noon_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _midnight_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _noon_turn(string _prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k);
		void _midnight_turn(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k);


		// attitude for GPS Block IIA
		void _attitude_GPSIIA(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_GPSIIA(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _midnight_turn_GPSIIA(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _midnight_turn_GPSIIA(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		// attitude for GPS Block IIR
		void _attitude_GPSIIR(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_GPSIIR(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		// attitude for GPS Block IIR-M
		void _attitude_GPSIIRM(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);

		// attitude for GPS Block IIF
		void _attitude_GPSIIF(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_GPSIIF(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _midnight_turn_GPSIIF(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _midnight_turn_GPSIIF(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		//atttude for GPS Block III
		void _attitude_GPSIII(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_GPSIII(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		// attitude for Galileo IOV
		void _attitude_GAL1(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_GAL1(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _noon_turn_GAL1(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _noon_turn_GAL1(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		// attitude for Galileo FOC
		void _attitude_GAL2(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_GAL2(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _noon_turn_GAL2(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _noon_turn_GAL2(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		// Continuous yaw steering attitude modes of BDS satellites
		void _cys_cast(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _cys_cast(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);
		void _cys_secm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _cys_secm(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);
		// attitude for BeiDou
		void _attitude_BDS(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_BDS(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _cys_qzs(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _cys_qzs(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);
		void _switch_qzs1(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _switch_qzs1(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		// attitude for QZSS
		void _attitude_QZS(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_QZS(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		// attitude for GLO
		void _attitude_GLO(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _attitude_GLO(string antype, string prn, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, ColumnVector & i, ColumnVector & j, ColumnVector & k);

		void _midnight_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _noon_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _midnight_turn_GLOM(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k);
		void _noon_turn_GLOM(string prn, double _beta, double _mi, ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun, double R, ColumnVector & i, ColumnVector & j, ColumnVector & k);


		void _yaw2ijk(t_gsatdata& satdata, double& yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k);
		void _yaw2ijk(ColumnVector& xsat, ColumnVector& vsat, ColumnVector& xsun, double& yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k);


		double _orb_angle(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun);
		double _beta(ColumnVector & xsat, ColumnVector & vsat, ColumnVector & xsun);

		double sign(double a, double b);

		void _set_beta0(double beta);
		double _get_beta0();


		map<string, double> _windUpTime;
		map<string, double> _windUpSum;
		t_gephplan          _ephplan;
		t_gallobj*          _gallobj;
		GRDMPFUNC           _grad_mf;
		ATTITUDES		   _attitudes;

		map<string, double>  _last_beta;
		map<string, double>  _last_yaw;
		map<string, t_gtime> _last_epo;

		double _beta0 = 0;
	};

} // namespace

#endif //  GPPPMODEL_H
