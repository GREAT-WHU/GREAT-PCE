/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "gdata/gsatdata.h"
#include "gmodels/gephplan.h"
#include "gutils/gtypeconv.h"
#include "gutils/gsysconv.h"
#include "gutils/gmatrixconv.h"

using namespace std;

namespace gnut {
	t_gsatdata::t_gsatdata()
	{
	}
	// -----------
	t_gsatdata::t_gsatdata(string site, string sat, const t_gtime& t)
		: t_gobsgnss(site, sat, t),
		_satcrd(0.0, 0.0, 0.0),
		_satpco(0.0, 0.0, 0.0),
		_dloudx(0.0, 0.0, 0.0),
		_drate(0.0),
		_clk(0.0),
		_dclk(0.0),
		_ele(0.0),
		_azi(0.0),
		_rho(0.0),
		_eclipse(false),
		_mfH(0.0),
		_mfW(0.0),
		_mfG(0.0),
		_wind(0.0),
		_low_prec(false),
		_slipf(false),
		_beta_val(999),
		_orb_angle_val(999),
		_yaw(999)

	{
		id_type(t_gdata::SATDATA);
		id_group(t_gdata::GRP_OBSERV);
	}


	// -----------
	t_gsatdata::t_gsatdata(const t_gobsgnss& obs)
		: t_gobsgnss(obs),
		_satcrd(0.0, 0.0, 0.0),
		_satpco(0.0, 0.0, 0.0),
		_clk(0.0),
		_ele(0.0),
		_azi(0.0),
		_rho(0.0),
		_eclipse(false),
		_mfH(0.0),
		_mfW(0.0),
		_mfG(0.0),
		_wind(0.0),
		_low_prec(false),
		_slipf(false),
		_beta_val(999),
		_orb_angle_val(999),
		_yaw(999)
	{
		id_type(t_gdata::SATDATA);
		id_group(t_gdata::GRP_OBSERV);
	}


	// -----------
	t_gsatdata::~t_gsatdata()
	{

	}

	// -----------
	void t_gsatdata::addpco(const t_gtriple& pco)
	{
		gtrace("t_gsatdat::addpco");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_satpco = pco;
		_gmutex.unlock();
		return;
	}

	// -----------
	void t_gsatdata::addcrd(const t_gtriple& crd)
	{
		gtrace("t_gsatdat::addcrd");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_satcrd = crd;
		_gmutex.unlock();
		return;
	}

	void t_gsatdata::addcrdcrs(const t_gtriple& crd)
	{
		_satcrdcrs = crd;
		return;
	}

	// -----------
	void t_gsatdata::addvel(const t_gtriple& vel)
	{
		gtrace("t_gsatdat::addvel");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_satvel = vel;
		_gmutex.unlock();
		return;
	}

	// -----------
	void t_gsatdata::addvel_crs(const t_gtriple& vel)
	{
		gtrace("t_gsatdat::addvel");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_satvel_crs = vel;
		_gmutex.unlock();
		return;
	}

	//----------
	int t_gsatdata::addprd(t_gallnav* gnav, bool corrTOT, bool msk_health)
	{
		gtrace("t_gsatdat::addprd");

		//_gmutex.lock();

		_low_prec = false;

		int irc = this->_addprd(gnav, corrTOT, msk_health);

		//_gmutex.unlock(); 
		return irc;
	}


	int t_gsatdata::addprd(const t_gtime T_sat, t_gallnav* gnav, bool corrTOT, bool msk_health)
	{
		gtrace("t_gsatdat::addprd");

		_gmutex.lock();

		_low_prec = false;

		int irc = this->_addprd(T_sat, gnav, corrTOT, msk_health);

		_gmutex.unlock(); return irc;
	}

	//----------
	int t_gsatdata::addprd_nav(t_gallnav* gnav, bool corrTOT, bool msk_health)
	{
		gtrace("t_gsatdat::addprd_nav");

		_gmutex.lock();

		_low_prec = true;

		int irc = this->_addprd(gnav, corrTOT, msk_health);

		_gmutex.unlock(); return irc;
	}


	// Compute rho, ele, azi ...
	// -----------------------
	int t_gsatdata::cmpVal(t_gtriple& xyz)
	{
		gtrace("t_gsatdata::cmpVal");

		if (double_eq(_satcrd[0], 0.0) ||
			double_eq(_satcrd[1], 0.0) ||
			double_eq(_satcrd[2], 0.0)) {
			cerr << "Satellite position has not been calculated" << endl;
			return -1;
		}


		if (double_eq(xyz[0], 0.0) ||
			double_eq(xyz[1], 0.0) ||
			double_eq(xyz[2], 0.0)) {
			cerr << "Station position has not been calculated" << endl;
			return -1;
		}

		t_gtriple neu_sat;
		t_gtriple ell_sit;
		xyz2ell(xyz, ell_sit, false);
		t_gtriple xyz_rho = _satcrd - xyz;

		xyz2neu(ell_sit, xyz_rho, neu_sat);

		// Correct Earth rotation
		ColumnVector xRec(3);
		double rho0 = (_satcrd.crd_cvect() - xyz.crd_cvect()).norm_Frobenius();
		double dPhi = OMEGA * rho0 / CLIGHT;

		xRec(1) = xyz[0] * cos(dPhi) - xyz[1] * sin(dPhi);
		xRec(2) = xyz[1] * cos(dPhi) + xyz[0] * sin(dPhi);
		xRec(3) = xyz[2];

		double tmp = (_satcrd.crd_cvect() - xRec).norm_Frobenius();
		_rho = tmp;

		double NE2 = neu_sat[0] * neu_sat[0] + neu_sat[1] * neu_sat[1];
		double ele = acos(sqrt(NE2) / _rho);
		if (neu_sat[2] < 0.0) {
			ele *= -1.0;
		}

		if (sqrt(NE2) / _rho > 1.0)
			_ele = 0.0;
		else _ele = ele;

		double azi = atan2(neu_sat[1], neu_sat[0]);
		if (azi < 0) azi += 2 * G_PI;
		_azi = azi;

		return 1;
	}



	// -----------
	void t_gsatdata::addclk(double clk)
	{
		gtrace("t_gsatdat::addclk");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_clk = clk;
		_gmutex.unlock();
	}

	void t_gsatdata::addnavclk(double d)
	{
		gtrace("t_gsatdat::addnavclk");
		_gmutex.lock();
		_navclk = d;
		_gmutex.unlock();
	}

	// -----------
	void t_gsatdata::addreldelay(double rel)
	{
		gtrace("t_gsatdat::reldelay");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_reldelay = rel;
		_gmutex.unlock();
	}

	void t_gsatdata::addreccrd(const t_gtriple& crd)
	{
		_reccrd = crd;
	}
	void t_gsatdata::addreccrdcrs(const t_gtriple& crd)
	{
		_reccrdcrs = crd;
	}
	void t_gsatdata::addsat2reccrs(const t_gtriple& crd)
	{
		_sat2reccrs = crd;
	}

	void t_gsatdata::addTRS2CRS(const Matrix& rotmat, const Matrix& drdxpole, const Matrix& drdypole, const Matrix& drdut1)
	{
		_rotmat = rotmat;
		_drdxpole = drdxpole;
		_drdypole = drdypole;
		_drdut1 = drdut1;
	}

	void t_gsatdata::addSCF2CRS(const Matrix& scf2crs, const Matrix& scf2trs)
	{
		_scf2crs = scf2crs;
		_scf2trs = scf2trs;
	}


	void t_gsatdata::addorbfunct(const Matrix& orbfunct)
	{
		_orbfunct = orbfunct;
	}

	void t_gsatdata::adddrate(double drate)
	{
		_drate = drate;
	}

	void t_gsatdata::adddloudx(const t_gtriple& unit)
	{
		_dloudx = unit;
	}

	void t_gsatdata::addsatindex(int idx)
	{
		_satindex = idx;
	}

	void t_gsatdata::addrecTime(const t_gtime& recTime)
	{
		_TR = recTime;
	}
	void t_gsatdata::addsatTime(const t_gtime& satTime)
	{
		_TS = satTime;
	}

	// -----------
	void t_gsatdata::addele(double ele)
	{
		gtrace("t_gsatdat::addele");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_ele = ele;
		_gmutex.unlock();
	}

	void t_gsatdata::addele_leo(double ele)
	{
		gtrace("t_gsatdat::addele");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_ele_leo = ele;
		_gmutex.unlock();
	}


	// -----------
	void t_gsatdata::addazi(double azi)
	{
		gtrace("t_gsatdat::addazi");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_azi = azi;
		_gmutex.unlock();
	}

	void t_gsatdata::addazi_sat(double azi_sat)
	{
		gtrace("t_gsatdat::addazi");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_azi_sat = azi_sat;
		_gmutex.unlock();
	}

	void t_gsatdata::addnadir(double nadir)
	{
		gtrace("t_gsatdat::addazi");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_nadir = nadir;
		_gmutex.unlock();
	}


	// -----------
	void t_gsatdata::addrho(double rho)
	{
		gtrace("t_gsatdat::addrho");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_rho = rho;
		_gmutex.unlock();
	}

	// -----------
	void t_gsatdata::addmfH(const double& mfH)
	{
		gtrace("t_gsatdat::addmfH");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_mfH = mfH;
		_gmutex.unlock();
	}

	// -----------
	void t_gsatdata::addmfW(const double& mfW)
	{
		gtrace("t_gsatdat::addmfW");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_mfW = mfW;
		_gmutex.unlock();
	}

	// -----------
	void t_gsatdata::addmfG(const double& mfG)
	{
		gtrace("t_gsatdat::addmfG");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_mfG = mfG;
		_gmutex.unlock();
	}


	// -----------
	t_gtriple t_gsatdata::satcrd() const
	{
		gtrace("t_gsatdat::satcrd");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif

		return _satcrd;
	}

	t_gtriple t_gsatdata::satcrdcrs() const
	{
		return _satcrdcrs;
	}
	// -----------
	t_gtriple t_gsatdata::satpco() const
	{
		gtrace("t_gsatdat::satcrd");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif

		return _satpco;
	}

	// -----------
	t_gtriple t_gsatdata::satvel() const
	{
		gtrace("t_gsatdat::satvel");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif

		return _satvel;
	}

	// -----------
	double t_gsatdata::clk() const
	{
		gtrace("t_gsatdat::clk");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _clk;
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::dclk()
	{
		gtrace("t_gsatdat::dclk");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _dclk;
		_gmutex.unlock();
		return tmp;
	}


	// -----------
	double t_gsatdata::reldelay() const
	{
		gtrace("t_gsatdat::reldelay");
		_gmutex.lock();
		double tmp = _reldelay;
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::drate()
	{
		return _drate;
	}

	// -----------
	Matrix  t_gsatdata::orbfunct()
	{
		return _orbfunct;
	}
	// -----------
	int	  t_gsatdata::satindex()
	{
		return _satindex;
	}
	// -----------
	t_gtriple t_gsatdata::reccrd()
	{
		return _reccrd;
	}
	t_gtriple t_gsatdata::reccrdcrs()
	{
		return _reccrdcrs;
	}
	t_gtriple t_gsatdata::sat2reccrs()
	{
		return _sat2reccrs;
	}
	// -----------
	Matrix t_gsatdata::rotmat()
	{
		return _rotmat;
	}

	// -----------
	Matrix t_gsatdata::drdxpole()
	{
		return _drdxpole;
	}

	// -----------
	Matrix t_gsatdata::drdypole()
	{
		return _drdypole;
	}

	// -----------
	Matrix t_gsatdata::drdut1()
	{
		return _drdut1;
	}
	Matrix t_gsatdata::scf2crs()
	{
		return _scf2crs;
	}
	Matrix t_gsatdata::scf2trs()
	{
		return _scf2trs;
	}
	// -----------
	t_gtime t_gsatdata::recTime()
	{
		gtrace("t_gsatdat::satTime");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		t_gtime tmp = _TR;
		_gmutex.unlock();
		return tmp;
	}

	// -----------
	t_gtime t_gsatdata::satTime() const
	{
		gtrace("t_gsatdat::satTime");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		t_gtime tmp = _TS;
		_gmutex.unlock();
		return tmp;
	}



	// -----------
	double t_gsatdata::ele() const
	{
		gtrace("t_gsatdat::ele");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _ele;
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::ele_leo() const
	{
		gtrace("t_gsatdat::ele");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _ele_leo;
		_gmutex.unlock();
		return tmp;
	}

	// -----------
	double t_gsatdata::ele_deg() const
	{
		gtrace("t_gsatdat::ele_deg");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _ele * 180.0 / G_PI;
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::ele_leo_deg() const
	{
		gtrace("t_gsatdat::ele_deg");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _ele_leo * 180.0 / G_PI;
		_gmutex.unlock();
		return tmp;
	}

	// -----------
	double t_gsatdata::azi()
	{
		gtrace("t_gsatdat::azi");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _azi;
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::azi_sat()
	{
		gtrace("t_gsatdat::azi_sat");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _azi_sat;
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::nadir()
	{
		gtrace("t_gsatdat::azi_sat");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _nadir;
		_gmutex.unlock();
		return tmp;
	}


	// -----------
	double t_gsatdata::rho()
	{
		gtrace("t_gsatdat::rho");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		double tmp = _rho;
		_gmutex.unlock();
		return tmp;
	}


	// valid
	// ----------
	bool t_gsatdata::valid()
	{
		gtrace("t_gsatdat::valid");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		bool tmp = this->_valid();
		_gmutex.unlock();
		return tmp;
	}

	// -----------
	void t_gsatdata::addslip(const bool& flag)
	{
		gtrace("t_gsatdat::addslip");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		_slipf = flag;
		_gmutex.unlock();
	}

	// -----------
	bool t_gsatdata::islip()
	{
		gtrace("t_gsatdat::isSlip");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		bool tmp = _slipf;
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::beta()
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif

		_gmutex.lock();
		double tmp = _b();
		_gmutex.unlock();
		return tmp;
	}

	double t_gsatdata::orb_angle()
	{
#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif

		_gmutex.lock();
		double tmp = _orb_angle();
		_gmutex.unlock();
		return tmp;
	}

	// clean data
	// ----------
	void t_gsatdata::clear()
	{
		gtrace("t_gsatdat::clear");

#ifdef BMUTEX
		boost::mutex::scoped_lock lock(_mutex);
#endif
		_gmutex.lock();
		this->_clear();
		_gmutex.unlock();
		return;
	}

	//----------
	int t_gsatdata::_addprd(t_gallnav* gnav, bool corrTOT, bool msk_health)
	{
		gtrace("t_gsatdat::_addprd");

		string satname(_satid);

		GSYS gs = this->gsys();

		GOBSBAND b1, b2;
		b1 = b2 = BAND;

		// automatic selection of two bands for IF LC
		set<GOBSBAND> bands = _band_avail_code();
		auto itBAND = bands.begin();
		if (corrTOT) {
			if (bands.size() < 1) {
				cout << "t_gsatdata bands.size() < 2" << endl;
				if (_log) _log->comment(2, "gsatdata", "At least two bands are necessary for TOT correction in sat pos/clk calculation!");
				return -1;
			}
			else if(bands.size() < 2)
			{
				b1 = *itBAND;
			}
			else
			{
				b1 = *itBAND;
				itBAND++;
				b2 = *itBAND;
			}
		}
		else
		{
			b1 = *itBAND;
		}


		double P3 = this->P3(b1, b2);
        if (double_eq(P3, 0.0)) P3 = this->obs_C(t_gband(b1, GOBSATTR::ATTR));


		//test for observations availability
		if (gnav == 0) {
			if (_log)
				_log->comment(2, "gsatdata", " satellite " + satname
					+ _epoch.str_ymdhms("  t_gallnav pointer is not available "));
			return -1;
		}

		//test for observations availability
		if (double_eq(P3, 0.0) && corrTOT) {
			if (_log)
				_log->comment(2, "gsatdata", " satellite " + satname
					+ _epoch.str_ymdhms(" P3 = 0!"));
			return -1;
		}

		double xyz[3] = { 0.0, 0.0, 0.0 };
		double vel[3] = { 0.0, 0.0, 0.0 };
		double var[3] = { 0.0, 0.0, 0.0 };
		double clk = 0.0;
		double dclk = 0.0;
		double clkrms = 0.0;

		if (satname.substr(0, 1) != "G" &&
			satname.substr(0, 1) != "R" &&
			satname.substr(0, 1) != "E" &&
			satname.substr(0, 1) != "J" &&
			//       satname.substr(0,1) != "S" &&       
			satname.substr(0, 1) != "C")
		{
			if (_log) _log->comment(2, "gsatdata", " satelite " + satname + _epoch.str_ymdhms(" Undefined satellite system! "));
			return -1;
		}

		t_gtime epoT(t_gtime::GPS);
		double satclk = 0.0;
		double satclk2 = 1.0;
		int cnt = 0;

		if (corrTOT) {
			while (fabs(satclk - satclk2) > 1.e-3 / CLIGHT) {
				satclk2 = satclk;
				epoT = _epoch - P3 / CLIGHT - satclk;

				int irc = gnav->clk(satname, epoT, &clk, &clkrms, &dclk, msk_health);

				if (irc < 0 || cnt++ > 25) {
					if (_log) _log->comment(2, "gsatdata", " satelite " + satname
						+ _epoch.str_ymdhms(" clocks not calculated (irc|iter) for epoch: "));
					return -1;
				}
				satclk = clk;
			}
		}
		else {
			epoT = _epoch;
			int irc = gnav->clk(satname, epoT, &satclk, &clkrms, &dclk, msk_health);
			if (irc < 0) {
				if (_log)
					_log->comment(2, "gsatdata", " satelite " + satname
						+ _epoch.str_ymdhms(" clocks not calculated for epoch "));
				return -1;
			}
		}

		int irc = 0;
		irc = gnav->pos(satname, epoT, xyz, var, vel, msk_health);

		if (irc < 0) {
			if (_log)
				_log->comment(2, "gsatdata", " satelite " + satname
					+ _epoch.str_ymdhms(" coordinates not calculated for epoch "));
			return -1;
		}

		t_gtriple txyz(xyz);
		t_gtriple tvel(vel);

		// relativistic correction
		if (gs != GLO ||
			(gs == GLO && gnav->id_type() == t_gdata::ALLPREC)) {
			double rel = 2.0 * (txyz[0] * vel[0] + txyz[1] * vel[1] + txyz[2] * vel[2]) / CLIGHT / CLIGHT;
			shared_ptr<t_geph> eph = gnav->find(satname, epoT);
			if (gs == BDS && eph && gnav->id_type() == t_gdata::ALLRTCM) {
				shared_ptr <t_gnavbds> gnavb = dynamic_pointer_cast <t_gnavbds>(eph);
				double clk0 = 0.0, dclk0 = 0.0, clkrms0 = 0.0;
				gnavb->clk(epoT, &clk0, &clkrms0, &dclk0, msk_health);
				rel = gnavb->rel();
			}
			else if (gs == GAL && eph && gnav->id_type() == t_gdata::ALLRTCM) {
				shared_ptr <t_gnavgal> gnave = dynamic_pointer_cast <t_gnavgal>(eph);
				double clk0 = 0.0, dclk0 = 0.0, clkrms0 = 0.0;
				gnave->clk(epoT, &clk0, &clkrms0, &dclk0, msk_health);
				rel = gnave->rel();
			}

			if (rel == 0.0) {
				if (_log)
					_log->comment(0, "gsatdata", " satelite " + satname
						+ _epoch.str_ymdhms(" relativity correction not calculated for epoch "));
				return -1;
			}

			satclk -= rel;
			_TS = epoT;
			_reldelay = rel * CLIGHT;
		}


		_satcrd = txyz;
		_satvel = tvel;
		_clk = satclk * CLIGHT;
		_dclk = dclk * CLIGHT;
		return 1;
	}

	// update sat clk
	int t_gsatdata::_addprd(const t_gtime& T_sat, t_gallnav* gnav, bool corrTOT, bool msk_health)
	{
		gtrace("t_gsatdat::_addprd");

		string satname(_satid);

		GSYS gs = this->gsys();

		if (gnav == 0) {
			if (_log)
				_log->comment(2, "gsatdata", " satellite " + satname
					+ _epoch.str_ymdhms("  t_gallnav pointer is not available "));
			return -1;
		}

		double xyz[3] = { 0.0, 0.0, 0.0 };
		double vel[3] = { 0.0, 0.0, 0.0 };
		double var[3] = { 0.0, 0.0, 0.0 };
		double satclk = 0.0;
		double dclk = 0.0;
		double clkrms = 0.0;

		if (satname.substr(0, 1) != "G" &&
			satname.substr(0, 1) != "R" &&
			satname.substr(0, 1) != "E" &&
			satname.substr(0, 1) != "J" &&
			satname.substr(0, 1) != "C")
		{
			if (_log) _log->comment(2, "gsatdata", " satelite " + satname + _epoch.str_ymdhms(" Undefined satellite system! "));
			return -1;
		}

		int irc = 0;
		irc = gnav->clk(satname, T_sat, &satclk, &clkrms, &dclk, msk_health);
		if (irc < 0) {
			if (_log)
				_log->comment(2, "gsatdata", " satelite " + satname
					+ T_sat.str_ymdhms(" clocks not calculated for epoch "));
			return -1;
		}

		irc = gnav->pos(satname, T_sat, xyz, var, vel, msk_health);

		if (irc < 0) {
			if (_log)
				_log->comment(2, "gsatdata", " satelite " + satname
					+ T_sat.str_ymdhms(" coordinates not calculated for epoch "));
			return -1;
		}

		t_gtriple txyz(xyz);
		t_gtriple tvel(vel);

		_satcrd = txyz;
		_satvel = tvel;
		_TS = T_sat;

		return 1;
	}

	int t_gsatdata::_correction(double* xyz, double* vxyz, double& clk, double& dclk, double* xyz_corr, double* vxyz_corr, double& clk_corr, double& dclk_corr)
	{
		ColumnVector pos(xyz, 3), vel(vxyz, 3), rao_pos(xyz_corr, 3), rao_vel(vxyz_corr, 3), out(0.0, 3);
		rao2xyz(pos, vel, rao_pos, out);
		for (size_t i = 0; i < 3; i++)
		{
			xyz[i] -= out(i + 1);
		}
		rao2xyz(pos, vel, rao_vel, out);
		for (size_t i = 0; i < 3; i++)
		{
			vxyz[i] -= out(i + 1);
		}
		clk += clk_corr / CLIGHT;
		dclk += dclk_corr / CLIGHT;
		return 1;
	}

	// clean internal function
	// ----------
	void t_gsatdata::_clear()
	{
		gtrace("t_gsatdat::_clear");

		t_gobsgnss::_clear();
		_ele = 0.0;
		_azi = 0.0;
		_rho = 0.0;
		_clk = 0.0;
	}



	// clean internal function
	// ----------
	bool t_gsatdata::_valid() const
	{
		gtrace("t_gsatdat::_valid");

		if (_rho == 0.0 ||
			t_gobsgnss::_valid()) return false;

		return true;
	}

	// reset eclipsing flag
	// ---------------------
	void  t_gsatdata::setecl(bool ecl)
	{
		_eclipse = ecl;
	}

	// ---------------------
	bool t_gsatdata::ecl()
	{
		gtrace("t_gsatdat::ecl");

		return _eclipse;
	}

	// -------------------------------------
	void t_gsatdata::addecl(map<string, t_gtime>& lastEcl)
	{
		gtrace("t_gsatdat::addecl");

		if (fabs(_b()) < EPS0_GPS && fabs(_orb_angle()) < EPS0_GPS) {
			_eclipse = true;
			lastEcl[_satid] = _epoch;
			return;
		}
		else {
			auto itLast = lastEcl.find(_satid);
			if (itLast != lastEcl.end()) {
				double tdiff = _epoch.diff(itLast->second);
				if (abs(tdiff) <= POST_SHADOW) _eclipse = true;
				else _eclipse = false;
			}
		}
	}

	// add postfit residuals
	//------------------
	void t_gsatdata::addres(RESIDTYPE restype, GOBSTYPE type, double res)
	{
		_gmutex.lock();

		if (restype == RESIDTYPE::RES_ORIG) {
			if (type == TYPE_C) _code_res_orig.push_back(res);
			if (type == TYPE_L) _phase_res_orig.push_back(res);
		}

		if (restype == RESIDTYPE::RES_NORM) {
			if (type == TYPE_C) _code_res_norm.push_back(res);
			if (type == TYPE_L) _phase_res_norm.push_back(res);
		}

		_gmutex.unlock();
	}


	// get postfit residuals
	//------------------
	vector<double> t_gsatdata::residuals(RESIDTYPE restype, GOBSTYPE type)
	{
		_gmutex.lock();

		vector<double> res;

		if (restype == RESIDTYPE::RES_ORIG) {
			if (type == TYPE_C) res = _code_res_orig;
			else if (type == TYPE_L) res = _phase_res_orig;
		}

		if (restype == RESIDTYPE::RES_NORM) {
			if (type == TYPE_C) res = _code_res_norm;
			else if (type == TYPE_L) res = _phase_res_norm;
		}

		_gmutex.unlock(); return res;
	}

	// clean residuals
	// ---------------   
	void t_gsatdata::clear_res(RESIDTYPE restype)
	{
		if (restype == RESIDTYPE::RES_ORIG) {
			_code_res_orig.clear();
			_phase_res_orig.clear();
		}

		if (restype == RESIDTYPE::RES_NORM) {
			_code_res_norm.clear();
			_phase_res_norm.clear();
		}
	}

	// add postfit residuals
	//------------------
	void t_gsatdata::addwind(const double& wind)
	{
		_gmutex.lock();

		_wind = wind;

		_gmutex.unlock();
	}


	// get stored wind up
	//------------------
	double t_gsatdata::wind()
	{
		return _wind;
	}

	// get hydrostatic mapping factor
	//------------------
	double t_gsatdata::mfH()
	{
		

		return _mfH;
	}

	// get wet mapping factor
	//------------------
	double t_gsatdata::mfW()
	{
		

		return _mfW;
	}

	// get tropo gradient mapping factor
	//------------------
	double t_gsatdata::mfG()
	{
		

		 return _mfG;
	}

	// Sun elevation relative to orbital plane
	double t_gsatdata::_b()
	{
		gtrace("t_gsatdat::beta");

		// test if already calculated
		if (!double_eq(_beta_val, 999)) return _beta_val;

		if (_satcrd.zero()) return 999;

		double beta = 0.0;
		double dt = 300;

		double dmjd = _epoch.dmjd();
		t_gephplan eph;
		ColumnVector Sun = eph.sunPos(dmjd, false).crd_cvect();   //ICRF
		double gmt = eph.gmst(dmjd);
		double gmt_dt = eph.gmst(dmjd + dt / 86400.0);

		ColumnVector Satcrd = _satcrd.crd_cvect();
		ColumnVector Satvel = _satvel.crd_cvect();
		ColumnVector Satcrd_dt = Satcrd + Satvel * dt;

		t_geop80 eop;
		Matrix prec = eop.precMatrix(dmjd);
		Matrix nut = eop.nutMatrix(dmjd);

		// ITRF -> ICRF
		Satcrd = prec.i() * nut.i() * rotZ(-gmt) * Satcrd;
		Satcrd_dt = prec.i() * nut.i() * rotZ(-gmt_dt) * Satcrd_dt;

		ColumnVector n = crossproduct(Satcrd, Satcrd_dt);

		n /= n.NormFrobenius();

		ColumnVector nSun = Sun / Sun.NormFrobenius();

		double cosa = DotProduct(nSun, n);

		beta = G_PI / 2.0 - acos(cosa);

		_beta_val = beta;

		return beta;
	}

	// Orbit angle
	double t_gsatdata::_orb_angle()
	{
		gtrace("t_gsatdat::orb_angle");

		if (!double_eq(_orb_angle_val, 999)) return _orb_angle_val;

		if (_satcrd.zero()) return 999;

		double mi = 0.0;
		double dt = 30;

		double dmjd = _epoch.dmjd();
		t_gephplan eph;
		ColumnVector Sun = eph.sunPos(dmjd, false).crd_cvect();
		double gmt = eph.gmst(dmjd);
		double gmt_dt = eph.gmst(dmjd + dt / 86400.0);

		ColumnVector Satcrd = _satcrd.crd_cvect();
		ColumnVector Satvel = _satvel.crd_cvect();
		ColumnVector Satcrd_dt = Satcrd + Satvel * dt;

		t_geop80 eop;
		Matrix prec = eop.precMatrix(dmjd);
		Matrix nut = eop.nutMatrix(dmjd);

		// ITRF -> ICRF
		Satcrd = prec.i() * nut.i() * rotZ(-gmt) * Satcrd;
		Satcrd_dt = prec.i() * nut.i() * rotZ(-gmt_dt) * Satcrd_dt;

		ColumnVector n = crossproduct(Satcrd, Satcrd_dt);

		ColumnVector es = Satcrd / Satcrd.NormFrobenius();

		ColumnVector eSun = Sun / Sun.NormFrobenius();

		ColumnVector p = crossproduct(Sun, n);
		p /= p.NormFrobenius();

		double E = acos(DotProduct(es, p));

		double SunSat = acos(DotProduct(es, eSun));

		if (SunSat > G_PI / 2) {
			if (E <= G_PI / 2) mi = G_PI / 2 - E;
			else            mi = G_PI / 2 - E;
		}
		else {
			if (E <= G_PI / 2) mi = G_PI / 2 + E;
			else            mi = E - G_PI - G_PI / 2;
		}

		_orb_angle_val = mi;

		return mi;
	}

}
