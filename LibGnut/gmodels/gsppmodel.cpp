
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.

-*/

#include "gmodels/gsppmodel.h"

//#define DEBUG_CORR
//#define DEBUG_GLO

namespace gnut {

	// Constructors
	// -------------------
	t_gsppmodel::t_gsppmodel()
	{
		gtrace("t_gsppmodel::t_gsppmodel");

		_tropoModel = make_shared<t_gtropo>();
		_gallbias = 0;
	}

	t_gsppmodel::t_gsppmodel(string site, t_glog* glog, t_gsetbase* settings)
		: _observ(OBSCOMBIN::IONO_FREE)
	{
		gtrace("t_gsppmodel::t_gsppmodel");

		_tropoModel = make_shared<t_gtropo>();
		_gallbias = 0;

		_settings = settings;
		_site = site;
		_log = glog;
		_phase = false;

		set<GSYS> systems = GNSS_SUPPORTED();
		for (set<GSYS>::iterator it = systems.begin(); it != systems.end(); it++) {
			_maxres_C[*it] = dynamic_cast<t_gsetgnss*>(_settings)->maxres_C(*it);
			_maxres_L[*it] = dynamic_cast<t_gsetgnss*>(_settings)->maxres_L(*it);
		}

		_maxres_norm = dynamic_cast<t_gsetproc*>(_settings)->max_res_norm();
		_tropo_mf    = dynamic_cast<t_gsetproc*>(_settings)->tropo_mf();
		_trpModStr   = dynamic_cast<t_gsetproc*>(_settings)->tropo_model();

		_resid_type = dynamic_cast<t_gsetproc*>(_settings)->residuals();
		_observ     = dynamic_cast<t_gsetproc*>(_settings)->obs_combin();
		_cbiaschar  = dynamic_cast<t_gsetproc*>(_settings)->cbiaschar();

		_band_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(_settings)->band_index(gnut::GPS);
		_band_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(_settings)->band_index(gnut::GAL);
		_band_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(_settings)->band_index(gnut::GLO);
		_band_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(_settings)->band_index(gnut::BDS);
		_band_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(_settings)->band_index(gnut::QZS);

		_freq_index[gnut::GPS] = dynamic_cast<t_gsetgnss*>(_settings)->freq_index(gnut::GPS);
		_freq_index[gnut::GAL] = dynamic_cast<t_gsetgnss*>(_settings)->freq_index(gnut::GAL);
		_freq_index[gnut::GLO] = dynamic_cast<t_gsetgnss*>(_settings)->freq_index(gnut::GLO);
		_freq_index[gnut::BDS] = dynamic_cast<t_gsetgnss*>(_settings)->freq_index(gnut::BDS);
		_freq_index[gnut::QZS] = dynamic_cast<t_gsetgnss*>(_settings)->freq_index(gnut::QZS);

		if (_trpModStr == TROPMODEL::SAASTAMOINEN)  _tropoModel = make_shared<t_saast>();
		else if (_trpModStr == TROPMODEL::DAVIS)         _tropoModel = make_shared<t_davis>();
		else if (_trpModStr == TROPMODEL::HOPFIELD)      _tropoModel = make_shared<t_hopf>();
	}

	// Destructor
	// ---------------------
	t_gsppmodel::~t_gsppmodel()
	{
		gtrace("t_gsppmodel::~t_gsppmodel");
	}


	// Outliers detection based on chi2 testing of normalized residuals
	// ----------
	int t_gsppmodel::outlierDetect_chi(vector<t_gsatdata>& data,
		SymmetricMatrix& Qx,
		const SymmetricMatrix& Qsav,
		const ColumnVector& v)
	{
		gtrace("t_gsppmodel::outlierDetect_chi");

		vector<t_gsatdata>::iterator it;
		vector<t_gsatdata>::iterator itMaxV;

		double maxV = 0.0;

		int ii = 0;
		for (it = data.begin(); it != data.end(); it++) {
			++ii;

			if (maxV == 0.0 || v(ii)*v(ii) > maxV) {

				maxV = v(ii)*v(ii);
				itMaxV = it;
			}

			if (_phase) {
				++ii;
				if (maxV == 0.0 || v(ii)*v(ii) > maxV) {

					maxV = v(ii)*v(ii);
					itMaxV = it;
				}
			}
		}

		if (maxV > 5.024) { // chi2(100) = 2.706; chi2(050) = 3.841; chi2(025) = 5.024; chi2(010) = 6.635; chi2(005) = 7.879;

			if (_log)
				_log->comment(2, "gpppmodel", _site + " outlier " + itMaxV->sat()
					+ " size:" + int2str(data.size())
					+ " v_norm: " + dbl2str(maxV)
					+ " " + itMaxV->epoch().str_hms());
			data.erase(itMaxV);
			Qx = Qsav;
			return 1;
		}

		return 0;
	}

	// Outliers detection
	// ----------
	int t_gsppmodel::outlierDetect(vector<t_gsatdata>& data,
		SymmetricMatrix& Qx,
		const SymmetricMatrix& Qsav,
		const ColumnVector& v)
	{
		gtrace("t_gsppmodel::outlierDetect");

		vector<t_gsatdata>::iterator itMaxVcodeGPS;
		vector<t_gsatdata>::iterator itMaxVcodeGLO;
		vector<t_gsatdata>::iterator itMaxVcodeGAL;
		vector<t_gsatdata>::iterator itMaxVcodeBDS;
		vector<t_gsatdata>::iterator itMaxVcodeQZS;

		vector<t_gsatdata>::iterator itMaxVphaseGPS;
		vector<t_gsatdata>::iterator itMaxVphaseGLO;
		vector<t_gsatdata>::iterator itMaxVphaseGAL;
		vector<t_gsatdata>::iterator itMaxVphaseBDS;
		vector<t_gsatdata>::iterator itMaxVphaseQZS;

		double maxVcodeGPS, maxVcodeGLO, maxVcodeGAL, maxVcodeBDS, maxVcodeQZS;
		maxVcodeGPS = maxVcodeGLO = maxVcodeGAL = maxVcodeBDS = maxVcodeQZS = 0.0;

		double maxVphaseGPS, maxVphaseGLO, maxVphaseGAL, maxVphaseBDS, maxVphaseQZS;
		maxVphaseGPS = maxVphaseGLO = maxVphaseGAL = maxVphaseBDS = maxVphaseQZS = 0.0;

		// deviding multi-freq residual vector into single-freq vectors
		vector<ColumnVector> v_frqs = _devide_res(v);

		// find maximal code/phase residuals for individual GNSS   
		for (unsigned int i = 0; i < v_frqs.size(); i++) 
		{
			double maxres = _maxres(v_frqs[i], GPS, false, data, itMaxVcodeGPS);
			if (maxres > maxVcodeGPS) maxVcodeGPS = maxres;
			maxres = _maxres(v_frqs[i], GLO, false, data, itMaxVcodeGLO);
			if (maxres > maxVcodeGLO) maxVcodeGLO = maxres;
			maxres = _maxres(v_frqs[i], GAL, false, data, itMaxVcodeGAL);
			if (maxres > maxVcodeGAL) maxVcodeGAL = maxres;
			maxres = _maxres(v_frqs[i], BDS, false, data, itMaxVcodeBDS);
			if (maxres > maxVcodeBDS) maxVcodeBDS = maxres;
			maxres = _maxres(v_frqs[i], QZS, false, data, itMaxVcodeQZS);
			if (maxres > maxVcodeQZS) maxVcodeQZS = maxres;

			maxres = _maxres(v_frqs[i], GPS, true, data, itMaxVphaseGPS);
			if (maxres > maxVphaseGPS) maxVphaseGPS = maxres;
			maxres = _maxres(v_frqs[i], GLO, true, data, itMaxVphaseGLO);
			if (maxres > maxVphaseGLO) maxVphaseGLO = maxres;
			maxres = _maxres(v_frqs[i], GAL, true, data, itMaxVphaseGAL);
			if (maxres > maxVphaseGAL) maxVphaseGAL = maxres;
			maxres = _maxres(v_frqs[i], BDS, true, data, itMaxVphaseBDS);
			if (maxres > maxVphaseBDS) maxVphaseBDS = maxres;
			maxres = _maxres(v_frqs[i], QZS, true, data, itMaxVphaseQZS);
			if (maxres > maxVphaseQZS) maxVphaseQZS = maxres;
		}

#ifdef DEBUG
		cout << "Max res range: " << maxVcodeGPS << " " << itMaxVcodeGPS->sat() << endl;
		cout << "Max res phase: " << maxVphaseGPS << " " << itMaxVphaseGPS->sat() << endl;
#endif
		//Only detect the outliers for the constellations with maximum outliers
		//The GLONASS is set as the basic reference
		double maxvc = maxVcodeGLO, maxvp = maxVphaseGLO;
		GSYS   maxsys = GLO;
		if (maxVcodeGPS > maxvc || maxVphaseGPS > maxvp) { maxvc = maxVcodeGPS; maxvp = maxVphaseGPS; maxsys = GPS; }
		if (maxVcodeGAL > maxvc || maxVcodeGAL > maxvp) { maxvc = maxVcodeGAL; maxvp = maxVcodeGAL;  maxsys = GAL; }
		if (maxVcodeBDS > maxvc || maxVphaseBDS > maxvp) { maxvc = maxVcodeBDS; maxvp = maxVphaseBDS; maxsys = BDS; }
		if (maxVcodeQZS > maxvc || maxVphaseQZS > maxvp) { maxvc = maxVcodeQZS; maxvp = maxVphaseQZS; maxsys = QZS; }


		if (maxsys == GLO && _check_outl(false, maxVcodeGLO, itMaxVcodeGLO, data)) { data.erase(itMaxVcodeGLO);  Qx = Qsav; return 1; }
		if (maxsys == GLO && _check_outl(true, maxVphaseGLO, itMaxVphaseGLO, data)) { data.erase(itMaxVphaseGLO); Qx = Qsav; return 1; }

		if (maxsys == GPS && _check_outl(false, maxVcodeGPS, itMaxVcodeGPS, data)) { data.erase(itMaxVcodeGPS);  Qx = Qsav; return 1; }
		if (maxsys == GPS && _check_outl(true, maxVphaseGPS, itMaxVphaseGPS, data)) { data.erase(itMaxVphaseGPS); Qx = Qsav; return 1; }

		if (maxsys == GAL && _check_outl(false, maxVcodeGAL, itMaxVcodeGAL, data)) { data.erase(itMaxVcodeGAL);  Qx = Qsav; return 1; }
		if (maxsys == GAL && _check_outl(true, maxVphaseGAL, itMaxVphaseGAL, data)) { data.erase(itMaxVphaseGAL); Qx = Qsav; return 1; }

		if (maxsys == BDS && _check_outl(false, maxVcodeBDS, itMaxVcodeBDS, data)) { data.erase(itMaxVcodeBDS);  Qx = Qsav; return 1; }
		if (maxsys == BDS && _check_outl(true, maxVphaseBDS, itMaxVphaseBDS, data)) { data.erase(itMaxVphaseBDS); Qx = Qsav; return 1; }

		if (maxsys == QZS && _check_outl(false, maxVcodeQZS, itMaxVcodeQZS, data)) { data.erase(itMaxVcodeQZS);  Qx = Qsav; return 1; }
		if (maxsys == QZS && _check_outl(true, maxVphaseQZS, itMaxVphaseQZS, data)) { data.erase(itMaxVphaseQZS); Qx = Qsav; return 1; }

		return 0;
	}


	// Outliers detection
	// ----------
	int t_gsppmodel::outlierDetect(vector<t_gsatdata>& data,
		SymmetricMatrix& Qx,
		const SymmetricMatrix& Qsav)
	{
		gtrace("t_gsppmodel::outlierDetect");

		vector<t_gsatdata>::iterator itMaxVcodeNORM = data.end();
		vector<t_gsatdata>::iterator itMaxVphaseNORM = data.end();

		vector<t_gsatdata>::iterator itMaxVcodeORIG = data.end();
		vector<t_gsatdata>::iterator itMaxVphaseORIG = data.end();

		vector<t_gsatdata>::iterator itDataErase = data.end();

		double maxVcodeNORM = 0.0;
		double maxVphaseNORM = 0.0;

		double maxVcodeORIG = 0.0;
		double maxVphaseORIG = 0.0;

		// find maximal code/phase residuals
		maxVcodeNORM = _maxres(false, data, itMaxVcodeNORM, RESIDTYPE::RES_NORM);
		maxVphaseNORM = _maxres(true, data, itMaxVphaseNORM, RESIDTYPE::RES_NORM);

		maxVcodeORIG = _maxres(false, data, itMaxVcodeORIG, RESIDTYPE::RES_ORIG);
		maxVphaseORIG = _maxres(true, data, itMaxVphaseORIG, RESIDTYPE::RES_ORIG);

#ifdef DEBUG
		if (itMaxVcodeNORM != data.end()) cout << "Max res range norm: " << fixed << setprecision(3) << maxVcodeNORM << " " << itMaxVcodeNORM->sat() << endl;
		if (itMaxVphaseNORM != data.end()) cout << "Max res phase norm: " << fixed << setprecision(3) << maxVphaseNORM << " " << itMaxVphaseNORM->sat() << endl;
		if (itMaxVcodeORIG != data.end()) cout << "Max res range orig: " << fixed << setprecision(3) << maxVcodeORIG << " " << itMaxVcodeORIG->sat() << endl;
		if (itMaxVphaseORIG != data.end()) cout << "Max res phase orig: " << fixed << setprecision(3) << maxVphaseORIG << " " << itMaxVphaseORIG->sat() << endl;
		int ooo; cin >> ooo;
#endif

		if (_check_outl(true, maxVphaseNORM, itMaxVphaseNORM, maxVphaseORIG, itMaxVphaseORIG, itDataErase, data)) { 
			auto it = find(_outlier_sat.begin(), _outlier_sat.end(), itDataErase->sat());
			if (it != _outlier_sat.end()) {
				_outlier_sat.push_back(itDataErase->sat());
			}
			data.erase(itDataErase);
			Qx = Qsav; return 1; 
		}
		if (_check_outl(false, maxVcodeNORM, itMaxVcodeNORM, maxVcodeORIG, itMaxVcodeORIG, itDataErase, data)) { 
			auto it = find(_outlier_sat.begin(), _outlier_sat.end(), itDataErase->sat());
			if (it == _outlier_sat.end()) {
				_outlier_sat.push_back(itDataErase->sat());
			}
			data.erase(itDataErase);
			Qx = Qsav; return 1; 
		}

		return 0;
	}

	int t_gsppmodel::outlierDetect(vector<t_gsatdata>& data, SymmetricMatrix & Qx, const SymmetricMatrix &Qsav, vector<t_gsatdata>::iterator & exception_sat)
	{
		gtrace("t_gsppmodel::outlierDetect");

		vector<t_gsatdata>::iterator itMaxVcodeNORM = data.end();
		vector<t_gsatdata>::iterator itMaxVphaseNORM = data.end();

		vector<t_gsatdata>::iterator itMaxVcodeORIG = data.end();
		vector<t_gsatdata>::iterator itMaxVphaseORIG = data.end();

		vector<t_gsatdata>::iterator itDataErase = data.end();

		double maxVcodeNORM = 0.0;
		double maxVphaseNORM = 0.0;

		double maxVcodeORIG = 0.0;
		double maxVphaseORIG = 0.0;

		// find maximal code/phase residuals
		maxVcodeNORM = _maxres(false, data, itMaxVcodeNORM, RESIDTYPE::RES_NORM);
		maxVphaseNORM = _maxres(true, data, itMaxVphaseNORM, RESIDTYPE::RES_NORM);

		maxVcodeORIG = _maxres(false, data, itMaxVcodeORIG, RESIDTYPE::RES_ORIG);
		maxVphaseORIG = _maxres(true, data, itMaxVphaseORIG, RESIDTYPE::RES_ORIG);

		if (_check_outl(true, maxVphaseNORM, itMaxVphaseNORM, maxVphaseORIG, itMaxVphaseORIG, itDataErase, data))
		{
			exception_sat = itDataErase;
			Qx = Qsav;
			return 1;
		}
		if (_check_outl(false, maxVcodeNORM, itMaxVcodeNORM, maxVcodeORIG, itMaxVcodeORIG, itDataErase, data))
		{
			exception_sat = itDataErase;
			Qx = Qsav;
			return 1;
		}
		return 0;
	}

	// model for computed range value 
	// ----------
	double t_gsppmodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& gsatdata, t_gobs& gobs, bool com)
	{
		gtrace("t_gsppmodel::cmpObs");

		// Cartesian coordinates to ellipsodial coordinates
		t_gtriple xyz, ell;

		string strEst = dynamic_cast<t_gsetgen*>(_settings)->estimator();
		bool isFLT = (strEst == "FLT");
		if (isFLT)
		{
			if (param.getCrdParam(_site, xyz) < 0)
			{
				xyz = _grec->crd_arp(epoch);
			}
		}
		else {
			if (param.getCrdParam(_site, xyz) < 0)
			{
				xyz = _grec->crd(epoch);
			}
			xyz += _grec->eccxyz(epoch);
		}

		xyz2ell(xyz, ell, false);

		t_gtriple satcrd = gsatdata.satcrd();
		ColumnVector cSat = satcrd.crd_cvect();

		string sat = gsatdata.sat();
		string rec = gsatdata.site();
		t_gtime epo = gsatdata.epoch();

		// Tropospheric wet delay correction
		double trpDelay = 0;
		trpDelay = tropoDelay(epoch, param, ell, gsatdata);
		if (fabs(trpDelay) > 50)
		{
			if(_log) _log->logError("t_gsppmodel", "cmpObs", "trpDelay > 50");
			return -1;
		}

		// idx for par correction
		int i = -1;

		// Receiver clock correction 
		double clkRec = 0.0;
		i = param.getParam(_site, par_type::CLK, "");
		if (i >= 0)
		{
			clkRec = param[i].value();
		}
		else
		{
			if (_log) _log->logError("t_gsppmodel", "cmpObs", _site + " ! warning:  Receiver Clock is not included in parameters!");
		}
		
		// system time offset	
		double isb_offset = isbCorrection(param, sat , rec , gobs);

		double ifb = 0.0;
		i = param.getParam(_site, par_type::IFB_GPS, "");
		if (i >= 0 && gobs.is_code() &&_freq_index[gsatdata.gsys()][ gobs.band()] == FREQ_3) {
			ifb = param[i].value();
		}

		i = param.getParam(_site, par_type::IFB_GAL, "");
		if (i >= 0 && gobs.is_code() && _freq_index[gsatdata.gsys()][gobs.band()] == FREQ_3) {
			ifb = param[i].value();
		}
		i = param.getParam(_site, par_type::IFB_BDS, "");
		if (i >= 0 && gobs.is_code() && _freq_index[gsatdata.gsys()][gobs.band()] == FREQ_3) {
			ifb = param[i].value();
		}
		i = param.getParam(_site, par_type::IFB_QZS, "");
		if (i >= 0 && gobs.is_code() && _freq_index[gsatdata.gsys()][gobs.band()] == FREQ_3) {
			ifb = param[i].value();
		}

		// Return value
		return gsatdata.rho() +
			clkRec -
			gsatdata.clk() +
			trpDelay +
			isb_offset + 
            ifb;
	}
	
	double t_gsppmodel::cmpObsD(t_gtime & epoch, t_gallpar & param, t_gsatdata& gsatdata, t_gobs & gobs)
	{
		gtrace("t_gsppmodel::cmpObsD");

		// Cartesian coordinates to ellipsodial coordinates
		t_gtriple xyz, ell;
		ColumnVector cRec(3), vRec(3);

		if (param.getCrdParam(_site, xyz) > 0) {
			cRec = xyz.crd_cvect();
		}
		else {
			xyz = _grec->crd_arp(epoch);
			cRec = xyz.crd_cvect();
		}
		xyz2ell(xyz, ell, false);

		int i = param.getParam(_site, par_type::VEL_X, "");
		int j = param.getParam(_site, par_type::VEL_Y, "");
		int k = param.getParam(_site, par_type::VEL_Z, "");
		int l = param.getParam(_site, par_type::CLK_RAT, "");
		vRec(1) = param[i].value();
		vRec(2) = param[j].value();
		vRec(3) = param[k].value();
		double dclk_Rec = param[l].value();
		double dclk = gsatdata.dclk();
		t_gtriple satcrd = gsatdata.satcrd();
		t_gtriple satvel = gsatdata.satvel();
		ColumnVector cSat = satcrd.crd_cvect();
		ColumnVector vSat = satvel.crd_cvect();
		ColumnVector e = (cSat - cRec) / gsatdata.rho();


		double res = dotproduct(e, vSat - vRec) +
			OMEGA / CLIGHT * (vSat(2)*cRec(1) + cSat(2)*vRec(1)
				- vSat(1)*cRec(2) - cSat(1)*vRec(2));


		return res +
			dclk_Rec -
			dclk;
	}

	// Compute troposperic delay
	// -----------
	double t_gsppmodel::tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple ell, t_gsatdata& satdata)
	{
		gtrace("t_gsppmodel::tropoDelay");

		if (_tropoModel == 0)
		{
			if (_log)
			{
				_log->comment(0, "gppp", "Tropo Model setting is not correct. Default used! Check config.");
			}
			else
			{
				cerr << "gppp - Tropo Model setting is not correct. Default used! Check config.\n";
			}
			_tropoModel = make_shared<t_saast>();
		}

		double ele = satdata.ele();

		double delay = 0.0;
		double zwd = 0.0;
		double zhd = 0.0;

		if (abs(ell[2]) > 1E4)
		{
			return 0.0;
		}

		int i;
		i = param.getParam(_site, par_type::TRP, "");
		if (i >= 0)
		{
			zwd = param[i].value();
			zhd = param[i].apriori();
		}
		else
		{
			if (_tropoModel != 0)
			{
				zwd = _tropoModel->getZWD(ell, epoch);
				zhd = _tropoModel->getZHD(ell, epoch);
			}
		}

		if (_tropo_mf == ZTDMPFUNC::GMF)
		{
			double gmfh, gmfw, dgmfh, dgmfw;
			t_gmf mf;
			mf.gmf(epoch.mjd(), ell[0], ell[1], ell[2], G_PI / 2.0 - ele,
				gmfh, gmfw, dgmfh, dgmfw);
			//      delay = gmfh * _tropoModel->getZHD(ell, epoch) + gmfw * zwd;
			delay = gmfh * zhd + gmfw * zwd;

		}
		else if (_tropo_mf == ZTDMPFUNC::COSZ)
		{
			double mf = 1 / sin(ele);
			//      delay = 1/sin(ele) * _tropoModel->getZHD(ell, epoch) +1/sin(ele) * zwd;
			delay = mf * zhd + mf * zwd;
		}

		return delay;
	}

	// Compute ionospheic delay
	// -----------
	double t_gsppmodel::ionoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& gsatdata, t_gobs& gobs)
	{
		gtrace("t_gsppmodel::ionoDelay");

		double iono_param = 0.0;
		double iono_apriori = 0.0;
		//double iono_model   = 0.0;
		double iono_delay = 0.0;

		double mf = 1 / sqrt(1.0 - pow(R_SPHERE / (R_SPHERE + 450000.0) * sin(G_PI / 2.0 - gsatdata.ele()), 2));

		double f1 = gsatdata.frequency(t_gsys::band_priority(gsatdata.gsys(), FREQ_1));
		double fk = gsatdata.frequency(gobs.band());
		double alfa = 0.0;

		if (gobs.is_phase())
		{
			alfa = -(f1 * f1) / (fk * fk);
		}
		if (gobs.is_code())
		{
			alfa = (f1 * f1) / (fk * fk);
		}

		// ionosphere slant delay parameter
		int i = param.getParam(_site, par_type::SION, gsatdata.sat());
		if (i >= 0)
		{
			iono_delay = alfa * param[i].value();
		}

		// ionosphere vertical delay paremeter
		i = param.getParam(_site, par_type::VION, gsatdata.sat());
		if (i >= 0)
		{
			iono_apriori = alfa * mf * param[i].apriori();
			iono_param = alfa * mf * param[i].value();
			iono_delay = iono_apriori + iono_param;
		}

		return iono_delay;
	}

	double t_gsppmodel::isbCorrection(t_gallpar & param, string & sat, string & rec, t_gobs & gobs)
	{
		double isb_offset = 0.0;

		auto gsys = t_gsys::sat2gsys(sat);

		switch (gsys)
		{
		case GPS:
		{
			break;
		}
		// GLONASS system time offset
		case GLO:
		{
			int idx_isb = param.getParam(rec, par_type::GLO_ISB, "");
			if (idx_isb >= 0)
			{
				isb_offset = param[idx_isb].value();
				//cout << "GLO Offset: " << sat << " " << glonass_offset << " " << gsatdata.epoch().str_hms() << endl;
			}

			int idx_ifb = param.getParam(rec, par_type::GLO_IFB, sat);
			if (idx_ifb >= 0)
			{
				isb_offset += param[idx_ifb].value();
			}

			break;
		}
		case GAL:
		{
			// Galileo system time offset
			int i = param.getParam(rec, par_type::GAL_ISB, "");
			if (i >= 0)
			{
				isb_offset = param[i].value();
				//cout << "GAL Offset: " << sat << " " << galileo_offset << " " << gsatdata.epoch().str_hms() << endl;       
			}
			break;
		}
		case BDS:
		{
			// BaiDou system time offset
			int i = param.getParam(rec, par_type::BDS_ISB, "");
			if (i >= 0)
			{
				isb_offset = param[i].value();
				//cout << "BDS ISB: " << sat << " " << param[i].value() << " " << gsatdata.epoch().str_hms() << endl;
			}
			break;
		}
		// QZSS system time offset
		case QZS:
		{
			int i = param.getParam(_site, par_type::QZS_ISB, "");
			if (i >= 0)
			{
				isb_offset = param[i].value();
			}
			break;
		}
		default:
			throw logic_error("can not support such sys : " + t_gsys::gsys2str(gsys));
		}

		return isb_offset;
	}

	void t_gsppmodel::reset_observ(OBSCOMBIN observ)
	{
		_observ = observ;
	}

	void t_gsppmodel::setrec(shared_ptr<t_gobj> rec)
	{
		_grec = rec;
	}

	// Find maximal residual
	double t_gsppmodel::_maxres(const ColumnVector& v, GSYS gs, bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA)
	{
		unsigned int inc = 2;

		if (v.Nrows() == data.size()) {      // code data only 
			if (phase) return 0.0;
			inc = 1;
		}

		vector<t_gsatdata>::iterator it;
		int ii = 1;
		if (phase) ii = 2;
		double maxres = 0.0;

		for (it = data.begin(); it != data.end(); it++) {
			if (it->gsys() != gs) {
				ii += inc;
				continue;
			}

			if (maxres == 0.0 || abs(v(ii)) > maxres) {
				maxres = abs(v(ii));
				itDATA = it;
			}

			ii += inc;
		}

		return maxres;
	}


	// Find maximal residual
	double t_gsppmodel::_maxres(bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA, RESIDTYPE res_type, GSYS gs)
	{

		vector<t_gsatdata>::iterator it;

		double maxres = 0.0;

		for (it = data.begin(); it != data.end(); it++) {
			if (it->gsys() != gs && gs != GNS) continue;

			vector<double> res;
			if (phase) res = it->residuals(res_type, TYPE_L);
			else      res = it->residuals(res_type, TYPE_C);

			for (auto itRES = res.begin(); itRES != res.end(); itRES++) {
				if (maxres == 0.0 || fabs(*itRES) > maxres) {
					maxres = fabs(*itRES);
					itDATA = it;
				}
			}

		}

		return maxres;
	}

	// check maximal residual   
	bool t_gsppmodel::_check_outl(bool phase, double& maxres, vector<t_gsatdata>::iterator& itData, vector<t_gsatdata>& data)
	{
		map<GSYS, double> map_res;
		if (phase) map_res = _maxres_L;
		else      map_res = _maxres_C;

		GSYS gs;
		if (itData != data.end()) gs = itData->gsys();
		else return false;

		if ((maxres > map_res[gs] && _resid_type == RESIDTYPE::RES_ORIG) ||
			(maxres > _maxres_norm && _resid_type == RESIDTYPE::RES_NORM)) {

			if (phase) _logOutl(true, itData->sat(), data.size(), maxres, itData->ele_deg(), itData->epoch(), _resid_type);
			else      _logOutl(true, itData->sat(), data.size(), maxres, itData->ele_deg(), itData->epoch(), _resid_type);

			return true;
		}
		return false;
	}

	// check maximal residual   
	bool t_gsppmodel::_check_outl(bool phase, double& maxresNORM, vector<t_gsatdata>::iterator& itDataNORM,
		double& maxresORIG, vector<t_gsatdata>::iterator& itDataORIG,
		vector<t_gsatdata>::iterator& itDataErase, vector<t_gsatdata>& data)
	{
		map<GSYS, double> map_res;
		if (phase) map_res = _maxres_L;
		else      map_res = _maxres_C;

		GSYS gs;
		if (itDataORIG != data.end()) gs = itDataORIG->gsys();
		else return false;

		if (_resid_type == RESIDTYPE::RES_ORIG)
		{
			if (maxresORIG > map_res[gs])
			{
				itDataErase = itDataORIG;
				if (phase) _logOutl(true, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), _resid_type);
				else      _logOutl(false, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), _resid_type);
				return true;
			}
		}
		else if (_resid_type == RESIDTYPE::RES_NORM)
		{
			if (maxresNORM > _maxres_norm)
			{
				itDataErase = itDataNORM;
				if (phase) _logOutl(true, itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), _resid_type);
				else      _logOutl(false, itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), _resid_type);
				return true;
			}
		}
		else if (_resid_type == RESIDTYPE::RES_ALL)
		{
			if (maxresORIG > map_res[gs] || maxresNORM > _maxres_norm)
			{
				if (itDataORIG == itDataNORM)
				{
					itDataErase = itDataORIG;
					if (phase) _logOutl(true, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RESIDTYPE::RES_ORIG);
					else      _logOutl(false, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RESIDTYPE::RES_ORIG);
					return true;
				}
				else
				{
					double ratioORIG = maxresORIG / map_res[gs];
					double ratioNORM = maxresNORM / _maxres_norm;
					if (ratioNORM >= ratioORIG)
					{
						itDataErase = itDataNORM;
						if (phase) _logOutl(true, itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), RESIDTYPE::RES_NORM);
						else      _logOutl(false, itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), RESIDTYPE::RES_NORM);
						return true;
					}
					else
					{
						itDataErase = itDataORIG;
						if (phase) _logOutl(true, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RESIDTYPE::RES_ORIG);
						else      _logOutl(false, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RESIDTYPE::RES_ORIG);
						return true;
					}
				}
			}
		}

		return false;
	}

	// logging outlier   
	void t_gsppmodel::_logOutl(bool phase, string prn, int data_size, double maxres, double ele, t_gtime epo, RESIDTYPE resid_type)
	{
		string obsType = "";
		string resType = "";
		if (phase) obsType = "phase";
		else      obsType = "range";

		if (resid_type == RESIDTYPE::RES_NORM) resType = "Norm residual";
		if (resid_type == RESIDTYPE::RES_ORIG) resType = "Orig residual";
		if (resid_type == RESIDTYPE::RES_ALL)  resType = "All residual";

		ostringstream os;
		os << _site << " outlier (" << resType << ": " << obsType << ") " << prn
			<< " size:" << fixed << setw(2) << data_size
			<< " v: " << fixed << setw(16) << right << setprecision(3) << maxres
			<< " ele: " << fixed << setw(6) << setprecision(2) << ele
			<< " " << epo.str_hms();
		if (_log)
			_log->comment(1, "gppp", os.str());

	}

	// devide multi-freq residuals vector into single-freq
	vector<ColumnVector> t_gsppmodel::_devide_res(const ColumnVector& v_orig)
	{
		vector<ColumnVector> vec;

		unsigned int k = 1;
		if (_observ == OBSCOMBIN::RAW_DOUBLE) k = 2;

		if (k == 1) {
			vec.push_back(v_orig);
			return vec;
		}

		ColumnVector v_L1(v_orig.Nrows() / k);
		ColumnVector v_L2(v_orig.Nrows() / k);

		int i = 1;
		int j = 1;
		while (i <= v_orig.Nrows() - 1)
		{
			v_L1(j) = v_orig(i);
			v_L2(j) = v_orig(i + 1);
			j += 1;
			i += k;
		}

		vec.push_back(v_L1);
		vec.push_back(v_L2);

		return vec;
	}


} // namespace
