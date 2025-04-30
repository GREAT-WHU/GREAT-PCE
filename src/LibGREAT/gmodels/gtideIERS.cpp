/**
 * @file         gtideIERS.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        precise model for computing correction
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gmodels/gtideIERS.h"
#include "gutils/gconst.h"
#include "gutils/gsysconv.h"


namespace great
{
	t_gtideIERS::t_gtideIERS(t_glog* l) : t_gtide(l)
	{

	}

	t_gtriple t_gtideIERS::tide_solid(const t_gtime & epo, t_gtriple & xyz, Matrix& rot_trs2crs, t_gnavde* nav_planet)
	{
		const double h20 = 0.6026e0;
		const double dh2 = -0.0006e0;
		const double l20 = 0.0831e0;
		const double dl2 = 0.0002e0;
		const double h3 = 0.293e0;
		const double l3 = 0.015e0;


		ColumnVector sun_pos(0.0,3),moon_pos(0.0,3);

		// get sun pos in J2000
		t_gtime epo_tt = epo;
		epo_tt.tsys(t_gtime::TT);

		nav_planet->get_pos(epo_tt.dmjd(), "SUN", sun_pos);
		nav_planet->get_pos(epo_tt.dmjd(), "MOON", moon_pos);

		// change to TRS, unit: km
		sun_pos = rot_trs2crs.t() * sun_pos;
		//t_gtriple sun_pos_trs(sun_pos);

		moon_pos = rot_trs2crs.t() * moon_pos;
		//t_gtriple moon_pos_trs(moon_pos);

		double rSun = sqrt(DotProduct(sun_pos, sun_pos));
		double rMoon = sqrt(DotProduct(moon_pos, moon_pos));

		// unit vector for site
		ColumnVector unit_site_pos = xyz.unitary();
		// get BLH
		t_gtriple blh;
		xyz2ell(xyz, blh, false);
		double lat = blh[0];
		double lon = blh[1];
		double height = blh[2];
		double colat = G_PI / 2 - blh[0];

		// calculate
		ColumnVector dx(0.0,3);
		for (int i = 0; i < 2;i++)
		{
			ColumnVector unit_pos(0.0,3);
			double GM_body = 0.0, R_body,dotl = 0.0, scale = 0.0;
			switch (i)
			{
			case 0:
				unit_pos = sun_pos / rSun;
				dotl = DotProduct(unit_pos, unit_site_pos);
				GM_body = SUN_GM;
				R_body = rSun;
				break;
			case 1:
				unit_pos = moon_pos / rMoon;
				dotl = DotProduct(unit_pos, unit_site_pos);
				GM_body = MOON_GM;
				R_body = rMoon;
				break;
			}
			
			scale = GM_body / EARTH_GM * EARTH_R * pow((EARTH_R / R_body), 3);

			double h2 = h20 + dh2 * (3.0 * sin(lat) * sin(lat) - 1.0) / 2.0;
			double l2 = l20 + dl2 * (3.0 * sin(lat) * sin(lat) - 1.0) / 2.0;
			double usite_part = (h2 * (1.5 * dotl * dotl - 0.50) - 3.0 * l2 * dotl * dotl) * scale;
			double ubody_part = 3.e0 * l2 * dotl * scale;

			dx = dx + usite_part * unit_site_pos + ubody_part * unit_pos;

			scale = scale * (EARTH_R / R_body);
			ubody_part = (7.5 * dotl * dotl - 1.5)*l3;
			usite_part = (h3 * ( 2.5 * dotl * dotl - 1.5) * dotl - ubody_part * dotl)*scale;
			ubody_part = ubody_part * scale;
			dx = dx + usite_part * unit_site_pos + ubody_part * unit_pos;
		}
		return t_gtriple(dx);
	}

	t_gtriple t_gtideIERS::tide_solid(const t_gtime& epo, t_gtriple& xyz, Matrix& rot_trs2crs, ColumnVector& sun_pos, ColumnVector& moon_pos)
	{
		const double h20 = 0.6026e0;
		const double dh2 = -0.0006e0;
		const double l20 = 0.0831e0;
		const double dl2 = 0.0002e0;
		const double h3 = 0.293e0;
		const double l3 = 0.015e0;

		// get sun pos in J2000
		t_gtime epo_tt = epo;
		epo_tt.tsys(t_gtime::TT);

		// change to TRS, unit: km
		sun_pos = rot_trs2crs.t() * sun_pos;
		//t_gtriple sun_pos_trs(sun_pos);

		moon_pos = rot_trs2crs.t() * moon_pos;
		//t_gtriple moon_pos_trs(moon_pos);

		double rSun = sqrt(DotProduct(sun_pos, sun_pos));
		double rMoon = sqrt(DotProduct(moon_pos, moon_pos));

		// unit vector for site
		ColumnVector unit_site_pos = xyz.unitary();
		// get BLH
		t_gtriple blh;
		xyz2ell(xyz, blh, false);
		double lat = blh[0];
		double lon = blh[1];
		double height = blh[2];
		double colat = G_PI / 2 - blh[0];

		// calculate
		ColumnVector dx(0.0, 3);
		for (int i = 0; i < 2; i++)
		{
			ColumnVector unit_pos(0.0, 3);
			double GM_body = 0.0, R_body, dotl = 0.0, scale = 0.0;
			switch (i)
			{
			case 0:
				unit_pos = sun_pos / rSun;
				dotl = DotProduct(unit_pos, unit_site_pos);
				GM_body = SUN_GM;
				R_body = rSun;
				break;
			case 1:
				unit_pos = moon_pos / rMoon;
				dotl = DotProduct(unit_pos, unit_site_pos);
				GM_body = MOON_GM;
				R_body = rMoon;
				break;
			}

			scale = GM_body / EARTH_GM * EARTH_R * pow((EARTH_R / R_body), 3);

			double h2 = h20 + dh2 * (3.0 * sin(lat) * sin(lat) - 1.0) / 2.0;
			double l2 = l20 + dl2 * (3.0 * sin(lat) * sin(lat) - 1.0) / 2.0;
			double usite_part = (h2 * (1.5 * dotl * dotl - 0.50) - 3.0 * l2 * dotl * dotl) * scale;
			double ubody_part = 3.e0 * l2 * dotl * scale;

			dx = dx + usite_part * unit_site_pos + ubody_part * unit_pos;

			scale = scale * (EARTH_R / R_body);
			ubody_part = (7.5 * dotl * dotl - 1.5) * l3;
			usite_part = (h3 * (2.5 * dotl * dotl - 1.5) * dotl - ubody_part * dotl) * scale;
			ubody_part = ubody_part * scale;
			dx = dx + usite_part * unit_site_pos + ubody_part * unit_pos;
		}
		return t_gtriple(dx);
	}

	t_gtriple t_gtideIERS::tide_pole()
	{
		return t_gtriple();
	}

	t_gtriple t_gtideIERS::tide_pole_pod(const t_gtime& epo, double xpole, double ypole, t_gtriple & xyz)
	{
		t_gtriple blh;
		xyz2ell(xyz, blh, false);
		double lat = blh[0];
		double lon = blh[1];
		double height = blh[2];
		double colat = G_PI / 2 - blh[0];

		double mean_xpole = 0.0, mean_ypole = 0.0;
		getMeanPole(epo.dmjd(), mean_xpole, mean_ypole);

		mean_xpole = xpole - mean_xpole;
		mean_ypole = -(ypole - mean_ypole);

		t_gtriple dxi(0.0,0.0,0.0);

		dxi[1] = 9.0*cos(colat)     *(mean_xpole*sin(lon) - mean_ypole * cos(lon));
		dxi[0] = 9.0*cos(2.0*colat)*(mean_xpole*cos(lon) + mean_ypole * sin(lon));
		dxi[2] = -32.0*sin(2.0*colat)*(mean_xpole*cos(lon) + mean_ypole * sin(lon));
		
		t_gtriple dx(0.0, 0.0, 0.0);
		//enu neu
		neu2xyz(blh, dxi, dx);

		//!!rotation matrix from east - north - radial to x - y - z
		//dx = dx + dxi / 1.e6;
		dx = dx / 1.e6;
		return dx;
	}
	t_gtriple t_gtideIERS::load_ocean(const t_gtime & epoch, const string & site, const t_gtriple & xRec)
	{
		// calcute angular argument
		const double dtr = 0.174532925199e-1;
		const double speed[11] = { 1.40519e-4,1.45444e-4,1.37880e-4,1.45842e-4,0.72921e-4,0.67598e-4 ,
								   0.72523e-4,0.64959e-4,0.053234e-4,0.026392e-4,0.003982e-4 };
		
		const double angfac[4][11] = { {2.0,0.0,2.0,2.0,1.0,1.0,-1.0,1.0,0.0,0.0,2.0},
									   {-2.0,0.0,-3.0,0.0,0.0,-2.0,0.0,-3.0,2.0,1.0,0.0},
									   {0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0},
									   {0.0,0.0,0.0,0.0,0.25,-0.25,-0.25,-0.25,0.0,0.0,0.0} };

	
		if (!_gotl) {//lvhb modified in 2020924
			if (_log) {
				_log->comment(1, "No ocean load file!");
			}
			return t_gtriple();
		}
		int capd = epoch.doy() + 365 * (epoch.year() - 1975) + (epoch.year() - 1973) / 4;
		double capt = (27392.500528e0 + 1.000000035e0 * capd) / 36525.0;
		double h0 = (279.69668e0 + (36000.768930485e0 + 3.03e-4 * capt)*capt)*dtr;
		double s0 = (((1.9e-6 * capt - 0.001133)*capt + 481267.88314137)*capt + 270.434358)*dtr;
		double p0 = (((-1.2e-5 * capt - 0.010325)*capt + 4069.0340329577)*capt + 334.329653)*dtr;

		double angle[11] = { 0.0 };
		double sod = epoch.sod() + epoch.dsec();
		for (int i = 0; i < 11; i++)
		{
			angle[i] = speed[i] * sod + angfac[0][i] * h0 + angfac[1][i] * s0 + angfac[2][i] * p0 + angfac[3][i] * G_PI * 2;
			angle[i] = fmod(angle[i], G_PI * 2);
			if (angle[i] < 0) angle[i] = angle[i] + G_PI * 2;
		}
		
		//get oceanload data
		t_gtriple ell;
		xyz2ell(xRec, ell, true); 
		Matrix coef_oceanload;
		int isexist = _gotl->data(coef_oceanload, ell[1],ell[0]);
		if (isexist < 0)
		{
			coef_oceanload.resize(6, 11);
			coef_oceanload = 0.0;
			if (_log)
			{
				_log->comment(1, site + " :Not found in ocean load file!!!!");
			}
		}
		//compute NEU for 11 components
		t_gtriple dx;
		for (int i = 0; i < 11; i++)
		{
			dx[0] = dx[0] - coef_oceanload(3,i + 1)*cos(angle[i] - coef_oceanload(6,i + 1)*D2R);
			dx[1] = dx[1] - coef_oceanload(2,i + 1)*cos(angle[i] - coef_oceanload(5,i + 1)*D2R);
			dx[2] = dx[2] + coef_oceanload(1,i + 1)*cos(angle[i] - coef_oceanload(4,i + 1)*D2R);
		}
		t_gtriple blh, xyz;
		xyz2ell(xRec, blh, false);
		neu2xyz(blh, dx, xyz);
		xyz = xyz / 1.e3;
		return xyz;
	}
	t_gtriple t_gtideIERS::load_atmosph()
	{
		return t_gtriple();
	}

	t_gtriple t_gtideIERS::tide_freq(const string & site,const t_gtriple & xRec,double gast)
	{
		t_gtriple blh;
		xyz2ell(xRec, blh, false);
		double lat = blh[0];
		double lon = blh[1];
		double height = blh[2];
		
		double dotl = -0.0000253*sin(lat)*cos(lat)*sin(gast + lon);
		t_gtriple xyz = xRec * dotl / xRec.norm();

		return xyz;
	}

	void t_gtideIERS::getMeanPole(double mjd, double& xpm, double& ypm)
	{
		// 51544 is the mjd2000
		// 55197 is the mjd2010
		double dt = (mjd - 51544) / 365.25;

		if (mjd > 55197)
		{
			xpm = 23.513 + 7.6141 * dt;
			ypm = 358.891 - 0.6287 * dt;
		}
		else
		{
			xpm = 55.974 + (1.8243 + (0.18413 + 0.007024 * dt) * dt) * dt;
			ypm = 346.346 + (1.7896 - (0.10729 - 0.000908 * dt) * dt) * dt;
		}

		//to seconds
		xpm = xpm * (1e-3);
		ypm = ypm * (1e-3);
	}
}

