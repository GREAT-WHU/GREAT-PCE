/**
 * @file         gtideIERS.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        precise model for computing correction
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GTIDEIERS_H
#define GTIDEIERS_H

#include "gexport/ExportLibGREAT.h"
#include "gmodels/gtide.h"
#include "gdata/gnavde.h"

namespace great
{

	/** @brief The class for gtide of IERS model.*/
	class LibGREAT_LIBRARY_EXPORT t_gtideIERS : public t_gtide
	{
	public:
		/** @brief constructor.*/
		t_gtideIERS(t_glog* l);
		virtual ~t_gtideIERS() {};

		/** @brief solid earth tides(Mark sure the pos in J2000).*/
		t_gtriple tide_solid(const t_gtime& epo, t_gtriple& xyz, Matrix& rot_trs2crs, t_gnavde* nav_planet);
		t_gtriple tide_solid(const t_gtime& epo, t_gtriple& xyz, Matrix& rot_trs2crs, ColumnVector& sun_pos, ColumnVector& moon_pos);

		/** @brief pole tides.*/
		t_gtriple tide_pole();
		/** @brief pole tides.*/
		t_gtriple tide_pole_pod(const t_gtime& epo, double xpole, double ypole, t_gtriple& xyz);
		/** @brief ocean tide loading.*/
		t_gtriple load_ocean(const t_gtime& epoch, const string& site, const t_gtriple& xRec) override;
		/** @brief atmospheric tide loading.*/
		t_gtriple load_atmosph() override;
		/** @brief get frequency of tide.*/
		t_gtriple tide_freq(const string& site, const t_gtriple & xRec, double gast);
		/** @brief get mean pole.*/
		void getMeanPole(double mjd, double& xpm, double& ypm);

	protected:


	}; 
}// namespace

#endif // !GTIDEPOD_H
