/**
 * @file         greldelay.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        reldelay model class
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#include "gmodels/greldelay.h"
#include "gutils/gconst.h"

double great::t_reldelay_model::reldelay(t_gtriple& crd_site,t_gtriple& vel_site, t_gtriple& crd_sat, t_gtriple& vel_sat)
{ 
	double reldelay = 0.0;
	//reldelay = 2.0*(DotProduct(crd_sat.crd_cvect(), vel_sat.crd_cvect()) - DotProduct(crd_site.crd_cvect(), vel_site.crd_cvect())) / CLIGHT;
	reldelay = 2.0 * (DotProduct(crd_sat.crd_cvect(), vel_sat.crd_cvect())) / CLIGHT;
	ColumnVector xsat = crd_sat.crd_cvect();
	ColumnVector xsite = crd_site.crd_cvect();

	double r = xsite.norm_Frobenius() + xsat.norm_Frobenius();
	ColumnVector xsat2site = (xsite - xsat);
	double r_site2sat = xsat2site.norm_Frobenius();

	reldelay += 2.0 * GM_CGCS / CLIGHT / CLIGHT * log((r + r_site2sat) / (r - r_site2sat));

	return reldelay;
}
