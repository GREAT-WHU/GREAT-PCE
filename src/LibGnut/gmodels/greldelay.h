/**
 * @file         greldelay.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        reldelay model class
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GRELDELAY_H
#define GRELDELAY_H 

#include "gexport/ExportLibGnut.h"
#include "gutils/gtriple.h"

using namespace std;
using namespace gnut;

namespace great
{
	/** @brief class for t_reldelay_model. */
	class LibGnut_LIBRARY_EXPORT t_reldelay_model
	{
	public:
		/** @brief default constructor. */
		t_reldelay_model() {};
		virtual ~t_reldelay_model() {};

		// mapping function
		double reldelay(t_gtriple& crd_site, t_gtriple& vel_site, t_gtriple& crd_sat, t_gtriple& vel_sat);

	protected:
	};
}

#endif