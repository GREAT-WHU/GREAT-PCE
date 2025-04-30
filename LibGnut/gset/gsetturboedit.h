/**
 * @file         gsetturboedit.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#ifndef GSETTURBOEDIT_H
#define GSETTURBOEDIT_H

#define XMLKEY_TURBOEDIT             "turboedit"

#include <map>
#include <string>
#include <iostream>

#include "gset/gsetbase.h"
#include "gexport/ExportLibGnut.h"

using namespace std;
using namespace gnut;

namespace great
{


	class LibGnut_LIBRARY_EXPORT t_gsetturboedit : public virtual t_gsetbase
	{
	public:
		/** @brief default constructor. */
		t_gsetturboedit();

		/** @brief default destructor. */
		virtual ~t_gsetturboedit() {};

		bool liteMode();

		bool isAmbOutput();

		bool isEphemeris();

		bool checkPC(double& pc_limit);

		bool checkMW(double& mw_limit);

		bool checkGF(double& gf_limit, double& gf_rms_limit);

		bool checkSingleFreq(double& sf_limit);

		bool checkGap(int& gap_limit);

		int smoothWindows();

		bool checkShort(int& short_limit);

		bool checkStatistics(double& minPercent, int& minMeanNprn, int& maxMeanNamb);

		/** @brief settings help. */
		void help();

	protected:
		
		double _defaultPCLimit;
		double _defaultMWLimit;
		double _defaultGFLimit;
		double _defaultGFRmsLimit;
		double _defaultSingleFreqLimit;
		int    _defaultGapArcLimit;
		int    _defaultShortArcLimit;
		double _defaultMinPercent;
		int    _defaultMinMeanNprn;
		int    _defaultMaxMeanNamb;

	};
}
#endif // !SETTURBOEDIT_H
