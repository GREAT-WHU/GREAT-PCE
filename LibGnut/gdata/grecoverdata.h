/**
 * @file         grecoverdata.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Storage all recover observ equation/ parameter information
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#ifndef GRECOVERDATA_H
#define GRECOVERDATA_H

#include <set>
#include "gutils/gtime.h"
#include "gmodels/gpar.h"
#include "gexport/ExportLibGnut.h"

using namespace gnut;

namespace great
{
	/** @brief class for grcover_head. */
	class LibGnut_LIBRARY_EXPORT t_grcover_head
	{
	public:
		/** @brief default constructor. */
		t_grcover_head();
		~t_grcover_head();
		/** @brief get begin time. */
		t_gtime get_beg_time() const;
		/** @brief get end time. */
		t_gtime get_end_time() const;

	public:

		double interval;        // interval
		double sigma0 = 0;      // sigma0
		set<string>  sat_list;  // list of sat
		set<string>  site_list; // list of site
		set<t_gtime> time_list; // list of time
	};
	/** @brief class for grcover_data. */
	class LibGnut_LIBRARY_EXPORT t_grecover_data : public t_gdata
	{

	public:
		/** @brief default constructor. */
		t_grecover_data();
		virtual ~t_grecover_data();

	public:
		virtual t_gtime get_recover_time() const = 0;
		virtual string convert2strline() const = 0;
	};
	/** @brief class for grcover_data. */
	class LibGnut_LIBRARY_EXPORT t_grecover_equation : public t_grecover_data
	{
	public:
		/** @brief default constructor. */
		t_grecover_equation(const t_gtime& time, const string site, string sat);
		t_grecover_equation(const t_grecover_equation& other);
		~t_grecover_equation();

		void operator =(const t_grecover_equation& other);

		t_gtime get_recover_time() const override;
		string convert2strline() const override;

		void set_recover_equation(const t_gobscombtype& obstype, const pair<double, double>& resinfo, int is_newamb);

	public:
		t_gtime time;
		string sat_name;
		string site_name;
		t_gobscombtype obstype;
		double weight;
		double resuidal;
		int is_newamb;
	};

	/** @brief class for grcover_par. */
	class LibGnut_LIBRARY_EXPORT t_grecover_par : public t_grecover_data
	{
	public:
		/** @brief default constructor. */
		t_grecover_par(const t_grecover_par& ohter);
		t_grecover_par(const t_gpar& par, double correct_value);
		~t_grecover_par();

		t_gtime get_recover_time() const override;
		string convert2strline() const override;

		bool operator<(const t_grecover_par&) const;

	public:
		t_gpar par;
		double correct_value;
	};

}

#endif // !GRECOVERDATA_H