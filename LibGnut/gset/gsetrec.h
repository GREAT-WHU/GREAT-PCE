/*
*
* @verbatim
	History

	@endverbatim
*
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

*
* @file        gsetrec.h
* @brief       implements receiver object setting class
* @author      Jan Dousa
* @version     1.0.0
* @date        2012-10-23
*
*/

#ifndef GSETREC_H
#define GSETREC_H

#define XMLKEY_REC "receiver"

#include <map>
#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gdata/grec.h"
#include "gset/gsetbase.h"
#include "gutils/gtypeconv.h"
#include "gutils/gtriple.h"
#include "gmodels/ggpt.h"

#define HSL_UNKNOWN -9999

using namespace std;

namespace gnut
{

	class LibGnut_LIBRARY_EXPORT t_gsetrec : public virtual t_gsetbase
	{
	public:
		/** @brief constructor */
		t_gsetrec();

		/** @brief destructor */
		~t_gsetrec();

		/** @brief settings check */
		void check();

		/** @brief settings help */
		void help();

		int get_crd_xyz(t_gtriple& xyz, string s);
		t_gtriple get_crd_xyz(string s);
		t_gtriple get_std_xyz(string s);
		set<string> objects();
		shared_ptr<t_grec> grec(string s, t_glog* glog = 0);

		t_gobj* obj(string s);
		string  rec(string s);
		string  ant(string s);

		/**
		 * @brief get the List of recevier names
		 * @return set<string> : List of recevier names
		 */
		virtual set<string> recs();
		set<string> all_rec();

	protected:
		/** @brief get crd xyz */
		t_gtriple _get_crd_xyz(string s);
		t_gtriple _get_std_xyz(string s);

		/** @brief get ecc neu */
		t_gtriple _get_ecc_neu(string s);
		t_gtriple _get_ecc_xyz(string s);

		/** @brief get crd blh */
		t_gtriple _get_crd_blh(string s);

		/** @brief Global Pressure Temperature model */
		t_gpt _ggpt;

		set<string> _objects();                        // get all objects names

		string _rec;                                   // default receiver name
		string _ant;                                   // default antenna name
		string _id;                                    // receiver id
		string _name;                                  // receiver name
		string _desc;                                  // receiver description
		string _domes;                                 // receiver monumentation domes
		t_gtime _beg;                                  // default begin time
		t_gtime _end;                                  // default end time
		double _X;                                     // receiver X-coordinate [m] 
		double _Y;                                     // receiver Y-coordinate [m] 
		double _Z;                                     // receiver Z-coordinate [m]
		double _dX;                                    // receiver X-eccentricity [m]
		double _dY;                                    // receiver Y-eccentricity [m] 
		double _dZ;                                    // receiver Z-eccentricity [m]
		double _dE;                                    // receiver E-eccentricity [m]
		double _dN;                                    // receiver N-eccentricity [m]
		double _dU;                                    // receiver U-eccentricity [m]   
		double _ZTD;
		double _ZTD_sig;
		string _CRD_sig;
		bool   _overwrite;

	private:
	};

} // namespace

#endif
