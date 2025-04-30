/**
*
* @verbatim
History
*
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     gallprec.h
* @brief    Purpose: implementation of SP3 ephemeris container
*
*
* @author   JD
* @version  1.0.0
* @date     2011-01-10
*
*/

#ifndef GALLPREC_H
#define GALLPREC_H

#include "gall/gallnav.h"
#include "gdata/gephprec.h"
#include "gutils/gtriple.h"
#include "gmodels/gpoly.h"


using namespace std;

namespace gnut
{

	typedef map<string, double>       t_map_dat;     // single data for a single satellite
	typedef map<t_gtime, t_map_dat>       t_map_epo;     // all data for a single satellite
	typedef map<string, t_map_epo>       t_map_prn;     // all data for all satellites
	/**
	*@brief Class for t_gallprec derive from t_gallnav
	*/
	class LibGnut_LIBRARY_EXPORT t_gallprec : public t_gallnav
	{
	public:
		/** @brief clk type. */
		enum  clk_type {
			AS,
			AR,
			UNDEF
		};


	public:
		/** @brief default constructor. */
		t_gallprec();
		virtual ~t_gallprec();

		typedef map<string, shared_ptr<t_gephprec> > t_map_sp3;
		/**
		 * @brief inherited from gallnav to fix healthy
		 *
		 * @param sat
		 * @param t
		 * @return true
		 * @return false
		 */
		virtual bool health(string sat, const t_gtime& t);
		/**
		 * @brief pos
		 *
		 * @param sat
		 * @param t
		 * @param xyz
		 * @param var
		 * @param vel
		 * @param chk_mask
		 * @return int
		 */
		virtual int pos(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true); // [m]
		/**
		 * @brief  GNAV quality
		 *
		 * @param sat
		 * @param t
		 * @param xyz
		 * @param var
		 * @param vel
		 * @param chk_mask
		 * @return int
		 */
		virtual int nav(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true); // [m] GNAV quality
		/**
		 * @brief clk
		 *
		 * @param sat
		 * @param t
		 * @param clk
		 * @param var
		 * @param dclk
		 * @param chk_mask
		 * @return int
		 */
		virtual int clk(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL, bool chk_mask = true); // [s]
		virtual int clk_cdr(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL, double* ifcb = NULL, bool chk_mask = true);
		virtual int clk_sp3(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL); // [s]

		// ALTERNATIVES for direct interpolations - using gephprec!
		virtual int pos_int(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL); // [m]
		virtual int clk_int(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL); // [s]
		// ALTERNATIVES for direct interpolations !
		virtual int pos_alt(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL); // [m]
		virtual int clk_alt(string sat, const t_gtime& t, double* clk, double* var = NULL, double* dclk = NULL); // [s]
		int intv(string sat);
        int intv();
		
		bool posnew(string sat, const t_gtime& t, double  xyz[3]);
		void add_interval(string sat, int intv);
        void add_interval(int intv);
		void add_agency(string producer);
		void add_clk_interval(double intv);
		void add_orb_interval(double intv);
		int addvel(string sat, const t_gtime& ep, double xyzt[4], double var[4]);  // [m], [m]
		int addclk(string sat, const t_gtime& ep, double  clk[3], double var[3]);  // [s]
		int add_pos_vel(string sat, const t_gtime& ep, t_gtriple xyz, double t, t_gtriple dxyz, double dt,double xyzt[4], double var[4],int obs_num=0);  // [s]
		int addpos(string sat, const t_gtime& ep, t_gtriple xyz, double t, t_gtriple dxyz, double dt);  // [m], [usec]
		int add_delta_pos_vel(string sat, const t_gtime& ep, const int& iod, const t_gtriple& dxyz, const t_gtriple& dvxyz);
		int add_delta_clk(string sat, const t_gtime& ep, const int& iod, const double& dt, const double& dot_dt, const double& dot_dot_dt);
		
		int get_clk_correction(const string& sat, const t_gtime& t, int iod, double& clk, double& dclk);

		bool corr_avali();
		bool corr_avali(t_gtime now);
		/**
		 * @brief use clkrnx
		 *
		 * @param b
		 */
		void use_clkrnx(bool b) { _clkrnx = b; }
		/**
		 * @brief use clksp3
		 *
		 * @param b
		 */
		void use_clksp3(bool b) { _clksp3 = b; }
		/**
		 * @brief use clknav
		 *
		 * @param b
		 */
		void use_clknav(bool b) { _clknav = b; }
		/**
		 * @brief use posnav
		 *
		 * @param b
		 */
		void use_posnav(bool b) { _posnav = b; }
		/**
		 * @brief clean all sp3
		 */
		void clean_all() { _mapsp3.clear(); _mapprec.clear(); }
		/**
		 * @brief clean outer
		 *
		 * @param beg
		 * @param end
		 */
		void clean_outer(const t_gtime& beg = FIRST_TIME, const t_gtime& end = LAST_TIME);
		t_gtime get_beg();
		t_gtime get_end();
		string get_agency();
		vector<string> get_sats();
		vector<string> get_sat3();
		string get_data_type();
		string get_sat_type();
		int get_sp3_size();
		void set_beg(t_gtime beg);
		void set_end(t_gtime end);
		void set_sat(vector<string>sat);
		void set_sat3(vector<string> sat);
		void set_data_type(string type);
		void set_agency(string agency);
		void set_sat_type(string type);
		bool get_pos_vel(string sat, t_gtime epoch, double pos[3], double vel[3],int& obsnum);
		virtual set<string> satellites();                                      // get all satellites
		virtual unsigned int nepochs(const string& prn);                     // get number of epochs
		virtual t_gtime beg_data(string prn = "");                            // get first t_gephprec epoch
		virtual t_gtime end_data(string prn = "");                            // get  last t_gephprec epoch
		// IMPROVE beg_time/end_time to distinguish GALLNAV/GALLPREC - t_satview !
		virtual t_gtime beg_time(string prn = "") { return beg_gnav(prn); }    // get first t_gephprec epoch
		virtual t_gtime end_time(string prn = "") { return end_gnav(prn); }    // get  last t_gephprec epoch
		// =========================================================
		virtual t_gtime beg_clk(string prn = "");                            // get first precise clocks epoch
		virtual t_gtime end_clk(string prn = "");                            // get last precise clocks epoch
		virtual t_gtime beg_prec(string prn = "");                            // get first precise polynomials epoch
		virtual t_gtime end_prec(string prn = "");                            // get last precise polynomials epoch
        virtual set<string> clk_objs();										// get all clk satellites and receivers
		virtual set<t_gtime> clk_epochs();	


		set<clk_type> get_clk_type() const;




	protected:
		virtual shared_ptr<t_geph> _find(string sat, const t_gtime& t);  // find appropriate t_geph element
		virtual int         _get_crddata(string sat, const t_gtime& t);  // fill PT,X,Y,Z vectors
		virtual int         _get_clkdata(string sat, const t_gtime& t);  // fill CT,C vectors

		virtual int _get_delta_pos_vel(const string& sat, const t_gtime& t);
		virtual int _get_delta_pos_vel(const string& sat, const t_gtime& t, int iod, t_gtime& tRef, t_map_dat& orbcorr);

		virtual int _get_delta_clk(const string& sat, const t_gtime& t);
		virtual int _get_delta_clk(const string& sat, const t_gtime& t,int iod, t_gtime& tRef, t_map_dat& clkcorr);

		t_map_sat         _mapprec;     // map of sp3 polynomials
		t_map_prn         _mapsp3;      // precise orbits&clocks (SP3) - full discrete data sets
		t_map_prn         _mapclk;      // precise clocks (CLOCK-RINEX) - full discrete data sets

	private:
		t_map_sp3         _prec;        // CACHE: single SP3 precise ephemeris for all satellites
		unsigned int      _degree_sp3;  // polynom degree for satellite sp3 position and clocks
		double            _sec;         // default polynomial units
		t_gtime           _ref;         // selected reference epoch for crd data/polynomials
		t_gtime           _clkref;      // selected reference epoch for clk data/polynomials
		bool              _clkrnx;      // true: use               clk from Rinex Clocks
		bool              _clksp3;      // true: use alternatively clk from sp3 (~15min!)
		bool              _clknav;      // true: use alternatively nav (low-precise clocks)
		bool              _posnav;      // true: use alternatively nav (low-precise orbits)
		unsigned int	  _udclkInt;
		unsigned int	  _udorbInt;
		int               _intv;
        map<string,int>    _intvm;
		string            _agency;
		t_gtime           _tbeg;
		t_gtime           _tend;
		vector<string>    _sats;
		vector<string>    _sat3;
		string            _datatype;
		string            _sattype;
		string            _frame;
	  // CACHE for approximative solutions
		map<string, t_gtime>   _poly_beg;
		map<string, t_gtime>   _poly_end;
		map<string, t_gpoly>   _poly_x;
		map<string, t_gpoly>   _poly_y;
		map<string, t_gpoly>   _poly_z;
		
		// BEGIN OF TEMPORARY (ALTERNATIVE for direct interpolation)
		vector<double>    _PT;          // vector of time-difference (X for polynomials)
		vector<t_gtime>   _T;           // vector of full time       (X for polynomials)
		vector<double>    _X;           // vector of x-coordinate    (Y for polynomials)
		vector<double>    _Y;           // vector of y-coordinate    (Y for polynomials)
		vector<double>    _Z;           // vector of z-coordinate    (Y for polynomials)

		vector<double>    _CT;          // vector of time-difference (X for polynomials)
		vector<double>    _C;           // vector of clk correction  (Y for polynomials)



		vector<double>    _PTCorr;          // vector of time-difference for SSR Correction(X for polynomials)
		vector<t_gtime>   _TCorr;           // vector of full time for SSR Correction      (X for polynomials)
		vector<double>    _XCorr;           // vector of x-coordinate for SSR Correction   (Y for polynomials)
		vector<double>    _YCorr;           // vector of y-coordinate for SSR Correction   (Y for polynomials)
		vector<double>    _ZCorr;           // vector of z-coordinate for SSR Correction   (Y for polynomials)

		vector<double>    _CTCorr;          // vector of time-difference for SSR Correction(X for polynomials)
		vector<double>    _CCorr;           // vector of clk correction for SSR Correction (Y for polynomials)
	  // END OF TEMPORARY (ALTERNATIVE)
		// CLK TYPE
		set<clk_type> _clk_type_list;
	};

} // namespace

#endif
