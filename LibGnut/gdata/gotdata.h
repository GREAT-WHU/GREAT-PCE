/**
 * @file         gotdata.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        The class for storaging ocean_tide data.
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */

#ifndef  GOTDATA_H
#define  GOTDATA_H

#include <vector>

#include "gdata/gdata.h"
#include "gexport/ExportLibGnut.h"

using namespace std;
using namespace gnut;

namespace great
{

	/**
	*@brief	   Class for storaging one record data of ocean_tide file
	*
	* The class contains one line record data of ocean_tide include
	*		index of spheric harmoics and Doodson and coefficients
	*/
	class t_octrec
	{
	public:
		/** @brief default constructor. */
		t_octrec();

		/** @brief default destructor. */
		virtual ~t_octrec() {};

		int n;            ///< the index of spheric harmoics
		int m;            ///< the index of spheric harmoics
		int index[6];     ///< the index of Doodson
		double coeff[4];  ///< coefficients of the ocean tide(C_PROG, C_RETG, S_PROG, S_RETG)
	};

	/**
	*@brief	   Class for storaging ocean_tide data
	*
	* The class contains Doodson index and coeffcients of
	*		oceantide in different order and provide the
	*		interface to get data.
	*/
	class LibGnut_LIBRARY_EXPORT t_gocean_tide : public t_gdata
	{
	public:
		/** @brief default constructor. */
		t_gocean_tide();

		/** @brief default destructor. */
		virtual ~t_gocean_tide() {};

		/**
		* @brief add one oceanide record to _oct.
		* @param[in]   m		index m
		* @param[in]   n		index n
		* @param[in]   index	Doodson's index
		* @param[in]   coeff	coefficient(c_prog,c_retg,s_prog,s_cretg)
		*/
		void set_oct(const int& m, const int& n, int index[], double coeff[]);

		/**
		* @brief add a load coefficient value to vector _load.
		* @param[in]   load		load coefficient value
		*/
		void set_load(const double& load);


		/** @brief whether the ocean tide data is empty.
		* @return  bool
		*	@retval   true   ocean tide data is empty
		*   @retval   false  ocean tide data is existent
		*/
		bool is_empty();

		/**
		* @brief get all oceantide coefficient.
		* @return   all oceantide coefficient
		*/
		vector<t_octrec>& oct();

		/**
		* @brief get all load coefficient.
		* @return   all load coefficient
		*/
		vector<double>& load();


	protected:
		vector<t_octrec> _oct;		///< ocean_tide coefficient data
		vector<double> _load;		///< ocean_tide load coefficient data
	};
} //namespace

#endif