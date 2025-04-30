
/**
* @verbatim
	History
	2012-05-11  JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gpair.h
* @brief       Purpose: implements 2D coordinates representation (e.g. horizontal coordinates)
* @author      JD
* @version     1.0.0
* @date        2012-05-11
*
*/

#ifndef GPAIR_H
#define GPAIR_H


#include "gexport/ExportLibGnut.h"

#include <iostream>
#include <string.h>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"

using namespace std;

namespace gnut
{
	/** @brief class for t_gpair. */
	class LibGnut_LIBRARY_EXPORT t_gpair
	{

	public:
		/** @brief default constructor. */
		t_gpair();
		t_gpair(double x, double y);
		t_gpair(double crd[2]);
		explicit t_gpair(const ColumnVector& crd);
		virtual ~t_gpair();

		t_gpair& operator=(const t_gpair& other);       // assignment operator
		t_gpair    operator+(const t_gpair& other) const; //
		bool       operator==(const t_gpair& tr) const;   // equal operator
		bool       operator<(const t_gpair& tr) const;
		double& operator[](const size_t idx);          // get a reference of element
		double     operator[](const size_t idx) const;    // get value of element

		friend ostream& operator<<(ostream& os, const t_gpair& x);

		double        crd(int idx) const;                 // get single element
		void          set(int idx, double newValue);      // set single element
		void          set(const ColumnVector&);           // set array by ColumnVector
		void          set(double crd[2]);                 // set array by array
		double* crd_array();                        // get array
		ColumnVector  crd_cvect();                        // get ColumnVector
		t_gpair& crd_pair();                         // get pair
		ColumnVector  unitary();                          // get unit ColumnVector
		bool          zero();                             // true: zero elements, false: not zero elements

	protected:

	private:
		double         _crd[2];   ///< Two-dimensional coordinates

	};

} // namespace

#endif
