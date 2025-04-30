
/**
*
* @verbatim
	History
	2011-04-26  PV: created
	2012-04-06  JD: extracted matrix conversion utilities from old utils.h

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file        gmatrixconv.h
* @brief       Purpose: matrix conversion utilities
* @author      PV
* @version     1.0.0
* @date        2011-04-26
*
*/

#ifndef GMATRIXCONV_H
#define GMATRIXCONV_H

#include <ostream>   // need for NEWMATIO !
#include <vector>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "gexport/ExportLibGnut.h"

using namespace std;
namespace gnut
{

	// probably will be obsolete using gtriple everywhere !!!!!
	// -----
	LibGnut_LIBRARY_EXPORT ColumnVector array2colVec(double*);                     // convert array to columnVector

	LibGnut_LIBRARY_EXPORT void Matrix_remRC(SymmetricMatrix&, int row, int col);  // remove   r_th row and c_th column in SymMatrix
	LibGnut_LIBRARY_EXPORT void Matrix_remRC(DiagonalMatrix&, int row);  // remove   r_th row and c_th column in DiagonalMatrix
	LibGnut_LIBRARY_EXPORT void Matrix_rem(SymmetricMatrix&, vector<int>&);  // remove rows and columns stored in set
	LibGnut_LIBRARY_EXPORT void Matrix_addRC(SymmetricMatrix&, int row, int col);  // add zero r_th row and c_th column in SymMatrix
	LibGnut_LIBRARY_EXPORT void Matrix_addRC(Matrix&, int row, int col);  // add zero r_th row and c_th column in Matrix
	LibGnut_LIBRARY_EXPORT void Matrix_remR(Matrix&, int row);       // remove r_th row in Matrix (zhshen)
	LibGnut_LIBRARY_EXPORT void Matrix_remR(Matrix&, vector<int>&);  // remove rows row stored in set (zhshen)
	LibGnut_LIBRARY_EXPORT void Matrix_addRC(DiagonalMatrix&, int row);           // add zero r_th row and r_th column in DiagMatrix
	LibGnut_LIBRARY_EXPORT void Matrix_swap(SymmetricMatrix&, int a, int b);      // swap a_th and b_th row and column in SymMatrix
	LibGnut_LIBRARY_EXPORT void Matrix_cpRC(SymmetricMatrix, SymmetricMatrix&, int r, int c); // copy row and col
	LibGnut_LIBRARY_EXPORT void Vector_add(ColumnVector&, int row);           // add zero r_th row in Column Vector

	LibGnut_LIBRARY_EXPORT void indexing(const ColumnVector& v, const ColumnVector& v_sorted, vector<int>& index);

	LibGnut_LIBRARY_EXPORT DiagonalMatrix SR(DiagonalMatrix& D);
	LibGnut_LIBRARY_EXPORT Matrix rotX(double Angle);
	LibGnut_LIBRARY_EXPORT Matrix rotY(double Angle);
	LibGnut_LIBRARY_EXPORT Matrix rotZ(double Angle);

	LibGnut_LIBRARY_EXPORT SymmetricMatrix cov2corr(SymmetricMatrix& Q);
} // namespace

#endif