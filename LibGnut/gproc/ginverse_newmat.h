/**
 * @file         ginverse_newmat.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief	     header files of matrix inverse
 * @version      1.0
 * @date         2024-08-29
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */


#ifndef GINVERSE_NEWMAT_H
#define GINVERSE_NEWMAT_H

#include "gproc/ginverse.h"
#include "newmat/newmat.h"

namespace great
{
	template<class MyMatrix>
	/** @brief class for t_ginverse_newmat derive from t_ginverse. */
	class t_ginverse_newmat: public t_ginverse<MyMatrix>
	{
	public:
		t_ginverse_newmat(bool is_from_zero) :
			t_ginverse<MyMatrix>(is_from_zero)
		{
		}

		virtual void sovle_NEQ(int num, MyMatrix& NEQ, const MyMatrix& W, MyMatrix& X) override
		{
			Matrix neq(num, num);
			Matrix w(num, 1);
			for (int i = 0; i < num; i++) 
			{
				for (int j = 0; j < num; j++)
				{
					neq(i + 1, j + 1) = NEQ(i + this->_beg, j + this->_beg);
				}
				w(i + 1,1) = W(i + this->_beg,this->_beg);
			}

			// using Newmat Matrix Lib
			LowerTriangularMatrix L = Cholesky(neq);
			neq << L.i();
			neq << neq.t() * neq;
			w << neq * w;

			for (int i = 0; i < num; i++)
			{
				for (int j = 0; j < num; j++)
				{
					NEQ(i + this->_beg, j + this->_beg) = neq(i + 1, j + 1);
				}
				X(i + this->_beg,this->_beg) = w(i + 1,1);
			}
		}

		virtual void inverse_NEQ(int num, MyMatrix& NEQ) override
		{
			Matrix neq(num, num);
			for (int i = 0; i < num; i++)
			{
				for (int j = 0; j < num; j++)
				{
					neq(i + 1, j + 1) = NEQ(i + this->_beg, j + this->_beg);
				}
			}

			LowerTriangularMatrix L = Cholesky(neq);
			neq << L.i();
			neq << neq.t() * neq;

			for (int i = 0; i < num; i++)
			{
				for (int j = 0; j < num; j++)
				{
					NEQ(i + this->_beg, j + this->_beg) = neq(i + 1, j + 1);
				}
			}
		}

	};

}
#endif // !GINVERSE_NEWMAT_H
