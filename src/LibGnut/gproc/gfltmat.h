/**
*
* @verbatim
	History
	2011-01-10  PV: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gfltmat.h
* @brief      Purpose: implements gneq stacking class
*.
* @author     PV
* @version    1.0.0
* @date       2011-01-10
*
*/

#ifndef GFLTMAT_H
#define GFLTMAT_H

#include <map>
#include <vector>
#include <iostream>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"
#include "newmat/include.h"

#include "gmodels/gpar.h"
#include "gall/gallpar.h"
#include "gall/gallrslt.h"

namespace gnut
{

	// Po domluve muzeme spojit s t_gneq
	class LibGnut_LIBRARY_EXPORT t_gfltmat
	{
	public:
		//Constructor, destructor
		t_gfltmat() {};
		t_gfltmat(SymmetricMatrix Qp,
			SymmetricMatrix Qu,
			t_gallpar       xp,
			t_gallpar       xu);
		~t_gfltmat();

		SymmetricMatrix     Qp();
		void                Qp(const SymmetricMatrix& Qp);
		SymmetricMatrix     Qu();
		void                Qu(const SymmetricMatrix& Qu);
		DiagonalMatrix      Noise();
		void                Noise(const DiagonalMatrix& Noise);
		t_gallpar           xp();
		void                xp(const t_gallpar& xp);
		t_gallpar           xu();
		void                xu(const t_gallpar& xu);
		t_gtime             epo();
		void                epo(const t_gtime& t);
		set<string>         slips();
		void                slips(const set<string>& cs);
		vector<t_gsatdata>  data();
		void                data(const vector<t_gsatdata>& data);
		void                delParam(int i, int index);

	protected:
		SymmetricMatrix _Qp;
		SymmetricMatrix _Qu;
		DiagonalMatrix  _Noise;
		t_gallpar       _xp;
		t_gallpar       _xu;
		t_gtime         _epo;
		set<string>     _slips;
		vector<t_gsatdata> _data;
	};

} // namespace


#endif // GFLTMAT_H
