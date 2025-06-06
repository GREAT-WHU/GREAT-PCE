/**
*
* @verbatim
    History
    2011-01-10  PV: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gfltmat.cpp
* @brief      Purpose: implements gneq stacking class
*.
* @author     PV
* @version    1.0.0
* @date       2011-01-10
*
*/

#include <stdio.h>
#include <math.h>

#include "gproc/gfltmat.h"
#include "gutils/gmatrixconv.h"

using namespace std;

namespace gnut {

// t_gfltmat class
//=====================

t_gfltmat::t_gfltmat( SymmetricMatrix Qp,
                      SymmetricMatrix Qu,
                      t_gallpar       xp,
                      t_gallpar       xu )
{
   this->_Qp = Qp;
   this->_Qu = Qu;
   this->_xp = xp;
   this->_xu = xu;
}

t_gfltmat::~t_gfltmat()
{
   _slips.clear();
   _data.clear();
}


SymmetricMatrix t_gfltmat::Qp()
{
   return _Qp;
}

void t_gfltmat::Qp(const SymmetricMatrix& Qp)
{
   this->_Qp = Qp;
}

SymmetricMatrix t_gfltmat::Qu()
{
   return _Qu;
}

void t_gfltmat::Qu(const SymmetricMatrix& Qu)
{
   this->_Qu = Qu;
}

DiagonalMatrix t_gfltmat::Noise()
{
   return _Noise;
}

void t_gfltmat::Noise(const DiagonalMatrix& Noise)
{
   this->_Noise = Noise;
}


t_gallpar t_gfltmat::xp()
{
   return _xp;
}

void t_gfltmat::xp(const t_gallpar& xp)
{
   this->_xp = xp;
}

t_gallpar t_gfltmat::xu()
{
   return _xu;
}

void t_gfltmat::xu(const t_gallpar& xu)
{
   this->_xu = xu;
}

t_gtime t_gfltmat::epo()
{
   return _epo;
}

void t_gfltmat::epo(const t_gtime& epo)
{
   this->_epo = epo;
}
   
set<string> t_gfltmat::slips()
{
   return _slips;
}

void t_gfltmat::slips(const set<string>& cs)
{
   this->_slips = cs;
}
   
vector<t_gsatdata> t_gfltmat::data()
{
   return _data;
}

void t_gfltmat::data(const vector<t_gsatdata>& data)
{
   this->_data = data;
}   

void t_gfltmat::delParam(int i, int index)
{
   if (i<0 || index <= 0) return;
   
   Matrix_remRC(_Qu, index, index);  
   Matrix_remRC(_Qp, index, index);  

   _xu.delParam(i);
   _xu.reIndex();
   _xp.delParam(i);
   _xu.reIndex();
   
}

} // namespace
