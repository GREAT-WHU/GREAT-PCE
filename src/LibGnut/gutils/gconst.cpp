
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <iomanip>

#include "gutils/gconst.h"

using namespace std;

namespace gnut {

// static map initializer
// ---------
t_map_refr REFR_COEF()
{
  t_map_refr m;
   
  m["ESSEN"]    = { K1_ESS, K2_ESS, K3_ESS };
  m["SMITH"]    = { K1_SMW, K2_SMW, K3_SMW };
  m["BEVIS"]    = { K1_BEV, K2_BEV, K3_BEV };
  m["THAYER"]   = { K1_THA, K2_THA, K3_THA };
  m["RUEGER"]   = { K1_RUE, K2_RUE, K3_RUE };
  m["FOELSCHE"] = { K1_FOE, K2_FOE, K3_FOE };

  return m;
}

} // namespace
