/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <iostream>
#include <stdlib.h>

#include "gdata/gmonit.h" 
 
using namespace std;

namespace gnut {

/* ----------
 * constructor
 */
t_gmonit::t_gmonit( string id )
  : _moni_id(id) 
{}


/* ----------
 * destructor
 */
t_gmonit::~t_gmonit()
{}


/* ----------
 * basic class for monitoring purpose
 */
void t_gmonit::show( ostringstream& os, int verb )
{   
   os << _moni_id << " - method not implemented\n";
   return;
}

} // namespace
