
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include "gmodels/gmodel.h"

namespace gnut {   

// Destructor
// -----------------
t_gmodel::~t_gmodel()
{
   
}

// Set site name
// -----------
void t_gmodel::setSite(string site)
{
   this->_site = site;
}

} // namespace
