
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
 
-*/

#include <iostream>
#include <iomanip>

#include "gutils/gmutex.h"

namespace gnut {

// constructor
// ----------
t_gmutex::t_gmutex()
{
#ifdef USE_OPENMP
	omp_init_lock(&_mutex);
#endif // USE_OPENMP

}

t_gmutex::t_gmutex(const t_gmutex & Ohter)
{
#ifdef USE_OPENMP
	omp_init_lock(&this->_mutex);
#endif
}


// destructor
// ----------
t_gmutex::~t_gmutex()
{
#ifdef USE_OPENMP
	omp_destroy_lock(&this->_mutex);
#endif
}

t_gmutex t_gmutex::operator=(const t_gmutex & Other)
{
	return t_gmutex();
}


// ----------
void t_gmutex::lock()
{   
#ifdef USE_OPENMP
	omp_set_lock(&_mutex);
#else
	_mutex.lock();
#endif // USE_OPENMP
}
 

// ----------
void t_gmutex::unlock()
{
#ifdef USE_OPENMP
	omp_unset_lock(&_mutex);
#else
	_mutex.unlock();
#endif
}
   

} // namespace
