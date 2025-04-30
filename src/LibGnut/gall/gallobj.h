/**
*
* @verbatim
    History

*
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     gallobj.h
* @brief    container for all objects
*
* @author   JD
* @version  1.0.0
* @date     2012-12-04
*
*/

#ifndef GALLOBJ_H
#define GALLOBJ_H


#include <vector>
#include <memory>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gdata/grec.h"
#include "gdata/gtrn.h"
#include "gutils/gsysconv.h"
#include "gutils/gtypeconv.h"
#include "gall/gallotl.h"

using namespace std;

namespace gnut
{

	/**
	 *@brief Class for t_allobj derive from t_gdata
	 */
	class LibGnut_LIBRARY_EXPORT t_gallobj : public t_gdata
	{

	public:
		/** @brief default constructor. */
		t_gallobj();
		t_gallobj(t_gallpcv* pcv, t_gallotl* otl);
		virtual ~t_gallobj();

		typedef map<string, shared_ptr<t_gobj> > t_map_obj;          // allocated data of a single object

		void setPCV(t_gallpcv*  pcv);
		void setOTL(t_gallotl*  otl);
		/**
		 * @brief set/get single obj element
		 *
		 * @param obj
		 * @return int
		 */
		int add(set<string> objs);
		int add(shared_ptr<t_gobj> obj);
		shared_ptr<t_gobj> obj(string s);
		/**
		* @brief synchronize PCVs.
		*
		* @return    void
		*/
		void sync_pcvs();
		/**
		* @brief read satellite information.
		*
		* @param[in]  epo    epoch
		* @return    void
		*/
		void read_satinfo(t_gtime& epo);
		/**
		* @brief get all object elements
		*
		* @param[in]  id    id type
		* @return    all object elements
		*/
		virtual map<string, shared_ptr<t_gobj> > objects(t_gdata::ID_TYPE id = NONE);
		/**
		* @brief get object set
		*
		* @param[in]  id    id type
		* @return    object set
		*/
		virtual set<string> objects_set(t_gdata::ID_TYPE id = t_gdata::NONE);
		/**
		* @brief get all object elements.
		*
		* @param[in]  id
		* @return    all object
		*/
		virtual int obj_num();                                        // get object numbera

	private:
		t_map_obj      _mapobj;                                      // map of all objects
		t_gallpcv*     _gpcv;                                        // map of all PCV
		t_gallotl*     _gotl;                                        // map of all otl
		/** @brief alocate all t_gtrn objects for all GNSS. */
		void           _aloctrn();
	};

} // namespace

#endif
