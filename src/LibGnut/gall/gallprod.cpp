
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.


-*/

#include "gall/gallprod.h" 

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallprod::t_gallprod( )
  : t_gdata()
{
  gtrace("t_gallprod::constructor");
  id_type(  t_gdata::ALLPROD);
  id_group( t_gdata::GRP_PRODUCT );
}


// destructor
// ----------
t_gallprod::~t_gallprod()
{
  gtrace("t_gallprod::destructor");
  this->clear();
}


// clean
// --------------------
void t_gallprod::clear()
{
  gtrace("t_gallprod::clear");

  _gmutex.lock();
  _map_prod.clear();
  _gmutex.unlock();
}


// add product
// ---------------------
int t_gallprod::add(shared_ptr<t_gprod> prod, string site) // str only if prod->obj_pt = 0
{
  gtrace("gallprod::add");
   
  ID_TYPE type = prod->id_type();
   
  string id = prod->obj_id(); // if( id.empty() ){ return -1; } // no, use id instead!
   
  if( id.empty() ) id = site;
  if( id.empty() && ! site.empty() && id != site )
    cerr << "warning [gallprod] - adding product with both [and different] non-zero string and object id!\n";
   
  _gmutex.lock();
  _map_prod[id][type][prod->epoch()] = prod;  // allocate in heap

  _gmutex.unlock(); return 0;
}


// get product
// --------------------
shared_ptr<t_gprod> t_gallprod::get(const string& site, ID_TYPE type, const t_gtime& t)
{
  gtrace("gallprod::get(str)");
   
  _gmutex.lock(); 
   
  shared_ptr<t_gprod> obj_pt = _find(site, type, t);

  _gmutex.unlock(); return obj_pt;
}

// list of sites
// --------------------
set<string> t_gallprod::prod_sites()
{
  gtrace("t_gallprod::prod_sites");

  set<string> site_list;
  _gmutex.lock();
  t_map_prd::const_iterator it;  
  for( it = _map_prod.begin(); it != _map_prod.end(); ++it ){
     site_list.insert( it->first );
  }
  _gmutex.unlock(); return site_list;
}
      
// list of product types
// --------------------
set<t_gdata::ID_TYPE> t_gallprod::prod_types(const string& site)
{
  gtrace("t_gallprod::prod_types " + site );

  set<t_gdata::ID_TYPE> type_list;

  _gmutex.lock();

  t_map_id::const_iterator itTYPE;
  t_map_prd::const_iterator it = _map_prod.find( site );

  if( it == _map_prod.end() ) { _gmutex.unlock(); return type_list;}
   
  for( itTYPE = _map_prod[site].begin(); itTYPE != _map_prod[site].end(); ++itTYPE ){
     type_list.insert( itTYPE->first );
  }

  _gmutex.unlock(); return type_list;
}

// list of product types
// --------------------
set<t_gtime> t_gallprod::prod_epochs(const string& site, ID_TYPE type)
{
  gtrace("t_gallprod::prod_epochs " + site);

  set<t_gtime> epo_list;

  _gmutex.lock();

  t_map_epo::const_iterator itEPO;
  t_map_id::const_iterator itTYPE;
  t_map_prd::const_iterator it = _map_prod.find(site);

  if (it == _map_prod.end()) {
    _gmutex.unlock();
    return epo_list;
  } else {
    itTYPE = it->second.find(type);
    if (itTYPE == it->second.end()) {
      _gmutex.unlock();
      return epo_list;
    } else {
      for (itEPO = itTYPE->second.begin(); itEPO != itTYPE->second.end(); ++itEPO) {
        epo_list.insert(itEPO->first);
      }
    }
  }

  _gmutex.unlock();
  return epo_list;
}   

// clean data out of the interval
// ----------------------------   
void t_gallprod::clean_outer( const t_gtime& beg, const t_gtime& end )
{
  gtrace("t_gallprod::clean_outer");

  if( end < beg ) return;

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  

   // loop over all sites
   t_map_prd::iterator itKEY = _map_prod.begin();
   while( itKEY != _map_prod.end() ){
      string key = itKEY->first;

      t_map_id::iterator itID = itKEY->second.begin();
      while( itID != itKEY->second.end() ){
	 string id = t_gdata::type2str(itID->first);
	 
	 t_map_epo::iterator it;
	 t_map_epo::iterator itFirst = itID->second.begin();
	 t_map_epo::iterator itLast  = itID->second.end();	 
	 t_map_epo::iterator itBeg   = itID->second.lower_bound(beg);
	 t_map_epo::iterator itEnd   = itID->second.upper_bound(end);	 
	 
	 
	 // remove before BEGIN request
	 if( itBeg != itFirst ){
       
            it = itFirst;

	    // begin is last
	    if( itBeg == itLast ){
	       itBeg--;
	       if( _log ) _log->comment(2,"gallprod","products removed before: " + itFirst->first.str("%Y-%m-%d %H:%M:%S[%T] ") + key + " " + id);
	       
  	       while( it != itLast ) // itID->second.erase(it++);
	         if( (it->second).use_count() == 1 ){ it = itID->second.erase(it); }
	         else{ it++; }
	       
	    // begin is not last
	    }else{
	       if( _log ) _log->comment(2,"gallprod","products removed before: " + itBeg->first.str("%Y-%m-%d %H:%M:%S[%T] ") + key + " " + id);
	       
    	       while( it != itBeg ) // itID->second.erase(it++);
	         if( (it->second).use_count() == 1 ){ it = itID->second.erase(it); }
 	         else{ it++; }		       
	    }
	 }
	 
         // remove after END request
	 if( itEnd != itLast ){
	    if( _log ) _log->comment(2,"gallprod","products removed after : " + itEnd->first.str("%Y-%m-%d %H:%M:%S[%T] ") + key + " " + id);
	    it = itEnd;
	    
	    while( it != itLast ) // itID->second.erase(it++);
	      if( (it->second).use_count() == 1 ){ it = itID->second.erase(it); }
	      else{ it++; }
	 }
	 itID++;
      }
      itKEY++;
   }   
   
#ifdef BMUTEX   
  lock.unlock();
#endif
  _gmutex.unlock();
  return; 
}

// remove appropriate element t_gprod*
// ---------------------
void t_gallprod::rem(const string& site, ID_TYPE type, const t_gtime& t)
{
   gtrace("t_gallprod::rem(str)");
  
   _gmutex.lock();

   t_map_epo::iterator itEP;
   t_map_id::iterator itID;
   t_map_prd::iterator it = _map_prod.find( site );
   if( it != _map_prod.end() ){
      itID = it->second.find(type);
      if( itID != it->second.end() ){
	 itEP = itID->second.find(t);
	 if( itEP != itID->second.end() ) _map_prod[site][type].erase(itEP);
      }     
   }

   _gmutex.unlock();
   return;
}

// find appropriate element t_gprod*
// ---------------------
shared_ptr<t_gprod> t_gallprod::_find(const string& site, ID_TYPE type, const t_gtime& t)
{
  gtrace("t_gallprod::_find(str)");

  t_map_epo::iterator itEP;
  t_map_id::iterator itID;
  t_map_prd::iterator it = _map_prod.find(site);
  if (it != _map_prod.end()) {
    itID = it->second.find(type);
    if (itID != it->second.end()) {
      itEP = itID->second.find(t);
      if (itEP != itID->second.end())
        return itEP->second;
    }
  }

  shared_ptr<t_gprod> tmp;
  return tmp;
}

} // namespace
