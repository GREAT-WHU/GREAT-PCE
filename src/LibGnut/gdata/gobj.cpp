/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <stdio.h>
#include <math.h>

#include "gdata/gobj.h"
#include "gutils/gtypeconv.h"
#include "gutils/gsysconv.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gobj::t_gobj()
 :  _id(""),
    _name(""),
    _overwrite(false)
{
  gtrace("t_gobj::constructor");
  id_type(OBJ);

  
   _pcvnull = 0;
   //_allnav = nullptr;
}

// destructor
// ----------
t_gobj::~t_gobj()
{
  gtrace("t_gobj::destructor");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

//_mapantpcv.clear();
  _mappcv.clear();
  _mapeccxyz.clear();
  _mapeccneu.clear();   
  _mapant.clear();
  _mapcrd.clear();   
//_params.clear();
  _gmutex.unlock();
}


// set id
// ----------
void t_gobj::id( string str )
{
   gtrace("t_gobj::it(str): "+str);   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _id = str;
  _gmutex.unlock();
}


// get name
// ----------
string t_gobj::id() const
{
   gtrace("t_gobj::id");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _id;
  _gmutex.unlock(); 
  return tmp;
}

// set overwrite
// ----------
void t_gobj::overwrite(bool overwrite)
{
   gtrace("t_gobj::overwrite(bool)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _overwrite = overwrite; 
  _gmutex.unlock();
}
   
// get overwrite
// ----------
bool t_gobj::overwrite()
{
   gtrace("t_gobj::overwrite");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  _gmutex.unlock();   return _overwrite;
}   

// set name
// ----------
void t_gobj::name( string str )
{
   gtrace("t_gobj::name(str): "+str);
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
 _gmutex.lock();

 _name = str;   
   
 _gmutex.unlock(); return; 
}  

// get name
// ----------
string t_gobj::name() const
{
   gtrace("t_gobj::name");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _name;
  _gmutex.unlock(); return tmp;
}


// set domes
// ----------
void t_gobj::domes( string str )
{
   gtrace("t_gobj::domes(str)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  _domes = str;   
   
  _gmutex.unlock(); return; 
}  

// get DOMES
// ----------
string t_gobj::domes() const
{
   gtrace("t_gobj::domes");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _domes;
  _gmutex.unlock(); return tmp;
}


// set description
// ----------
void t_gobj::desc( string str )
{
   gtrace("t_gobj::desc(str)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  _desc = str;   
   
  _gmutex.unlock(); return; 
}   

// get description
// ----------
string t_gobj::desc() const
{
   gtrace("t_gobj::desc");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = _desc;
  _gmutex.unlock(); return tmp;
}

// set ecc offsets w.r.t. center of mass/reference point
// ----------
void t_gobj::eccxyz(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end )
{
   gtrace("t_gobj::eccxyz(...)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  _eccxyz(ecc, beg, end);
   
  _gmutex.unlock(); return; 
}
   

// set ecc offsets w.r.t. center of mass/reference point
// ----------
void t_gobj::_eccxyz(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end )
{
  t_gtriple zero(0.0, 0.0, 0.0);
  t_mapecc_xyz::iterator it =  _mapeccxyz.find(beg);

  if( end < beg ){     
    string lg = "Warning: " + _id + " not valid end time (end<beg) for eccentricity xyz (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";
    if( _log ) _log->comment(0,"gobj", lg);
    else               cerr << "gobj: " << lg << endl;
    return;
  }
   
  // begin record
  if( it == _mapeccxyz.end() ){     // not exists
    _mapeccxyz[beg] = ecc;
     
  }else{                         // record exists
    if( it->first == LAST_TIME ||
        it->second == zero ){
      
      _mapeccxyz[beg] = ecc;
    }else{
      string lg = "Warning: " + _id + " valid object record cannot be changed for eccentricity xyz";
      if( _log && !_overwrite ) _log->comment(1,"gobj", lg);
      return;
    }
  }

  // control end of record (with new beg search)
  it = _mapeccxyz.find(beg); it++;
   
  // beg was last in map (add final empty record)
  if( it == _mapeccxyz.end() ){
    _mapeccxyz[end] = zero;

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
        if( it->second == zero ) _mapeccxyz.erase(it); // remove obsolete empty record
        _mapeccxyz[end] = zero;
      
      }else{                                        // too close to next record
        string lg = "Warning: object " + _id + " 'eccxyz' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg);
        else               cerr << "gobj: " << lg << endl;
      }
    }else if(end != it->first){
      string lg = "Warning: object " + _id + " 'eccxyz' end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
      if( _log ) _log->comment(2,"gobj",lg);
    }
  }

  // remove duplicated empty records
  t_mapecc_xyz::iterator itNEW = _mapeccxyz.begin();
  t_mapecc_xyz::iterator itOLD = itNEW;
  while( itOLD != _mapeccxyz.end() ){    
    if( ++itNEW != _mapeccxyz.end() ){
      if( ( itNEW->second == zero && itOLD->second == zero ) )
        //          ( itNEW->first == LAST_TIME ) )
      {
          itNEW = _mapeccxyz.erase(itNEW);
      }
    }
    itOLD = itNEW;
  }
  return;
}


// get ecc offsets (>=t) w.r.t. center of mass/reference point
// ----------
t_gtriple t_gobj::eccxyz(const t_gtime& t) const
{
   gtrace("t_gobj::eccyz");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple eccxyz(0.0, 0.0, 0.0);
  eccxyz = this->_eccxyz(t);   
     
  _gmutex.unlock(); return eccxyz;
}

// get ecc offsets (>=t) w.r.t. center of mass/reference point
// ----------   
t_gtriple t_gobj::_eccxyz(const t_gtime& t) const
{
	t_gtriple tmp(0.0, 0.0, 0.0);
	t_mapecc_xyz::const_iterator it = _mapeccxyz.upper_bound(t);

	if (it != _mapeccxyz.begin()) { it--; tmp = it->second; return tmp; }
  
	if ((it == _mapeccxyz.begin() || tmp.zero()) && _mapeccneu.size()>0) {  // not found
    if(_mapeccneu.upper_bound(t) == _mapeccneu.begin()) return tmp;
    
		t_gtriple eccneu(0.0, 0.0, 0.0);
		eccneu = this->_eccneu(t);

		if (!eccneu.zero()) {
			t_gtriple xyz = this->_crd(t);
			t_gtriple ell; xyz2ell(xyz, ell, false);
			neu2xyz(ell, eccneu, tmp);
		}
	}

	return tmp;
}

// return validity for eccxyx at epoch t
void t_gobj::eccxyz_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const
{
  gtrace("t_gobj::eccxyz_validity");   
  
  t_mapecc_xyz::const_iterator it = _mapeccxyz.upper_bound(t);
  if( it != _mapeccxyz.begin()  && it != _mapeccxyz.end()){
    end = it->first;
    it--;
    beg = it->first;
  }
}

// set ecc offsets w.r.t. center of mass/reference point
// ----------
void t_gobj::eccneu(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end )
{
   gtrace("t_gobj::eccneu(...)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  _eccneu(ecc, beg, end);
   
  _gmutex.unlock(); return; 
}
   
   
// set ecc offsets w.r.t. center of mass/reference point
// ----------
void t_gobj::_eccneu(const t_gtriple& ecc, const t_gtime& beg, const t_gtime& end )
{
  t_gtriple zero(0.0, 0.0, 0.0);
  t_mapecc_neu::iterator it =  _mapeccneu.find(beg);

  if( end < beg ){
     string lg = "Warning: " + _id + " not valid end time (end<beg) for eccentricity neu (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";
     if( _log ) _log->comment(0,"gobj", lg );
     else               cerr << "gobj: " << lg << endl;
     return;
  }
   
  // begin record
  if( it == _mapeccneu.end() ){     // not exists
      _mapeccneu[beg] = ecc;
     
  }else{                         // record exists
    if( it->first == LAST_TIME ||
        it->second == zero ){

       _mapeccneu[beg] = ecc;
    }else{
       string lg = "Warning: " + _id + " valid object record cannot be changed for eccentricity xyz";
       if( _log && !_overwrite ) _log->comment(1,"gobj", lg );
       return;
    }
  }

  // control end of record (with new beg search)
  it = _mapeccneu.find(beg); it++;
   
  // beg was last in map (add final empty record)
  if( it == _mapeccneu.end() ){
    _mapeccneu[end] = zero;

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
        if( it->second == zero ) _mapeccneu.erase(it); // remove obsolete empty record
        _mapeccneu[end] = zero;
      
      }else{                                        // too close to next record
        string lg = "Warning: object " + _id + " 'eccneu' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
        else               cerr << "gobj: " << lg << endl;
      }
    }else if (end != it->first){
      string lg = "Warning: object " + _id + " 'eccneu' end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
      if( _log ) _log->comment(2,"gobj",lg );
    }
  }

  // remove duplicated empty records
  t_mapecc_neu::iterator itNEW = _mapeccneu.begin();
  t_mapecc_neu::iterator itOLD = itNEW;
  while( itOLD != _mapeccneu.end() ){    
    if( ++itNEW != _mapeccneu.end() ){
      if( ( itNEW->second == zero && itOLD->second == zero ) )
        //          ( itNEW->first == LAST_TIME ) )
      {
          itNEW = _mapeccneu.erase(itNEW);
      }
    }
    itOLD = itNEW;
  }
  return;
}


// get ecc offsets (>=t) w.r.t. center of mass/reference point (interface only)
// ----------
t_gtriple t_gobj::eccneu(const t_gtime& t) const
{
   gtrace("t_gobj::eccneu");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple eccneu(0.0, 0.0, 0.0);
  eccneu = this->_eccneu(t);

  _gmutex.unlock(); return eccneu;
}
   
// get ecc offsets (>=t) w.r.t. center of mass/reference point
// ---------
t_gtriple t_gobj::_eccneu(const t_gtime& t) const
{
	t_gtriple tmp(0.0, 0.0, 0.0);
	t_mapecc_neu::const_iterator it = _mapeccneu.upper_bound(t);

	if (it != _mapeccneu.begin()) { it--; tmp = it->second; return tmp; }

	if ((it == _mapeccneu.begin() || tmp.zero()) && _mapeccxyz.size()>0) { // not found
    if(_mapeccxyz.upper_bound(t) == _mapeccxyz.begin()) return tmp;

		// transformation from XYZ if not available NEU
		t_gtriple eccxyz = this->_eccxyz(t);

		if (!eccxyz.zero()) {
			t_gtriple xyz = this->_crd(t);
			xyz2neu(eccxyz, xyz, tmp);
		}
	}

	return tmp;
}

// return validity for eccneu at epoch t
void t_gobj::eccneu_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const
{
  gtrace("t_gobj::eccneu_validity");   
  
  t_mapecc_neu::const_iterator it = _mapeccneu.upper_bound(t);
  if( it != _mapeccneu.begin() && it != _mapeccneu.end()){
    end = it->first;
    it--;
    beg = it->first;
  }
}

// set object position
// ----------
void t_gobj::crd(const t_gtriple& crd, const t_gtriple& std, const t_gtime& beg, const t_gtime& end, bool overwrite)
{
  gtrace("t_gobj::crd(t_gtriple, epo)");

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  _crd(crd, std, beg, end,overwrite);
   
  _gmutex.unlock(); return; 
}

// set object position
// ----------
void t_gobj::crd(const t_gtriple& crd, const t_gtriple& std)
{
	gtrace("t_gobj::crd(t_gtriple, epo)");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	_mapcrd.begin()->second.first = crd;
	_mapcrd.begin()->second.second = std;

	_gmutex.unlock(); return;
}
   
// set object position
// ----------
void t_gobj::_crd(const t_gtriple& crd, const t_gtriple& std, const t_gtime& beg, const t_gtime& end ,bool overwrite)
{
  pair<t_gtriple, t_gtriple> val;
  val = make_pair(crd, std);
  
  pair<t_gtriple, t_gtriple> zero;
  t_gtriple nullCRD(0.0, 0.0, 0.0);
  t_gtriple nullSTD(0.0, 0.0, 0.0);
  zero = make_pair(nullCRD, nullSTD);
  
  t_mapcrd::iterator it =  _mapcrd.find(beg);
   if( end < beg ){   
      string lg = "Warning: " + _id + " not valid end time (end<beg) for coordinates (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";
      if( _log ) _log->comment(0,"gobj",lg );
      else               cerr << "gobj: " << lg << endl;
      return;
  }
  // begin record
  if( it == _mapcrd.end() ){     // not exists
      _mapcrd[beg] = val;
  }else{                         // record exists
    if( it->first == LAST_TIME ||
        it->second == zero || overwrite){
       _mapcrd[beg] = val;
    }else{
       string lg = "Warning: " + _id + " valid object record cannot be changed (coordinates)";
       if( _log && !_overwrite ) _log->comment(1,"gobj",lg );
       return;
    }
  }
  // control end of record (with new beg search)
  it = _mapcrd.find(beg); it++;
   
  // beg was last in map (add final empty record)
  if( it == _mapcrd.end() ){
    _mapcrd[end] = zero;

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
		// change by ZHJ (Here not set the zero according to the before time）
 		auto beg_it = _mapcrd.find(beg);
		if (beg_it == _mapcrd.begin()) {
			if (it->second == zero) _mapcrd.erase(it); // remove obsolete empty record
			_mapcrd[end] = zero;
		}
		else {
			beg_it--;
			_mapcrd[end] = beg_it->second;
		}
        
      }else{                                        // too close to next record
        string lg = "Warning: object " + _id + " 'crd' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
        else               cerr << "gobj: " << lg << endl;
      }
    }else if(end != it->first){
		if (overwrite) {

			bool insert_flag = false;
			// if can't find end insert 
			if (_mapcrd.find(end) == _mapcrd.end()){
				_mapcrd.insert(make_pair(end, zero));
				insert_flag = true;
			}
			auto end_pre = _mapcrd.find(end); end_pre--;
			auto beg_post = _mapcrd.find(beg); beg_post++;
			if (insert_flag) {
				// assign pre crd info to end
				_mapcrd[end] = end_pre->second;
			}
			
			// erase [beg_post,end_pre+1)
			_mapcrd.erase(beg_post, ++end_pre);
		}
		else {
			string lg = "Warning: object " + _id + " 'crd' end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
			if (_log) _log->comment(2, "gobj", lg);
		}
    }
  }

  // remove duplicated empty records
  t_mapcrd::iterator itNEW = _mapcrd.begin();
  t_mapcrd::iterator itOLD = itNEW;
  while( itOLD != _mapcrd.end() ){    
    if( ++itNEW != _mapcrd.end() ){
      if( ( itNEW->second == zero && itOLD->second == zero ) )
        //          ( itNEW->first == LAST_TIME ) )
      {
          itNEW = _mapcrd.erase(itNEW);
      }
    }
    itOLD = itNEW;
  }

  return;
}


// get object position (using eccentricity)
// ----------
t_gtriple t_gobj::crd_arp(const t_gtime& t) const
{
   gtrace("t_gobj::crd_arp");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple marker(0.0, 0.0, 0.0);
  marker = this->_crd(t);

  // applying eccentricity
  t_gtriple arp(0.0, 0.0, 0.0);
  arp = marker + this->_eccxyz(t);
   
  _gmutex.unlock(); return arp;
}
   
// interface for _crd(const t_gtime& t) const
// ----------
t_gtriple t_gobj::crd(const t_gtime& t) const
{
   gtrace("t_gobj::crd");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple marker(0.0, 0.0, 0.0);
  marker = this->_crd(t);

  _gmutex.unlock(); return marker;
}

bool t_gobj::get_recent_crd(const t_gtime & t, const double& ref_std,t_gtriple& crd, t_gtriple& std)
{

	if (_mapcrd.size() == 0) {
		return false;
	}
	auto find_iter = _mapcrd.upper_bound(t);
	
	if (find_iter != _mapcrd.begin()) {
		find_iter--;
	}

	// find by time reverse traverse
	while (find_iter != _mapcrd.begin())
	{
		if (!find_iter->second.first.zero() && find_iter->second.second.norm() < ref_std) {
			crd = find_iter->second.first;
			std = find_iter->second.second;
			return true;
		}
		find_iter--;
	}

	// find_iter is the begin of mapcrd
	if (find_iter->second.second.norm() < ref_std) {
		crd = find_iter->second.first;
		std = find_iter->second.second;
		return true;
	}
	else {
		return false;
	}
}

bool t_gobj::get_adjacent_crd(const t_gtime& t, const double& ref_std, t_gtriple& crd, t_gtriple& std)
{

    if (_mapcrd.size() == 0) {
        return false;
    }
    auto find_iter = _mapcrd.lower_bound(t - 86400 * 7);

    if (find_iter != _mapcrd.begin()) {
        find_iter--;
    }

    // find by time reverse traverse
    while (find_iter != _mapcrd.begin())
    {
        if (!find_iter->second.first.zero() && find_iter->second.second.norm() < ref_std) {
            crd = find_iter->second.first;
            std = find_iter->second.second;
            return true;
        }
        find_iter--;
    }

    // find_iter is the begin of mapcrd
    if (find_iter->second.second.norm() < ref_std) {
        crd = find_iter->second.first;
        std = find_iter->second.second;
        return true;
    }
    else {
        return false;
    }
}

// interface for _std(const t_gtime& t) const
// ----------
t_gtriple t_gobj::std(const t_gtime& t) const
{
   gtrace("t_gobj::std");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtriple marker(0.0, 0.0, 0.0);
  marker = this->_std(t);

  _gmutex.unlock(); return marker;
}

// get object position (without eccentricity correction)
// ----------
t_gtriple t_gobj::_crd(const t_gtime& t) const
{  
#ifdef DEBUG
   for( t_mapcrd::const_iterator itE = _mapcrd.begin(); itE != _mapcrd.end(); ++itE  )
     cerr << "OBJ-CRD: searching EPOCH: " << t.str_ymdhms() 
          << "   found: " << itE->first.str_ymdhms() << "\n";
#endif

  t_gtriple crd(0.0, 0.0, 0.0);
  t_mapcrd::const_iterator it = _mapcrd.upper_bound(t);
  if( it == _mapcrd.begin() ){ return crd; }   // not found
  it--;
  t_gtriple tmp = it->second.first;

  return tmp;
}

// get object position std
// ----------
t_gtriple t_gobj::_std(const t_gtime& t) const
{  
#ifdef DEBUG
   for( t_mapcrd::const_iterator itE = _mapcrd.begin(); itE != _mapcrd.end(); ++itE  )
     cerr << "OBJ-CRD: searching EPOCH: " << t.str_ymdhms() 
          << "   found: " << itE->first.str_ymdhms() << "\n";
#endif

  t_gtriple std(0.0, 0.0, 0.0);
  t_mapcrd::const_iterator it = _mapcrd.upper_bound(t);
  if( it == _mapcrd.begin() ){ return std; }   // not found
  it--;
  t_gtriple tmp = it->second.second;

  return tmp;
}   

// return validity for crd at epoch t
void t_gobj::crd_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const
{
  gtrace("t_gobj::crd_validity");   
  
  t_mapcrd::const_iterator it = _mapcrd.upper_bound(t);
  if( it != _mapcrd.begin() ){
    end = it->first;
    it--;
    beg = it->first;
  }
}

// set pcv element
// ----------
void t_gobj::pcv(shared_ptr<t_gpcv> pcv, const t_gtime& beg, const t_gtime& end)
{
   gtrace("t_gobj::pcv(...)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  this->_pcv(pcv, beg, end);
  _gmutex.unlock();
}


// set pcv element (protected)
// ----------
void t_gobj::_pcv(shared_ptr<t_gpcv> pcv, const t_gtime& beg, const t_gtime& end)
{
   gtrace("t_gobj::_pcv");   
  
  /* t_mapantpcv::iterator itant = _mapantpcv.find(ant); */
  /* if (itant == _mapantpcv.end()){ */
  /*    _mapantpcv[ant][beg] = pcv; */
  /* }else{ */
  /*   t_mappcv mappcv = itant->second; */
     
    t_mappcv::iterator it = _mappcv.find(beg);

    if( end < beg ){
      string lg = "Warning: " + _id + " not valid end time (end<beg) for PCV (beg:" + beg.str_ymdhms() + " -> " + end.str_ymdhms() + ")";
      if( _log ) _log->comment(0,"gobj",lg );
      else               cerr << "gobj: " << lg << endl;
      return;
    }
   
    // begin record
    if( it == _mappcv.end() ){     // not exists
      _mappcv[beg] = pcv; 
    
    // record exists
    }else if( it->first == LAST_TIME || !it->second ){
      _mappcv[beg] = pcv; 

    }else{
      string lg = "Warning: " + _id + " valid object record cannot be changed (PCV)";
      if( _log && !_overwrite ) _log->comment(1,"gobj",lg );
      return;
    }

    // control end of record (with new beg search)
    it = _mappcv.find(beg); it++;
     
    // beg was last in map (add final empty record)
    if( it == _mappcv.end() ){
      _mappcv[end] = _pcvnull;
    }else{                                            // process end according to next value
      if( end < it->first ){                          // only if end is smaller then existing
        if( fabs(it->first - end) > 3600 ){           // significantly smaller!
          if( it->second == 0 ) _mappcv.erase(it); // remove obsolete empty record
          _mappcv[end] = _pcvnull;
          
        }else{                                        // too close to next record
          string lg = "Warning: object " + _id + " 'pcv' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
          if( _log ) _log->comment(2,"gobj",lg );
        }
      }else if(end != it->first){
        string lg = "Warning: object " + _id + " 'pcv' end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
      }
    }

    // remove duplicated empty records
    t_mappcv::iterator itNEW = _mappcv.begin();
    t_mappcv::iterator itOLD = itNEW;
    while( itOLD != _mappcv.end() ){    
      if( ++itNEW != _mappcv.end() ){
        if( ( !itNEW->second && !itOLD->second ) )
          //            ( itNEW->first == LAST_TIME ) )
        {
            itNEW = _mappcv.erase(itNEW);
        }
      }
      itOLD = itNEW;
    }     
//  }
}


// get pcv element (>=t)
// ----------
shared_ptr<t_gpcv> t_gobj::pcv(const t_gtime& t) const
{  
   gtrace("t_gobj::pcv(t)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  shared_ptr<t_gpcv> tmp = _pcv(t);
  _gmutex.unlock(); return tmp;
}

// get pcv element (>=t) (protected)
// ----------
shared_ptr<t_gpcv> t_gobj::_pcv(const t_gtime& t) const
{  
   gtrace("t_gobj::_pcv(t)");   

  shared_ptr<t_gpcv> pcv;   
  string ant;      

  t_mapant::const_iterator it = _mapant.upper_bound(t);
  if( it == _mapant.begin() ) {
    ant = "";  // antenna not found
    string lg = "Warning: unknown PCO (no antenna found in the object " + _id + " ) " + t.str_ymdhms();
    if( _log ) _log->comment(2,"gobj",lg );
    return _pcvnull;
  }else{
    ant = (--it)->second;     
  }

  t_mappcv::const_iterator it2 = _mappcv.upper_bound(t);
  if( it2 == _mappcv.begin() ) {     
     ostringstream ostr;
     ostr << "Warning: unknown PCO ( antenna " << left << setw(20) << ant << " not found in ATX ) " << t.str_ymdhms();
     if( _log ) _log->comment(2,"gobj",ostr.str() );
     return _pcvnull;
  }else{ 
    pcv = (--it2)->second;
  }

  if(pcv->pcvkey().compare(ant) != 0){
     if(pcv->pcvkey().compare(0, 16, ant, 0, 16) == 0){
       string lg = "Warning: PCO Used without considering randome " + pcv->pcvkey() + " " + ant;
       if (_log) _log->comment(2, "gobj", lg);
     }
     else{
       ostringstream ostr;
       ostr << "Warning: unknown PCO ( changed antenna " << left << setw(20) << ant << " not found in ATX ) " << t.str_ymdhms();
       if (_log) _log->comment(2, "gobj", ostr.str());
#ifdef DEBUG     
       cout << "gobj RESET: ant [" << ant << "] [" << pcv->pcvkey() << " " << t.str_ymdhms(" epo:") << "]\n";
#endif     
       return _pcvnull;
    }
  }
   
  return pcv;
}

// set antenna name
// ----------
void t_gobj::ant(string ant, const t_gtime& beg, const t_gtime& end)
{
   gtrace("t_gobj::and(...)");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  _ant(ant, beg, end); 
   
  _gmutex.unlock(); return; 
}
   

// set antenna name
// ----------
void t_gobj::_ant(string ant, const t_gtime& beg, const t_gtime& end)
{  
  t_mapant::iterator it =  _mapant.find(beg);

  if( end < beg ){
    string lg = "Warning: " + _id + " not valid end time (end<beg) for antenna (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";
    if( _log ) _log->comment(0,"gobj",lg );
    else               cerr << "gobj: " << lg << endl;
    return;
  }

  // begin record
  if( it == _mapant.end() ){     // not exists
      _mapant[beg] = ant;
     
  }else{                         // last value
    if( it->first == LAST_TIME ||
        it->second.empty() ){

       _mapant[beg] = ant;
    }else{                      // record exists
       string lg = "Warning: " + _id + " valid object record cannot be changed (antenna)";
       if( _log && !_overwrite ) _log->comment(1,"gobj",lg );
       return;
    }
  }

  // control end of record (with new beg search)
  it = _mapant.find(beg); it++;

  // beg was last in map (add final empty record)
  if( it == _mapant.end() ){
    _mapant[end] = "";

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
        if( it->second.empty() ) { 
          _mapant.erase(it); // remove obsolete empty record
          _mapant[end] = "";
        }else {
          _mapant[end] = it->second;
        }
        
      }else{                                        // too close to next record
        string lg = "Warning: object " + _id + " 'obj' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
        if( _log ) _log->comment(2,"gobj",lg );
        else               cerr << "gobj: " << lg << endl;
      }
    }else if(end != it->first){
      string lg = "Warning: object " + _id + " 'ant' " + ant + " end cut and tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");
      if( _log ) _log->comment(2,"gobj",lg );
    }
  }

  // remove duplicated empty records
  t_mapant::iterator itNEW = _mapant.begin();
  t_mapant::iterator itOLD = itNEW;
  while( itOLD != _mapant.end() ){    
    if( ++itNEW != _mapant.end() ){
      if( ( itNEW->second.empty() && itOLD->second.empty() ) )
        //          ( itNEW->first == LAST_TIME ) )
      {
          itNEW = _mapant.erase(itNEW);
      }
    }
    itOLD = itNEW;
  }

  return;
}


// get antenna name (>=t)
// ----------
string t_gobj::ant(const t_gtime& t) const
{  
   gtrace("t_gobj::ant()");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  string tmp = this->_ant(t);
  _gmutex.unlock(); return tmp;
}

// get antenna name (>=t) (protected)
// ----------
string t_gobj::_ant(const t_gtime& t) const
{  
   gtrace("t_gobj::_ant(t)");   
   
  t_mapant::const_iterator it = _mapant.upper_bound(t);
  if( it == _mapant.begin() ) return "";  // not found
  it--;
  return it->second;
}

     
// get time tags
// ----------
/*vector<t_gtime> t_gobj::ecc_id() const
{
  vector<t_gtime> tmp;
  t_mapecc_xyz::const_iterator itMAP = _mapeccxyz.begin();
  while( itMAP != _mapecc.end() ){
    tmp.push_back( itMAP->first );
    itMAP++;
  }
  return tmp;
}
*/

// return validity for antenna at epoch t
void t_gobj::ant_validity(const t_gtime& t, t_gtime& beg, t_gtime& end) const
{
  gtrace("t_gobj::ant_validity");   
  
  t_mapant::const_iterator it = _mapant.upper_bound(t);
  if( it != _mapant.begin() && it != _mapant.end()){
    end = it->first;
    it--;
    beg = it->first;
  }
}

// get time tags
// ----------
vector<t_gtime> t_gobj::pcv_id() const
{      
   gtrace("t_gobj::pcv_id");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<t_gtime> tmp;
//  t_mapantpcv::const_iterator itant = _mapantpcv.find(ant);
   
//  if (itant == _mapantpcv.end()){
//     return tmp;
//  }else{
  t_mappcv::const_iterator itMAP = _mappcv.begin();
  while( itMAP != _mappcv.end() ){
    tmp.push_back( itMAP->first );
    itMAP++;
  }
  _gmutex.unlock(); return tmp;     
//  }
   
}

     
// get time tags
// ----------
vector<t_gtime> t_gobj::ant_id() const
{
   gtrace("t_gobj::ant_id");   
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  vector<t_gtime> tmp = this->_ant_id();
  _gmutex.unlock(); return tmp;   
}

// get time tags (protected)
// ----------
vector<t_gtime> t_gobj::_ant_id() const
{   
   gtrace("t_gobj::_ant_id");   
   
  vector<t_gtime> tmp;
  t_mapant::const_iterator itMAP = _mapant.begin();
  while( itMAP != _mapant.end() ){
    tmp.push_back( itMAP->first );
    itMAP++;
  }
  return tmp;
}

// get time tags
// ----------
vector<t_gtime> t_gobj::crd_id() const
{
	gtrace("t_gobj::crd_id");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	vector<t_gtime> tmp;
	t_mapcrd::const_iterator itMAP = _mapcrd.begin();
	while (itMAP != _mapcrd.end()) {
		tmp.push_back(itMAP->first);
		itMAP++;
	}
	_gmutex.unlock(); return tmp;
}

// set mappcv for all t_gobj
// -----------------------
void t_gobj::sync_pcv(t_gallpcv* pcvs)
{
  gtrace("t_gobj::sync_pcv");   
   
  if( pcvs == 0 ) return;

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<t_gtime> vant = this->_ant_id();

  if( vant.size() == 0 ){  // ant in obj is not set
    if( _log ) _log->comment(1,"gobj","Warning: " + _id + " ant name not set: " );
    else                                    cerr << "gobj - Warning: " + _id + " ant name not set:" << endl;

  }else{
    for( vector<t_gtime>::iterator it = vant.begin(); it != vant.end(); it++ ){

      // set pcv for all antennas
      string ant = this->_ant(*it);
      t_gtime epo = *it;  
   
      shared_ptr<t_gpcv> gpcv = pcvs->gpcv(ant, "*", epo); 
      if( gpcv != _pcvnull ) {
        this->_pcv(gpcv, epo);
      }else{
        ostringstream ostr;
        ostr << "Warning: unknown PCO ( antenna " << left << setw(20) << ant << " not found in ATX ) " << epo.str_ymdhms();
        if( _log && !ant.empty() ) _log->comment(1,"gobj",ostr.str() );
        else if( _log &&  ant.empty() ) _log->comment(2,"gobj",ostr.str() );    
      }     
   
#ifdef DEBUG
       cout << " gobj.cpp - sync for ant " << ant << "  pcv found [" << gpcv << "] " << it->str_ymdhms() << endl;
#endif
    }
  }

  _gmutex.unlock(); return;
}



// check consistency
// ----------
void t_gobj::compare(shared_ptr<t_gobj> gobj, const t_gtime& tt, string source)
{
  gtrace("t_gobj::compare");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif

  _gmutex.lock();
  
  string old, alt;
  t_gtriple trip_old, trip_alt;
  t_gtriple std_old, std_alt;  

  // Name
  old = trim(_name); alt = trim(gobj->name());

  if( old != alt && !alt.empty()){
    if( old.empty() )
    {
      _name = alt;
      //if( _log ) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (Name): " + alt + " (" + tt.str_ymdhms() + ")");    
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (Name): " + alt + " (" + tt.str_ymdhms() + ")");
    }
    else if( _overwrite ) 
    {
      _name = alt;
      // if( _log ) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source + " (Name): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source + " (Name): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
    }
    //else if( _log ) _log->comment(1,"gobj","Warning: object " + _id + " does not match " + source + " (Name): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
    else if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " does not match " + source + " (Name): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
  }
   
  // Domes
  old = trim(_domes); alt = trim(gobj->domes());
  
  if( old != alt && !alt.empty())
  {      
    if( old.empty() )
    {
      _domes = alt;
      //if( _log ) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (Domes): " + alt + " (" + tt.str_ymdhms() + ")");  
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (Domes): " + alt + " (" + tt.str_ymdhms() + ")");
    }
    else if( _overwrite ) 
    {
      _domes = alt;
      //if( _log ) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source + " (Domes): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source + " (Domes): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
    }
    //else if( _log ) _log->comment(1,"gobj","Warning: object " + _id + " does not match " + source + " (Domes): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
    else if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " does not match " + source + " (Domes): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
  }   
   
  // Antenna
  old = trim(_ant(tt)); alt = trim(gobj->ant(tt));
  t_gtime beg, end;
  gobj->ant_validity(tt, beg, end);
  
  if( old != alt && !alt.empty())
  {
    if( old.empty() )
    {
      _ant(alt, beg, end);
      this->ant_validity(tt, beg, end);
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (Antenna): " + alt
                               + " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");    
    }
    else if( _overwrite ) 
    {
      _ant(alt, beg, end);
       if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source + " (Antenna): " + old + " -> " + alt +
                                " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");
    }else if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " setting does not match " + source + " (Antenna): " + old + " -> " + alt + " (" + tt.str_ymdhms() + ")");
  }

  // Coordinates
  trip_old = _crd(tt); trip_alt = gobj->crd(tt);
  std_old  = _std(tt); std_alt  = gobj->std(tt);
  gobj->crd_validity(tt, beg, end);  
  
  t_gtriple diff = trip_old - trip_alt;
  double dist = diff.norm();
  if( (dist > 10 || std_alt < std_old) && !trip_alt.zero() ) {   // fixed max diff 10 m
    if( trip_old.zero() || std_alt < std_old){                   // add std_alt < std_old, 防止出现坐标为snx坐标，std仍旧是10/999的情况
      _crd(trip_alt, std_alt, beg, end,true);
      this->crd_validity(tt, beg, end);      
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (Coordinates): " + dbl2str(trip_alt[0]) + " "
                                                                                                + dbl2str(trip_alt[1]) + " "
                                                                                                + dbl2str(trip_alt[2]) +
                               " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");
    }else if( _overwrite ){
      _crd(trip_alt, std_alt, beg, end,true);
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source + " (Coordinates): " + dbl2str(trip_old[0]) + " "
                                                                                               + dbl2str(trip_old[1]) + " "
                                                                                               + dbl2str(trip_old[2]) + " -> "
                                                                                               + dbl2str(trip_alt[0]) + " "
                                                                                               + dbl2str(trip_alt[1]) + " "
                                                                                               + dbl2str(trip_alt[2]) +
                               " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");              

    }else if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " does not match " + source + " (Coordinates)" + " (" + tt.str_ymdhms() + ")");
  }

  // NEU eccentricity   
  trip_old = _eccneu(tt); trip_alt = gobj->eccneu(tt); 
  gobj->eccneu_validity(tt, beg, end);
  
  if( trip_old != trip_alt && !trip_alt.zero()){
    if( trip_old.zero() ){
      _eccneu(trip_alt, beg, end);
      this->eccneu_validity(tt, beg, end);
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (NEU Eccentricity): " + dbl2str(trip_alt[0]) + " "
                                                                                                     + dbl2str(trip_alt[1]) + " "
                                                                                                     + dbl2str(trip_alt[2]) +
                               " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");
    }else if( _overwrite )
    {
      _eccneu(trip_alt, beg, end);
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source + " (NEU Eccentricity): " + dbl2str(trip_old[0]) + " "
                                                                                                    + dbl2str(trip_old[1]) + " "
                                                                                                    + dbl2str(trip_old[2]) + " -> "
                                                                                                    + dbl2str(trip_alt[0]) + " "
                                                                                                    + dbl2str(trip_alt[1]) + " "
                                                                                                    + dbl2str(trip_alt[2]) +
                               " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");            

    }else if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " does not match " + source + " (NEU Eccentricity)" + " (" + tt.str_ymdhms() + ")");
  }

  // xyz eccentricity
  trip_old = _eccxyz(tt); trip_alt = gobj->eccxyz(tt);
  gobj->eccxyz_validity(tt, beg, end);
    
  if( trip_old != trip_alt && !trip_alt.zero()){
    if ( trip_old.zero() ) {
      _eccxyz(trip_alt, beg, end);
      this->eccxyz_validity(tt, beg, end);
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " completed by " + source + " (XYZ Eccentricity): " + dbl2str(trip_alt[0]) + " "
                                                                                                     + dbl2str(trip_alt[1]) + " "
                                                                                                     + dbl2str(trip_alt[2]) +
                               " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");
    }else if( _overwrite ){     
      _eccxyz(trip_alt, beg, end);
      if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " modified  by " + source+ " (XYZ Eccentricity): " + dbl2str(trip_old[0]) + " "
                                                                                                    + dbl2str(trip_old[1]) + " "
                                                                                                    + dbl2str(trip_old[2]) + " -> "
                                                                                                    + dbl2str(trip_alt[0]) + " "
                                                                                                    + dbl2str(trip_alt[1]) + " "
                                                                                                    + dbl2str(trip_alt[2]) +
                               " (" + beg.str_ymdhms() + "->" + end.str_ymdhms() + ")");            

    }else if( _log && _log->verb() >= 1) _log->comment(1,"gobj","Warning: object " + _id + " does not match " + source + " (XYZ Eccentricity)" + " (" + tt.str_ymdhms() + ")");
  }   
  
  _gmutex.unlock(); return;
}


// operator for sorting
// ----------
bool t_gobj::operator<(const t_gobj& t) const 
{ 
 return _name < t.name();
}
   

// operator for sorting
// ----------
bool t_gobj::operator==(const t_gobj& t) const
{
 return  _name == t.name();
}
   
} // namespace
