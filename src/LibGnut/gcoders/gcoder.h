/**
*
* @verbatim
    History
    -1.0    JD        2011-04-20    creat the file.
    -1.1    JD        2013-03-07    general settings supported
@endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file     gcoder.h
* @brief    The base class for decoding
*
* @author   PV
* @version  1.1.0
* @date     2013-03-29
*
*/

#ifndef GCODER_H
#define GCODER_H 
 

#include <set>
#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <memory>
#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gio/gio.h"
#include "gio/gnote.h"
#include "gdata/gdata.h"
#include "gutils/gsys.h"
#include "gutils/gmutex.h"
#include "gutils/gcommon.h"
#include "gset/gsetbase.h"
#include "gset/gsetgnss.h"
#include "gset/gsetgen.h"
#include "gcoders/gcoder_buffer.h"

#define BUFFER_INCREASE_FAC  1.5
#define DEFAULT_BUFFER_SIZE  4096
#define MAXIMUM_BUFFER_SIZE  240000 // 9000000
using	 namespace std;
using    namespace great;

namespace gnut {

    using namespace std;

    class t_gio;
    /**
    *@brief       Class for decoding
    */
    class LibGnut_LIBRARY_EXPORT t_gcoder {

     public:
         /**
         * @brief default constructor.
         */
        t_gcoder();
        /**
        * @brief override constructor 1.
        *
        * @param[in]  s            setbase control
        * @param[in]  version    version of the gcoder
        * @param[in]  sz        size of the buffer
        * @param[in]  id        string for reporting
        */
        t_gcoder(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE, string id = "gcoder" );
        /**
        * @brief override constructor 2.
        * @param[in]  beg        begin time
        * @param[in]  end        end time
        * @param[in]  s            setbase control
        * @param[in]  version    version of the gcoder
        * @param[in]  sz        size of the buffer
        * @param[in]  id        string for reporting
        */
        t_gcoder(t_gtime beg, t_gtime end, t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE, string id = "gcoder");
        t_gtime beg_epoch = t_gtime(0, 0); ///<        begin epoch
        t_gtime end_epoch = t_gtime(0, 0); ///<        end epoch
        t_gtime epoch = t_gtime(0, 0);     ///<        epoch
        virtual ~t_gcoder();

        typedef map<int, pair<GOBS,int> > t_vobstypes;
        /**
        * @brief clear the buffer
        * @return void
        */
        virtual void clear();

        virtual void   version(string s){ _version = s; }
        virtual string version(){ return _version; }

        virtual void    out_length(int len){ _out_len = len; }
        virtual int     out_length(){ return _out_len; }

        virtual void    out_sample(float smp){ _out_smp = smp; }
        virtual float   out_sample(){ return _out_smp; }

        virtual void    out_epoch(t_gtime epo){ _out_epo = epo; }
        virtual t_gtime out_epoch(){     return _out_epo; }
        /** @brief decode head. */
        virtual int decode_head( char* buff, int sz,           vector<string>& errmsg ){ return 0; } // = 0;
        /** @brief decode data. */
        virtual int decode_data( char* buff, int sz, int& cnt, vector<string>& errmsg ){ return 0; } // = 0;
        /** @brief encode head. */
        virtual int encode_head( char* buff, int sz,           vector<string>& errmsg ){ return 0; } // = 0;
        /** @brief encode data. */
        virtual int encode_data( char* buff, int sz, int& cnt, vector<string>& errmsg ){ return 0; } // = 0;
   
        void glog(t_glog* l){ _log = l; }             // set/get glog pointer
        t_glog* glog(){ return _log; }
        /** @brief set/get path. */
        void path(string s);
        string path(){ return _fname; }
        /** @brief set gio_ptr. */
        void add_gio(weak_ptr<t_gio> p){ _gio_ptr = p; }
        /** @brief add/get data. */
        int add_data( string data_id, t_gdata* data );
        /** @brief return data. */
        t_gdata* data( string data_id ){ return _data[data_id]; }
        /**
         * @brief brief set/get notes (messages/warning/errors)
         *
         * @return const vector<t_gnote>&
         */
        void            mesg(t_note m, string s);     // set/get notes (messages/warning/errors)
        vector<t_gnote> mesg();
        /** @brief set close with warning. */
        void close_with_warning(bool b){ _close_with_warning = b; }

        void pgm(string pgm) {_pgm = pgm;}
   
     protected:
        int   endpos(){ return _endpos; }
        int     size(){ return _endpos; }             // number of char elements in buffer (excl \0)
        char* buffer(){ return _buffer; }

        virtual int  _gset(t_gsetbase* s);            // get settings (NO MORE PUBLIC, only via CONSTRUCT)
        /** @brief init. */
        virtual void _init();
        /**
        * @brief add data.
        * @param[in]  id        data type
        * @param[in]  data        t_gdata
        * @return void
        */
        virtual void _add_data(string id,t_gdata* data){};
        /**
        * @brief sampling filter for epochs.
        * @param[in]  epo        current epoch
        * @return
        *        true:the epoch fits sampling
        *        false:the epoch does not fit sampling
        */
        virtual bool _filter_epoch(const t_gtime& epo);
        /**
        * @brief GNSS/sat filter.
        * @param[in]  prn        satellite prn
        * @return
        *        true:the epoch fits gnss & sat
        *        false:the epoch does not fit gnss & sat
        */
        virtual bool _filter_gnss(const string& prn);
        /**
        * @brief get single line from the buffer.
        * @param[in]  str        the content of the single line
        * @param[in]  from_pos    the position in the buffer
        * @return      int
        */
        int _getline(string& str, int from_pos = 0);
        /**
        * @brief get the buffer.
        * @param[in]  buff        buffer
        * @return      int
        */
        int _getbuffer(const char*& buff);
        /**
        * @brief get the buffer.
        * @param[in]  buff        buffer
        * @return      int
        */
        int _add2buffer( char* buff, int sz );
        /**
        * @brief remove from buffer.
        * @param[in]  bytes_to_eat    int
        * @return      int
        */
        int _consume( int bytes_to_eat );

        weak_ptr<t_gio>         _gio_ptr;
    // vector<string>          _err;                 // cummative error message default systems
        vector<t_gnote>         _notes;               // cummative notes message/warning/error
        string                  _fname;               // decoded file
        string                  _class;               // string for reporting
        bool                    _initialized;         // if initialized
        bool                    _gnss;                // if gnss definition is requested
        int                     _out_len;             // [min] encoder data batch length
        float                   _out_smp;             // [sec] encoder data batch sample
        t_gtime                 _out_epo;             // encoder refrence epoch
        string                  _version;             // format version
        string                  _initver;             // format initial version
        map<string,t_gdata*>    _data;                // data pointer
        t_glog*                 _log;                 // log pointer
        t_gsetbase*             _set;                 // set pointer
        char*                   _buffer;              // class buffer
        int                     _buffsz;              // size of buffer
        int                     _endpos;              // after last position
        int					   _begpos;				 // begin position of useful data
        int                     _irc;
        bool                    _close_with_warning;  // close with warnings (desctructor)
   

        gcoder_char_buffer	       _decode_buffer;

        // settings
        t_gtime     _beg;                             // default beg time
        t_gtime     _end;                             // default end time
        double      _int;                             // default interval
        int         _scl;                             // default scaling for decimation (if >= 1Hz)
        set<string> _sys;                             // default systems
        set<string> _rec;                             // default sites/receivers
   
        map<GSYS,set<string> > _sat;                   // default satellites
        map<GSYS,set<string> > _obs;                   // default observations (all mixed for single GNSS !!)
        map<GSYS,set<string> > _nav;                   // default navigation messages

      // ENCODING
        int _fill_buffer( char* buff, int sz );
        stringstream     _ss;
        long             _ss_position;
        string           _pgm;
        bool             _hdr;
   
    #ifdef BMUTEX
       boost::mutex            _mutex;               // buffer mutex
    #else
      t_gmutex                 _mutex;               // buffer mutex
    #endif
     private:

    };

} // namespace

#endif
