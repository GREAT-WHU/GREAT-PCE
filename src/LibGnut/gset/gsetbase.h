/**
*
* @verbatim
	History

  @endverbatim
*
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

* @file       gsetbase.h
* @brief      implements base setting class
* @author     Jan Dousa
* @version    1.0.0
* @date       2012-10-23
*
*/

#ifndef GSETBASE_H
#define GSETBASE_H

#include <set>
#include <map>
#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "pugixml/src/pugixml.hpp"

#define XMLKEY_ROOT "config" ///< The define of ROOT node

#ifdef WIN32
#pragma warning(disable:4250)
#endif

using namespace std;
using namespace pugi;

namespace gnut
{
	/**
	* @brief The base class for xml set
	*/
	class LibGnut_LIBRARY_EXPORT t_gsetbase : public pugi::xml_node
	{
	public:
		/// defalut constructor
		t_gsetbase();
		/// defalut destructor
		virtual ~t_gsetbase();

		/**
		* @brief  read configuration from file
		* @param[in] file file name
		*/
		int read(const string& file);

		/**
		* @brief  read configuration from istream
		* @param[in] is istream name
		*/
		int read_istream(istream& is);

		/**
		* @brief  write configuration in file
		* @param[in] file file name
		*/
		int write(const string& file);

		/**
		* @brief  write configuration in ostream
		* @param[in] os istream name
		*/
		int write_ostream(ostream& os);

		/**
		* @brief  set/get glog pointer (VIRTUAL for gsetout!)
		* @param[in] l The pointer of log file
		*/
		virtual void glog(t_glog* l) { _log = l; }

		/**
		 * @brief return the glog pointer (VIRTUAL for gsetout!)
		 * @return glog pointer
		 */
		t_glog* glog()const { return _log; }

		/// application name
		virtual string pgm();

		/// application info/version
		virtual string app();

		/**
		* @brief  get some information of the application
		* @param[out] pgm program name
		* @param[out] ver program version
		* @param[out] dat program running data
		* @param[out] tim program running time
		* @param[out] rev program revision number
		* @param[out] own program owner name
		*/
		virtual void app(const string& pgm, const string& ver,
			const string& dat, const string& tim,
			const string& rev, const string& own);

		/// application usage
		virtual void usage();

		/**
		* @brief  get command line arguments
		* @param[in] argc   argument number
		* @param[in] argv[] command line array of arguments
		*/
		virtual void arg(int argc, char *argv[],
			bool add = false, bool thin = false);

		/// settings help
		virtual void help() = 0;

		/// settings check
		virtual void check() = 0;

		/// thin execution
		bool thin();

		/// output the log cache
		void print_log();

	protected:

		/**
		 * @brief get the double value of ROOT.elem.attrib
		 * @param[in] elem    element name in xml
		 * @param[in] attrib  attribute name in xml
		 * @return double : The value of ROOT.elem.attribute
		 */
		double         _dblatt(const string& elem, const string& attrib);


		/**
		* @brief get the set<string> value of ROOT.elem.subelem
		* @param[in] elem    element name in xml
		* @param[in] subelem subelem name in xml
		* @return set<string> : The value of ROOT.elem.subelem
		*/
		set<string>    _setval(const string& elem, const string& subelem);
		set<string>    _setvals(const string& elem, const string& subelem);

		/**
		* @brief get the vector<string> value of ROOT.elem.subelem
		* @param[in] elem    element name in xml
		* @param[in] subelem subelem name in xml
		* @return vector<string> : The value of ROOT.elem.subelem
		*/
		vector<string> _vecval(const string& elem, const string& subelem);

		/**
		* @brief get default child node of a xml node
		* @param[in] node    the parent node
		* @param[in] n       child node name
		* @param[in] val     child node value
		* @param[in] reset   reset or not
		* @return xml_node : The child node
		*/
		xml_node _default_node(xml_node& node, const char* n, const char* val = "", bool reset = false);

		/**
		* @brief get default node of a xml node by the attribute or reset it with v
		* @param[in] node    the parent node
		* @param[in] n       attribute name
		* @param[in] v       reset value
		* @param[in] reset   reset or not
		* @return xml_node : The child node
		*/
		void _default_attr(xml_node& node, const char* n, const string& v, bool reset = false);

		/**
		* @brief get default node of a xml node by the attribute or reset it with v
		* @param[in] node    the parent node
		* @param[in] n       attribute name
		* @param[in] v       reset value
		* @param[in] reset   reset or not
		* @return xml_node : The child node
		*/
		void _default_attr(xml_node& node, const char* n, const bool& value, bool reset = false);

		/**
		* @brief get default node of a xml node by the attribute or reset it with v
		* @param[in] node    the parent node
		* @param[in] n       attribute name
		* @param[in] v       reset value
		* @param[in] reset   reset or not
		* @return xml_node : The child node
		*/
		void _default_attr(xml_node& node, const char* n, const int& value, bool reset = false);

		/**
		* @brief get default node of a xml node by the attribute or reset it with v
		* @param[in] node    the parent node
		* @param[in] n       attribute name
		* @param[in] v       reset value
		* @param[in] reset   reset or not
		* @return xml_node : The child node
		*/
		void _default_attr(xml_node& node, const char* n, const double& value, bool reset = false);

		/// report element issues
		void _add_log(string element, string msg);

		/// update glog (mask,verb)
		virtual void _upd_glog() {};

		/// XML settings header
		virtual void help_header();

		/// XML settings footer
		virtual void help_footer();

		xml_document     _doc;        ///< root document
		xml_parse_result _irc;        ///< result status

	protected:
		t_glog*        _log;          ///< pointer to log file
		string         _name;         ///< configuration name
		string         _delimiter;    ///< delimiter for writting nodes/elements

		string         _pgm;          ///< program name
		string         _ver;          ///< program version
		string         _rev;          ///< program revision
		string         _own;          ///< program owner
		string         _dat;          ///< program data
		string         _tim;          ///< program time

		set<string>    _set;          ///< program settings
		map<string, set<string> >    _chache_log;  ///< The chache of log pointer

		mutable t_gmutex        _gmutex;           ///< For the thread

	private:

	};

} // namespace

#endif
