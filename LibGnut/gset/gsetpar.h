/**
*

* @verbatim
  The format of this block:
	<parameters>
	     <GEO   sigCX    = "0.001"  sigCY    = "0.001" sigCZ  = "0.001"/>   
		 <ERP   sigXPOLE = "0.300"  sigYPOLE = "0.030" sigUT1 = "0.030"  sigDXPOLE = "0.0001" sigDYPOLE = "0.002"  sigDUT1 = "0.002"/>     
		 <STA   ID = "AIRA" sigCLK = "9000" sigZTD = "0.201"    sigION = "1.000"     sigPOS = "0.1_0.1_0.1"/>
		 <SAT   ID = "G01"  sigCLK = "5000" sigPOS = "10_10_10" sigVEL = "10_10_10"  sigECOM = "0.1_0.1_0.1_0.1_0.1_0.1_0.1_0.1_0.1"/>
	</parameters>
  @endverbatim

* @verbatim
	History
	 -1.0 jdhuang	 2019-04-07 creat the file.
	 -1.1 jdhuang    2019-05-12 add the ion par for site.
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file			gsetpar.h
* @brief		set pars from XML
*
* @author       jdhuang, Wuhan University
* @version		1.0.0
* @date			2019-04-07
*
*/

#ifndef GSETPARS_H
#define GSETPARS_H

#define XMLKEY_PARS        "parameters"
#define XMLKEY_PARS_ERP    "ERP"
#define XMLKEY_PARS_GEO    "GEO"
#define XMLKEY_PARS_SAT    "SAT"
#define XMLKEY_PARS_STA    "STA"
#define XMLKEY_PARS_AMB    "AMB"
#define XMLKEY_PARS_SRP    "SRP"

#include "gset/gsetbase.h"
#include "gexport/ExportLibGnut.h"
#include "gutils/gtriple.h"
using namespace gnut;
namespace great
{
	class LibGnut_LIBRARY_EXPORT t_gsetpar : public virtual t_gsetbase
	{
	public:
		t_gsetpar();
		virtual ~t_gsetpar();

		double sigCX();
		double sigCY();
		double sigCZ();

		double sigXpole();
		double sigYpole();
		double sigDxpole();
		double sigDypole();
		double sigUt1();
		double sigDut1();

		
		double sigAmb();
		double sigZtd(string rec);
		double sigSion(string rec);
		double sigVion(string rec);
		double sigTropPd(string rec);
		double sigIonoPd(string rec);
		double sigGRD(string rec);
		double sigGrdPd(string rec);

		double sigRecCLK(string rec);
		double sigSatCLK(string sat);
		double sigRB();

		std::map<std::string, double> sigRecPos(string sta_name);
		std::map<std::string, double> sigSatPos(string sat_name);
		std::map<std::string, double> sigSatVel(string sat_name);
		std::map<std::string, double> sigSatEcom(string sat_name);
		std::map<std::string, double> sigSatAbw(string sat_name);


		t_gtriple   sigRecPosTriple(string sta_name);


		void help();
		void check() {};
	private:

		// get child value
		string _child_value(const string& child);
		string _child_attribute_value(const string& child, const string& attribute);
		string _attribute_value(const string & index_name, const string & index_value, const string & attribute);
		void  split(const std::string& s, std::string delim, std::vector< std::string >& ret);
		// default settings
		
		// for GEO
		const double _default_sigCX = 0.001;
		const double _default_sigCY = 0.001;
		const double _default_sigCZ = 0.001;

		// for ERP
		const double _default_sigXPOLE = 0.300;
		const double _default_sigYPOLE = 0.300;
		const double _default_sigDXPOLE = 0.030;
		const double _default_sigDYPOLE = 0.030;
		const double _default_sigUT1  = 1E-4;
		const double _default_sigDUT1 = 2E-3;

		// for ZTD/ION
		const double _default_sigZTD = 0.201;
		const double _default_sigSION   = 9000;  // 1.000 -> 9000
		const double _default_sigGpsIfb = 9000;
		const double _default_sigVION   = 9000;  // 1.000 -> 9000
		const double _default_sigTropPd = 0.015;
		const double _default_sigIonoPd = 0.15;
		const double _default_sigGRD = 0.001;
		const double _default_sigGrdPd = 0.0003;

		// for AMB
		const double _default_sigAMB = 10000;

		// for CLK
		const double _default_sta_sigCLK = 9000;
		const double _default_sat_sigCLK = 5000;
		const double _default_sta_sigRB = 0.1;

		// for POS
		const double _default_sta_sigPOS = 0.1;
		const double _default_sat_sigPOS = 10;


		// for VEL
		const double _default_sat_sigVEL = 0.10;

		// for solar pars
		const double _default_sat_sigECOM = 0.10;
		const double _default_sat_sigABW_SY = 1E20;
		const double _default_sat_sigABW_AD = 0.01;
		const double _default_sat_sigABW_R = 0.1;
		const double _default_leo_sigDYN = 0.10;
	};
}

#endif