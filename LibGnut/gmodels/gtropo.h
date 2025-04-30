
/**
*
* @verbatim
    History
    2011-01-10 /JD: created

  @endverbatim
* Copyright (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
*
* @file       gtropo.h
* @brief      Purpose: implements troposphere model class
*.
* @author     JD
* @version    1.0.0
* @date       2011-01-10
*
*/

#ifndef GTROPO_H
#define GTROPO_H 


#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gprod/gprodcrd.h"
#include "gmodels/ggpt.h"
#include "gexport/ExportLibGnut.h"

using namespace std;

namespace gnut
{
	/** @brief class for t_gtropo. */
	class LibGnut_LIBRARY_EXPORT t_gtropo
	{
	public:
		/** @brief default constructor. */
		t_gtropo();
		virtual ~t_gtropo();

		virtual double getZHD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]

	protected:
		t_gpt        _gpt;
	};

	/** @brief class for t_saast derive from t_gtropo. */
	class LibGnut_LIBRARY_EXPORT t_saast : public t_gtropo {
	public:
		/** @brief default constructor. */
		t_saast() {}
		~t_saast() {}

		virtual double getSTD(const double& ele, const double& hel);   // ! Radians: elevation
		virtual double getZHD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
	};

	/** @brief class for t_saast derive from t_gtropo. */
	class LibGnut_LIBRARY_EXPORT t_davis : public t_gtropo {
	public:
		t_davis() {}
		~t_davis() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]  
	};

	/** @brief class for t_hopf derive from t_gtropo. */
	class LibGnut_LIBRARY_EXPORT t_hopf : public t_gtropo {
	public:
		t_hopf() {}
		~t_hopf() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};

	/** @brief class for t_baby derive from t_gtropo. */
	class LibGnut_LIBRARY_EXPORT t_baby : public t_gtropo {
	public:
		t_baby() {}
		~t_baby() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};

	/** @brief class for t_chao derive from t_gtropo. */
	class LibGnut_LIBRARY_EXPORT t_chao : public t_gtropo {
	public:
		t_chao() {}
		~t_chao() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};

	/** @brief class for t_ifad derive from t_gtropo. */
	class LibGnut_LIBRARY_EXPORT t_ifad : public t_gtropo {
	public:
		t_ifad() {}
		~t_ifad() {}

		virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
		virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
	};

} // namespace

#endif
