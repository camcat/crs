
/*   Copyright (C) 2015 by Camilla Cattania and Fahad Khalid.
 *
 *   This file is part of CRS.
 *
 *   CRS is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CRS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CRS.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef CRS_LOGLIHELIHOOD_H_
#define CRS_LOGLIHELIHOOD_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../geom/convert_geometry.h"
#include "../inp_out/print_output.h"
#include "../inp_out/write_csep_forecast.h"

#include "../util/util1.h"
#include "calculateDCFSperturbed.h"
#include "forecast_stepG.h"

int CRSforecast (double *LL, int Nsur, struct pscmp *DCFS, struct eqkfm *eqkfm_aft, struct eqkfm *eqkfm0, struct flags flags,
		struct crust crst, struct Coeff_LinkList *AllCoeff, struct Coeff_LinkList *AllCoeff_aseis, int NTScont, int Nm, int Na, int NgridT, double **focmec, int *fmzonelim, int NFM,
		struct catalog cat, double *times, double tstart, double t0, double t1, double dtstep, double Asig, double ta, double r0,
		double **all_gammas0, int multiple_input_gammas, int fromstart,
		char * print_cmb, char *print_forex, char *print_foret, char * printall_cmb, char *printall_forex, char *printall_foret, char *print_LL, int refresh);

int CRSLogLikelihood (double *LL, double *Ldum0_out, double *Nev, double *I, double *r_out, int Nsur, struct pscmp *DCFS,
		struct eqkfm *eqkfm_aft, struct eqkfm *eqkfm0, struct flags flags,
		struct crust crst, struct Coeff_LinkList *AllCoeff, struct Coeff_LinkList *AllCoeff_aseis, int NTScont, int Nm, int Na, int NgridT, double **focmec, int *fmzonelim, int NFM,
		struct catalog cat, double *times, double tstart, double tt0, double tt1, double tw, double Mag_main, double Asig, double ta, double r0, int fixr,
		double *gammas0, double **all_new_gammas, int fromstart, int refresh);

#ifdef _CRS_MPI
// Added for use in partition size calculation when using
// MPI enabled code.
inline int roundUpFrac(double x) {
	double integralPart, fractionalPart, result;

	fractionalPart = modf(x, &integralPart);

	if(integralPart != 0.0 && fractionalPart != 0.0) {
		result = ceil(integralPart + 0.5);
	}
	else {
		result = integralPart;
	}

	return result;
}
#endif


#endif /* CRS_LOGLIHELIHOOD_H_ */
