
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


#ifndef OKADADCFS_H_
#define OKADADCFS_H_

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../defines.h"
#include "../util/moreutil.h"

#include "../util/util1.h"
#include "pscokada.h"

#define resolve_DCFS(...) resolve_DCFS0(__VA_ARGS__, 1)
#define resolve_DCFS_nocap(...) resolve_DCFS0(__VA_ARGS__, 0)
//----------------top level functions----------------//
//----------------top level functions----------------//
int resolve_DCFS0(struct pscmp DCFS, struct crust crst, double *strikeRs, double *dipRs, double *rake, int optrake, int cap);
int okadaCoeff(float ****Coeffs_st, float ****Coeffs_dip, float ****Coeffs_open, struct eqkfm *eqkfm1, int NF, struct crust crst);
int okadaCoeff_mpi(float ****Coeffs_st, float ****Coeffs_dip, float ****Coeffs_open, struct eqkfm *eqkfm1, int NF, struct crust crst);
int okadaCoeff2DCFS(float ***Coeffs_st, float ***Coeffs_d,float ***Coeffs_open, struct pscmp DCFS, struct eqkfm *eqkfm1);
int isoDCFS(struct pscmp DCFS, struct eqkfm eqkfm1);

//----------------auxiliary functions----------------//
int choose_focmec(struct eqkfm eqkfm1, double *strike, double *dip, double *rake);
void patch_pos(struct eqkfm eqfm, int p, double *east, double *north, double *depth);
double *normal_vector(double, double);
double *slip_vector(double strikeR, double dipR, double rakeR);
double *opt_s(double *stress, double sigma, double *n, double *result);
double *sum_v(double *v1, double *v2, double *sum, int N);
double resolve_S(double **S, double strikeR, double dipR, double rakeR, double f, double *stress0, double sigma0, int opt_rake);
double resolve_n(double **S, double *n, double fric, double *stress0, double sigma0, double *slip_v);

#endif /* OKADADCFS_H_ */
