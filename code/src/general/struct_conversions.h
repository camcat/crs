
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


#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "../defines.h"
#include "../geom/coord_trafos.h"

#include "../util/util1.h"
#include "eqkfm_copy.h"
#include "mem_mgmt.h"


void eqk_filter(struct eqkfm **eqkfm1, int *Ntot, double Mag, double Depth);
int * cat_filter(struct catalog *cat, double Mag, double);
int *combine_eqkfm(struct eqkfm *eqkfm1, struct eqkfm *eqkfm2, int N1, int N2, double dt, double dM, double dR, int overwrite);
int *combine_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM);
double **union_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM, int***, int *);
double **union_cats2(struct catalog cat, struct pscmp *DCFS, int N2, int ***ind, int *tot);
double *timesfromeqkfm(struct eqkfm *eqkfm1, int N, int *);
double *magssfromeqkfm(struct eqkfm *eqkfm1, int N, int *);
double *timesfrompscmp(struct pscmp *DCFS, int N);
double *magsfrompscmp(struct pscmp *DCFS, int N);
void eqkfm2dist(struct eqkfm *eqkfm1, double *lats, double *lons, double *depths, int N, int Ntot, int);

