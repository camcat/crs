
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


#ifndef GR_H_
#define GR_H_

#include <stdlib.h>
#include <math.h>

#include "../util/util1.h"

double *assign_GRnorm(double *mags, int N, double b, int Minf);
int compare (const void * a, const void * b);
int bin_equnumber(double *v, int N, int Nbin, double **bin_c, double **norm_count);
double Mc_maxcurv(double *mags, int N);
double calculatebvalue(double *mags, int N, double Mc);

#endif /* GR_H_ */
