
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


#ifndef HELMSTETTER_H_
#define HELMSTETTER_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../geom/find_gridpoints.h"

#include "../util/util1.h"

double *smoothed_rate_Helmstetter(double *xgrid, double *ygrid, double dx, double dy, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord);
double *smoothed_rate_Helmstetter_nonuni(double *xgrid, double *ygrid, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord);
double *smoothed_rate_Helmstetter_cat(struct catalog cat, struct crust crst, double *weights, int ord);
double *fit_depth(double *zgrid, double dz, int Ngrid, double *zs, double *err, double *weights, int N);

#endif /* HELMSTETTER_H_ */
