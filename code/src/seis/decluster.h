
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


#ifndef DECLUSTER_H_
#define DECLUSTER_H_

#include <math.h>

#include "../defines.h"

#include "../util/util1.h"

int *decluster_catalog_rescalegrid(struct catalog cat, struct crust crst, double Mmain, double **time_missing, int d3);
int * decluster_catalog(struct catalog cat, double Mmain, double **time_missing, int d3);
void KG74(double M, double *D, double *T);

#endif /* DECLUSTER_H_ */
