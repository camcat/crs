
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


#ifndef CMBOPT_H_
#define CMBOPT_H_


#endif /* CMBOPT_H_ */

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../okada/okadaDCFS.h"
#include "../util/mscorr.h"

#include "../util/util1.h"
#include "../util/roots3.h"

void cmbopt(double sxx, double syy, double szz, double sxy, double syz, double szx, double p, double f, double st0, double di0, double ra0, double *cmb,double *st1, double *di1, double *ra1, double *st2, double *di2, double *ra2);
void DCFScmbopt(struct pscmp *DCFS, int, struct crust crst);
