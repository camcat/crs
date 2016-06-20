
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


#ifndef READ_TXTTEMPLATE_H_
#define READ_TXTTEMPLATE_H_
#include <math.h>
#include <stddef.h>

#include "../defines.h"
#include "../geom/convert_geometry.h"
#include "../util/moreutil.h"

#include "../util/util1.h"
#include "read_matrix.h"

int read_rate(struct crust crst, char *fname, double **bg_rate, double *r0, double *minmag);
int read_csep_template(char *fname, int *no_magbins, int *nlat, int *nlon, int *ndep, int *ng, double *dlat, double *dlon, double *ddep, double *dmag,
		double **lats, double **lons, double **deps, double **rate, double *minlat, double *maxlat, double *minlon, double *maxlon, double *mindep, double *maxdep,
		double *minmag, double *maxmag, int *uni);
#endif /* READ_TXTTEMPLATE_H_ */
