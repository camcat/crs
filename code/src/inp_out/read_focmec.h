
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
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../geom/find_gridpoints.h"
#include "../seis/WellsCoppersmith.h"

#include "../util/util1.h"
#include "../general/eqkfm_copy.h"
#include "read_matrix.h"
#include "read_csv.h"

int readmultiplefocmec(char **focmecfiles, int nofiles,
					   struct crust crst, double, double, double dDCFS,
					   struct tm reftime, double t0, double t1, double tfocmec,
					   double mag, double ***focmec, int **firstelements,
					   int *NFM, int *NFM_timesel, struct eqkfm **eqkfm,
					   int sel, int fm2);

int readfocmec(char *focmecfile, struct crust crst, double, double, double dDCFS, struct tm reftime,
		double t0, double t1, double tfocmec, double mag, double ***focmec,	int *NFM, int *NFM_timesel, struct eqkfm **eqkfm,int sel, int fm2);

void select_fm_time(double **focmec, int *NFM, double Tstart);
