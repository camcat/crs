
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


#ifndef READ_ZMAP_H_
#define READ_ZMAP_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../geom/find_gridpoints.h"
#include "../seis/GR.h"

#include "../util/util1.h"
#include "read_matrix.h"

int readZMAP (struct catalog *cat, struct eqkfm **eqfm, int *, char *file, struct crust crst, struct tm reftime, double t0s, double t1s,
		double t0c, double t1c, double Mmain, double tw, double border, double, double dDCFS, int findgridpoints);

int read_firstlineZMAP(char *file, struct tm reftime, double *time);

#endif /* READ_ZMAP_H_ */
