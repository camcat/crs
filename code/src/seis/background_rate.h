
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


#ifndef BACKGROUND_RATE_H_
#define BACKGROUND_RATE_H_

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "../defines.h"
#include "../general/mem_mgmt.h"
#include "../inp_out/read_matrix.h"
#include "../inp_out/read_zmap.h"
#include "../util/moreutil.h"

#include "../util/util1.h"
#include "decluster.h"
#include "GR.h"
#include "smoothed_rate_Helmstetter.h"

int background_rate(char *catfile, struct crust *crst_in, struct tm reftime,
		double Mmain, double *minmag, double *rate, double **rate_grid, double dR, double dZ, double min_smoothing, int ord);

#endif /* BACKGROUND_RATE_H_ */
