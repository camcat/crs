
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


#ifndef SOUMOD1_H_
#define SOUMOD1_H_

#endif /* SOUMOD1_H_ */

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "../general/mem_mgmt.h"
#include "../geom/coord_trafos.h"
#include "../inp_out/print_output.h"

#include "../util/util1.h"
//#include "nr.h"

int scale_to_mag(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double * slips, double *rakes);
int suomod1_taper(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, int top, int bottom, int right, int left);
int suomod1_resample(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double disc);

