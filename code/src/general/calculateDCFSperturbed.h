
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


#ifndef CALCULATEDCFSPERTURBED_H_
#define CALCULATEDCFSPERTURBED_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../okada/okadaDCFS.h"
#include "../seis/cmbopt.h"
#include "../seis/soumod1.h"
#include "../seis/WellsCoppersmith.h"
#include "../util/moreutil.h"

#include "../util/util1.h"
#include "mem_mgmt.h"

void calculateDCFSperturbed(double **DCFSrand, struct pscmp *DCFS, struct eqkfm *eqkfmAf,
							struct eqkfm *eqkfm0, struct flags flag,
							double *times, int Nmain, int NA, struct crust crst,
							struct Coeff_LinkList *AllCoeff, struct Coeff_LinkList *AllCoeff_aseis, int NTScont,
							double **focmec, int *fmzoneslim, int NFM,
							double tdata0, double tdata1,
							int refresh, int which_recfault);

void smoothen_DCFS(struct pscmp DCFS, int, int, int, int, int **);
void smoothen_vector(int NgridT, int nLat, int nLon, int nD, double *values, long *seed, int**, int);

#endif /* CALCULATEDCFSPERTURBED_H_ */
