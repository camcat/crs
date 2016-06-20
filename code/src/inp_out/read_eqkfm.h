
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


#ifndef READ_EQKFM_H_
#define READ_EQKFM_H_

#include <math.h>
#include <stdio.h>

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "../general/mem_mgmt.h"
#include "../general/setup.h"
#include "../general/struct_conversions.h"
#include "../seis/soumod1.h"
#include "../seis/WellsCoppersmith.h"
#include "../util/files.h"

#include "../util/util1.h"

int eqkfm_addslipmodels(struct eqkfm *eqfm1, struct slipmodels_list all_slipmodels,
						struct eqkfm **eqfm_comb, int N1,
						int *Ncomb, int **nfout, double dt, double dmag, double res,
						struct crust crst, struct flags flags);

int focmec2slipmodel(struct crust crst, struct eqkfm *eqfm1, double res, int refine, int taper);
int read_eqkfm(char *fname, char *cmbformat, struct eqkfm **eqfm1, int *NF_out, double *Mw, double mu);
int read_farfalle_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF_out);
int read_pscmp_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF2);

#endif /* READ_EQKFM_H_ */
