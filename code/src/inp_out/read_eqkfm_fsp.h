
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

int read_fsp_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF_out);
void track_position(long *pos, int NP, FILE* fin);
int next_separator(FILE * fin, char *string);
int find_key(FILE *fin, char *string, double *value);
int scan_nth(char *string, int n, double *result);
int read_slipvalues(FILE *fin, struct eqkfm *eqfm);


#endif /* READ_EQKFM_H_ */
