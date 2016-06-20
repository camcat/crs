
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


#ifndef READ_INPUTFILE_H_
#define READ_INPUTFILE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../defines.h"
#include "../util/error.h"

#include "../util/util1.h"

int read_inputfile(char *input_fname, char *outname, char *fore_template,
		char *catname, char ***focmeccat, char *background_rate_file, char *background_rate_cat, char *fixedmecfile, char *slipmodelfile, char *afterslipmodelfile,
		char *model_parameters_file, char *Logfile, struct tm *reftime,
		double *Tstart, double *Tend, double *tstartLL, double *tendLL, long *seed, int *num_fm);

int read_slipformecfiles(char *inputfile, char ***listfiles, int *nfiles);
int read_listslipmodel(char *input_fname, struct tm reftime, struct slipmodels_list *allslipmodels, double res, int is_afterslip, int *aseismic_log, double *t0log, int *flag_multisnap);

#endif /* READ_INPUTFILE_H_ */
