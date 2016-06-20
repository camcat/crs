
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


#ifndef FORECAST_STEPG_H_
#define FORECAST_STEPG_H_

#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "../defines.h"

#include "../util/util1.h"
#include "struct_conversions.h"

int rate_state_evolution(struct catalog cat, double *times, double **cmpdata, struct pscmp *DCFS, double tt0, double tt1,
		double dt, double Asig, double ta, int points[], double *NeX, double *NeT, double *Rate_end, double ** all_NeT, int N, int NTS, int Neq,
		double *gamma_init, double *back_rate, double *R, int last);

#endif /* FORECAST_STEPG_H_ */
