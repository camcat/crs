
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


#ifndef SETUP_TIME_H_
#define SETUP_TIME_H_

#include "../defines.h"

int timesteps_log(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2, double smallstepstime);

int timesteps_lin(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2, int ***allind);

int setup_aseismic_multi_linear(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2);

int setup_aseismic_single_linear(double t0, double t1, struct eqkfm **eqk_aft, int NA, int *Nfaults, int *L, double **times2);

int setup_aseismic_splines(double t0, double t1, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2, long *seed);

int setup_aseismic_single_log(double t0, double t1, double ts, struct eqkfm **eqk_aft,
						 int NA, int *Nfaults, int *L, double **times2, long *seed);

#endif /* SETUP_TIME_H_ */
