
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


#define DAY2SEC 86400;

#include <stdio.h>
#include <time.h>

#include "../defines.h"

#include "../util/util1.h"

int sum_DCFS(struct pscmp *DCFS, double **cmb, int N, int Ntot);
int sum_DCFSrand(double **DCFSrand, double **cmb, int TS, int N);
int print_rate(char *fname, struct crust crst, double Mc, double *rate);
int print_grid(char *fname, struct pscmp DCFS, struct crust, double *rate);
int print_slipmodel(char* filename, struct eqkfm *eqfm1, int NF);
int print_cat(char *fname, struct catalog cat);
