
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


#ifndef LIN_INTERP_EQKFM_H_
#define LIN_INTERP_EQKFM_H_

void fit_lin(double *times1, double *times2, int Nas, int L, int *ind, int NP, double **slipbefore_st, double **slip_after);
void lin_interp_eqkfm(struct eqkfm **eqkfm_aft, int NF, int L, double *times2, int *ind);

#endif /* LIN_INTERP_EQKFM_H_ */
