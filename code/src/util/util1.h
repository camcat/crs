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

int compare_valueindex (const void * a, const void * b);
float ran1();
int *iarray(long nl, long nh);
double *darray(long nl, long nh);
unsigned long *larray(long nl, long nh);
double **d2array(long nrl, long nrh, long ncl, long nch);
float **f2array(long nrl, long nrh, long ncl, long nch);
int **i2array(long nrl, long nrh, long ncl, long nch);
float ***f3array(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_iarray(int *v, long nl, long nh);
void free_darray(double *v, long nl, long nh);
void free_larray(unsigned long *v, long nl, long nh);
void free_d2array(double **m, long nrl, long nrh, long ncl, long nch);
void free_f2array(float **m, long nrl, long nrh, long ncl, long nch);
void free_i2array(int **m, long nrl, long nrh, long ncl, long nch);
void free_f3array(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
