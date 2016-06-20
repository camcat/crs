
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

//#include "../defines.h"
#include "error.h"

int search_string(char *str1, char *str2);
int scan_nth(char *string, int n, double *result);
int closest_element(double *v, int N, double value, double toll);
int *nth_index(int i, int Ndim, int *dim);
void copy_matrix( double **m1, double ***m2, int a, int b);
void copy_vector(double *m1, double **m2, int a);
char ***tmatrix(long nrl, long nrh, long ncl, long nch, long length);
void mysort(unsigned long n, double *old_arr, int **ind, double **arr);
void free_tmatrix(char ***m, long nrl, long nrh, long ncl, long nch, long length);
double **mtimesm3(double **m1, double **m2, double ***);
double * mtimesv(double **M, double *v, double *v2, int D1, int D2);
double norm(double *v1, int D);
double vdotv(double *v1, double *v2, int D);
void normv (double *v, int D);
void unitv (double *v, int D);
void statistics (double *, int, double *, double *);
double ***d3array(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3array(double ***t, long nrl, long nrh, long ncl, long nch,	long ndl, long ndh);
void intersect_lists(int *l1, int *l2, int **l3, int **, int **, int N1, int N2, int *N3);
void nearest_neighbours(int NP, int D1, int D2, int D3, int **nn);
void interp_nn(int NP, int D1, int D2, int D3, double *values, double **allvalues, int all6, int **);
double *** duplicate_d3array(double ***S, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double * duplicate_darray(double *v, long nrl, long nrh);
double min_v(double *v, int N);
double max_v(double *v, int N);
int **i2array_firstlevel(long nrh);
double **d2array_firstlevel(long nrh);
void free_i2array_firstlevel(int **m, long nrl, long nrh, long ncl, long nch);
void free_d2array_firstlevel(double **m, long nrl, long nrh, long ncl, long nch);




