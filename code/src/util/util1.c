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

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "util1.h"
#include "../defines.h"
#include "error.h"
#include <gsl/gsl_rng.h>

#define offset 1


int *iarray(long x1, long xN)
/* Returns an int vector: v[x1..xN] */
{
	int *v;

	v=(int *)malloc((size_t) ((xN-x1+1+offset)*sizeof(int)));
	if (!v) error_quit("Could not allocate memory (iarray)\n");
	return v-x1+offset;
}


unsigned long *larray(long x1, long xN)
/* Returns an unsigned long vector: v[x1..xN] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((xN-x1+1+offset)*sizeof(long)));
	if (!v) error_quit("Could not allocate memory (larray)\n");
	return v-x1+offset;
}

double *darray(long x1, long xN)
/* Returns a double vector: v[x1..xN] */
{
	double *v;

	v=(double *)malloc((size_t) ((xN-x1+1+offset)*sizeof(double)));
	if (!v) error_quit("Could not allocate memory (darray)\n");
	return v-x1+offset;
}

float **f2array(long x1, long xN, long y1, long yN)
/* Returns a float matrix: m[x1..xN][y1..yN] */
{
	long i, nx=xN-x1+1,ny=yN-y1+1;
	float **m;

	/* Set pointers to rows */
	m=(float **) malloc((size_t)((nx+offset)*sizeof(float*)));
	if (!m) error_quit("Could not allocate memory (matrix)\n");
	m += offset;
	m -= x1;

	/* Allocate rows and set pointers */
	m[x1]=(float *) malloc((size_t)((nx*ny+offset)*sizeof(float)));
	if (!m[x1]) error_quit("Could not allocate memory (matrix)\n");
	m[x1] += offset;
	m[x1] -= y1;

	for(i=x1+1;i<=xN;i++) m[i]=m[i-1]+ny;

	/* return pointer to array of pointers to rows */
	return m;
}

double **d2array(long x1, long xN, long y1, long yN)
/* Returns a double matrix: m[x1..xN][y1..yN] */
{
	long i, nx=xN-x1+1,ny=yN-y1+1;
	double **m;

	/* Set pointers to rows */
	m=(double **) malloc((size_t)((nx+offset)*sizeof(double*)));
	if (!m) {
		error_noquit("Could not allocate memory (matrix)\n");
		return NULL;
	}
	m += offset;
	m -= x1;

	/* Allocate rows and set pointers */
	m[x1]=(double *) malloc((size_t)((nx*ny+offset)*sizeof(double)));
	if (!m[x1]) {
		error_noquit("Could not allocate memory (matrix)\n");
		return NULL;
	}
	m[x1] += offset;
	m[x1] -= y1;

	for(i=x1+1;i<=xN;i++) m[i]=m[i-1]+ny;

	/* return pointer to array of pointers to rows */
	return m;
}

double ***d3array(long x1, long xN, long y1, long yN, long z1, long zN)
/* Returns a double 3tensor with range t[x1..xN][y1..yN][z1..zN] */
{
	long i,j,nx=xN-x1+1,ny=yN-y1+1,nz=zN-z1+1;
	double ***t;

	/* Set pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nx+offset)*sizeof(double**)));
	if (!t) error_quit("Could not allocate memory (f3array)\n");
	t += offset;
	t -= x1;

	/* Set pointers to rows and set pointers */
	t[x1]=(double **) malloc((size_t)((nx*ny+offset)*sizeof(double*)));
	if (!t[x1]) error_quit("Could not allocate memory (f3array)\n");
	t[x1] += offset;
	t[x1] -= y1;

	/* Allocate rows and set pointers */
	t[x1][y1]=(double *) malloc((size_t)((nx*ny*nz+offset)*sizeof(double)));
	if (!t[x1][y1]) error_quit("Could not allocate memory (f3array)\n");
	t[x1][y1] += offset;
	t[x1][y1] -= z1;

	for(j=y1+1;j<=yN;j++) t[x1][j]=t[x1][j-1]+nz;
	for(i=x1+1;i<=xN;i++) {
		t[i]=t[i-1]+ny;
		t[i][y1]=t[i-1][y1]+ny*nz;
		for(j=y1+1;j<=yN;j++) t[i][j]=t[i][j-1]+nz;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

int **i2array(long x1, long xN, long y1, long yN)
/* Returns a int matrix: m[x1..xN][y1..yN] */
{
	long i, nx=xN-x1+1,ny=yN-y1+1;
	int **m;

	/* Set pointers to rows */
	m=(int **) malloc((size_t)((nx+offset)*sizeof(int*)));
	if (!m) error_quit("Could not allocate memory (matrix)\n");
	m += offset;
	m -= x1;


	/* Allocate rows and set pointers */
	m[x1]=(int *) malloc((size_t)((nx*ny+offset)*sizeof(int)));
	if (!m[x1]) error_quit("Could not allocate memory (matrix)\n");
	m[x1] += offset;
	m[x1] -= y1;

	for(i=x1+1;i<=xN;i++) m[i]=m[i-1]+ny;

	/* return pointer to array of pointers to rows */
	return m;
}

 
float ***f3array(long x1, long xN, long y1, long yN, long z1, long zN){
/* Returns a float 3tensor with range t[x1..xN][y1..yN][z1..zN] */
 
	long i,j,nx=xN-x1+1,ny=yN-y1+1,nz=zN-z1+1;
	float ***t;

	/* Set pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nx+offset)*sizeof(float**)));
	if (!t) error_quit("Could not allocate memory (f3array)\n");
	t += offset;
	t -= x1;

	/* Set pointers to rows and set pointers */
	t[x1]=(float **) malloc((size_t)((nx*ny+offset)*sizeof(float*)));
	if (!t[x1]) error_quit("Could not allocate memory (f3array)\n");
	t[x1] += offset;
	t[x1] -= y1;

	/* Allocate rows and set pointers */
	t[x1][y1]=(float *) malloc((size_t)((nx*ny*nz+offset)*sizeof(float)));
	if (!t[x1][y1]) error_quit("Could not allocate memory (f3array)\n");
	t[x1][y1] += offset;
	t[x1][y1] -= z1;

	for(j=y1+1;j<=yN;j++) t[x1][j]=t[x1][j-1]+nz;
	for(i=x1+1;i<=xN;i++) {
		t[i]=t[i-1]+ny;
		t[i][y1]=t[i-1][y1]+ny*nz;
		for(j=y1+1;j<=yN;j++) t[i][j]=t[i][j-1]+nz;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

int **i2array_firstlevel(long xN)
/* Returns a int matrix: m[0..xN]; columns will be allocated later. */
{
	int **m;

	/* Set pointers to rows */
	m=(int **) malloc((size_t)((xN+1)*sizeof(int*)));

	/* set rows to NULL (so that later it will be recognized that memory should be allocated).*/
	for (int i=0; i<=xN; i++) m[i]=NULL;

	/* return pointer to array of pointers to rows */
	return m;
}

double **d2array_firstlevel(long xN)
/* Returns a double matrix: m[x1..xN];  columns will be allocated later. */
{
	double **m;

	/* Set pointers to rows */
	m=(double **) malloc((size_t)((xN+1)*sizeof(double*)));

	/* set rows to NULL (so that later it will be recognized that memory should be allocated).*/
	for (int i=0; i<=xN; i++) m[i]=NULL;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_i2array_firstlevel(int **m, long x1, long xN, long y1, long yN)
/* free an int matrix allocated by i2array_firstlevel() */
{
	for (int i=x1; i<=xN; i++) if (m[i]) free_iarray(m[i], y1, yN);
	free((char *) (m+x1-offset));
}

void free_d2array_firstlevel(double **m, long x1, long xN, long y1, long yN)
/* free an int matrix allocated by d2array_firstlevel() */
{
	for (int i=x1; i<=xN; i++) if (m[i]) free_darray(m[i], y1, yN);
	free((char *) (m+x1-offset));
}
void free_array(float *v, long x1, long xN)
/* free a float vector allocated with vector() */
{
	free((char *) (v+x1-offset));
}

void free_iarray(int *v, long x1, long xN)
/* free an int vector allocated with iarray() */
{
	free((char *) (v+x1-offset));
}


void free_larray(unsigned long *v, long x1, long xN)
/* free an unsigned long vector allocated with larray() */
{
	free((char *) (v+x1-offset));
}

void free_darray(double *v, long x1, long xN)
/* free a double vector allocated with darray() */
{
	free((char *) (v+x1-offset));
}

void free_f2array(float **m, long x1, long xN, long y1, long yN)
/* free a float matrix allocated by f2array() */
{
	free((char *) (m[x1]+y1-offset));
	free((char *) (m+x1-offset));
}

void free_d2array(double **m, long x1, long xN, long y1, long yN)
/* free a double matrix allocated by d2array() */
{
	free((char *) (m[x1]+y1-offset));
	free((char *) (m+x1-offset));
}

void free_i2array(int **m, long x1, long xN, long y1, long yN)
/* free an int matrix allocated by i2array() */
{
	free((char *) (m[x1]+y1-offset));
	free((char *) (m+x1-offset));
}

void free_f3array(float ***t, long x1, long xN, long y1, long yN,
	long z1, long zN)
/* free a float f3array allocated by f3tensor() */
{
	free((char *) (t[x1][y1]+z1-offset));
	free((char *) (t[x1]+y1-offset));
	free((char *) (t+x1-offset));
}

void free_d3array(double ***t, long x1, long xN, long y1, long yN,
	long z1, long zN)
/* free a double f3array allocated by f3array() */
{
	free((char *) (t[x1][y1]+z1-offset));
	free((char *) (t[x1]+y1-offset));
	free((char *) (t+x1-offset));
}


float ran1(){

  static int first_entry=1;

  if (first_entry){
	gsl_rng_set (global_rand, -1*global_seed);
	first_entry=0;
  }

  return gsl_rng_uniform_pos (global_rand);

}

