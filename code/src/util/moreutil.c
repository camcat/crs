
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


#include "moreutil.h"
#include "util1.h"
#include "../defines.h"
//#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define NSTACK 50
#define FREE_ARG char*
#define NR_END 1

int search_string(char *str1, char *str2){
/* Check if str2 is contained in str1.
 */

	int nchar=strlen(str2), off=0;

	while(off+nchar<=strlen(str1)){
		if (strncmp(str1+off,str2,nchar)==0) {
			return 1;
		}
		off++;
	}

	return 0;
}

int closest_element(double *v, int N, double value, double toll){
	/*
	 * Returns index of the element in v closes to value; returns -1 if the difference is larger than toll.
	 * v range [0...N-1]
	 */

	int i;
	double mindist=1e30;

	for (int n=0; n<N; n++){
		if(fabs(v[n]-value)<mindist){
			mindist=fabs(v[n]-value);
			i=n;
		}
	}

	if (mindist>toll) return -1;
	else return i;
}

int *nth_index(int i, int Ndim, int *dim){
//given linear index and N dimensions with indices [1...dim[0]], [1...dim[1]], ..., [1...dim[Ndim-1]], returns array with indices of each dimension.

	int f=1;
	int *res=iarray(0,Ndim-1);

	for (int d=0; d<Ndim; d++) {
		res[d]=((i-1)%(f*dim[d]))/f+1;
		f*=dim[d];
	}

	return res;
}


void copy_matrix( double **m1, double ***m2, int a, int b){
	/* m1: original matrix; *m2: new matrix. indices: [1...a, 1...b]
	 * if *m2==NULL, memory allocated. Otherwise, m2 must have correct no. of elements.
	 */

	if (!(*m2)) {
		*m2=d2array(1,a,1,b);
		if (!(*m2)) memory_error_quit;

	}

	for (int ns=1; ns<=a; ns++){
		for (int n=1; n<=b; n++) (*m2)[ns][n]=m1[ns][n];
	}

	return;
}

void copy_vector(double *m1, double **m2, int a){
	/* m1: original vector; *m2: new vector. indices: [1...a]
	 * if *m2==NULL, memory allocated. Otherwise, m2 must have correct no. of elements.
	 */

	if (!(*m2)) *m2=darray(1,a);

	for (int ns=1; ns<=a; ns++){
		(*m2)[ns]=m1[ns];
	}

	return;
}

void mysort(unsigned long n, double *old_arr, int **ind_out, double **arr_out){
/* 
*  old_arr[1...n],
*  *ind_out[1...n],
*  *arr_out[1...n].
*/
	int *ind= (ind_out)? *ind_out : NULL;
	double *arr= *arr_out;
	struct value_index *xi;

	if (!ind) ind=iarray(1,n);
	if (!arr) arr=darray(1,n);

	xi=(struct value_index *) malloc(sizeof(struct value_index) * n);

	for (int i=0; i<n; i++){
		xi[i].value=old_arr[i+1];
		xi[i].index=i+1;	//indices in old_arr
	}

	qsort(xi, n, sizeof(struct value_index), compare_valueindex);

	for (int i=0; i<n; i++){
		arr[i+1]=xi[i].value;
		ind[i+1]=xi[i].index;	//indices in old_arr
	}

	free(xi);

        if (ind_out) *ind_out=ind;
        if (arr_out) *arr_out=arr;

}


int compare_valueindex (const void * a, const void * b)
{

  struct value_index *A= (struct value_index *) a;
  struct value_index *B= (struct value_index *) b;

  if ( A->value  <  B->value ) return -1;
  if ( A->value  == B->value ) return 0;
  if ( A->value  >  B->value ) return 1;
}


double **mtimesm3(double **m1, double **m2, double ***m30){
//if m30==NULL, m2 will be ignored and new memory wll be allocated; otherwise, m30 is filled (NB memory needs to have been allocated previously!).

//indices: [1...3].

	double **m3;

	m3= (m30) ? *m30 : d2array(1,3,1,3);

	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++){
			m3[i][j]=0;
			for (int k=1; k<=3; k++){
				m3[i][j]+=m1[i][k]*m2[k][j];
			}
		}
	}

	return m3;
}

double *mtimesv(double **M, double *v, double *v2, int D1, int D2){
//if v2==NULL, v2 will be ignored and new memory wll be allocated; otherwise, vector v2 is filled (NB memory needs to have been allocated previously!).

	int i, j;
	double temp;
	double *v2int;

	v2int=(v2)? v2 : darray(1,D2);

	for (i=1; i<=D2; i++){
		temp=0;
		for (j=1;j<=D1; j++){
			temp+=M[i][j]*v[j];
		}
	v2int[i]=temp;
	}
	return v2int;
}

double norm(double *v1, int D){
//calculates vector norm.

	double r2;
	r2=vdotv(v1,v1,D);

	return sqrt(r2);
}

double vdotv(double *v1, double *v2, int D){
//dot product.

	int k;
	double temp2=0;
	for (k=1; k<=D; k++){
		temp2+=v1[k]*v2[k];
	}
	return temp2;
}

void normv (double *v, int D){
// normalizes vector v (so that the sum of its elements is 1).

	double sum=0.0;
	for (int k=1; k<=D; k++) sum+=v[k];
	for (int k=1; k<=D; k++) v[k]*=(1.0/sum);

	return;
}

void unitv (double *v, int D){
// normalizes vector v (so that its length is 1).

	double length2;
	length2=vdotv(v, v, D);
	for (int k=1; k<=D; k++){
			v[k]=v[k]/pow(length2,0.5);
		}
	return;
}


void nearest_neighbours(int NP, int D1, int D2, int D3, int **nn){
// points change along D1, then along D2, then along D3 (DX is no of points in each dimension).
	int P[4], P2[4];
	int Ds[4]={0, D1, D2, D3};
	int D1D2=D1*D2;
	int pt;

	for (int i=1; i<=NP; i++){
		//reshape linear array into 3x3 array: i -> (P1,P2,P3).
		P[1]=(i-1)%D1+1;
		P[2]=((i-1)%D1D2)/D1 +1;
		P[3]=(i-1)/D1D2 +1;

		pt=0;
		for (int p=1; p<=3; p++)
			for (int d=-1; d<=1; d+=2){
				pt+=1;
				P2[p]=P[p]+d;
				P2[p%3+1]=P[p%3+1];
				P2[(p+1)%3+1]=P[(p+1)%3+1];

				if (P2[p]==0 || P2[p]>Ds[p]) nn[i][pt]=0;
				else {
					nn[i][pt]=P2[1]+(P2[2]-1)*D1+(P2[3]-1)*D1D2;
			}
		}
	}
}

void interp_nn(int NP, int D1, int D2, int D3, double *values, double **allvalues, int all6, int **nn0){
	//if all6==1 gives all 6 values; if all6==0, only maximum and minimum (also including central point).

	int **nn;

	if (!nn0) {
		nn=i2array(1,NP,1,6);
		nearest_neighbours(NP,D1,D2,D3,nn);
	}
	else nn=nn0;

	#pragma omp parallel for
	for (int i=1; i<=NP; i++){
		if (all6==0){
			allvalues[i][1]=values[i];
			allvalues[i][2]=values[i];
		}
		for (int v=1; v<=6; v++){
			if (nn[i][v]==0) continue;
			else {
				if (all6==0){
					allvalues[i][1]=fmin(allvalues[i][1],0.5*(values[nn[i][v]]+values[i]));
					allvalues[i][2]=fmax(allvalues[i][2],0.5*(values[nn[i][v]]+values[i]));
				}
				else allvalues[i][v]=0.5*(values[nn[i][v]]+values[i]);
			}
		}
	}

	if (nn0== (int **)0) free_i2array(nn,1,NP,1,6);
}

double * duplicate_darray(double *v, long nrl, long nrh){
	double * vnew=darray(nrl, nrh);
	for (int r=nrl; r<=nrh; r++){
		vnew[r]=v[r];
	}
	return vnew;
}

double min_v(double *v, int N){

	double min=1/0.0;
	for (int i=0; i<N; i++) if (v[i]<min) min=v[i];
	return min;
}

double max_v(double *v, int N){
//v[0...N-1];

	double max=-1/0.0;
	for (int i=0; i<N; i++) if (v[i]>max) max=v[i];
	return max;
}

#undef NSTACK
#undef SWAP
