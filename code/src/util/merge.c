
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


#include <stddef.h>


#include "../util/util1.h"
#include "../util/moreutil.h"
#include "merge.h"
#include "../defines.h"


void merge_multiple(double **vs, int *lens, int N, double **sorted, int *len_fin, int ***indices) {
/* merge N lists.
 * vs[0...N-1][..], where vs[n] is a list: vs[n][0...l-1], where l=lens[n].
 * lens[0...N-1] contains the lengths of the arrays.
 * indices[n] contains the list of elements which originally belonged to vs[n].
 *
 * It merges recursively, i.e. by merging the first two, adding the third, then the fourth, etc.
 * Probably not the most efficient way (may be faster to loop over the arrays to find the next element in the final array), but easier to implement.
 * Indices is also found in a silly way, since this function was originally output boolean values instead (ind_bool variable).
 */

 int lentot=0, lenmax=0;
 double *tempsorted0=NULL, *tempsorted1;
 int **tempind, **tempind2, **indbool;
 int templen, c;

 for (int n=0; n<N; n++) {
	 lentot+=lens[n];
	 lenmax=MAX(lenmax,lens[n]);
 }

 if (indices && !(*indices)) *indices=i2array(0,N-1,0,lenmax-1);
 if (sorted && !(*sorted)) *sorted=darray(0,lentot-1);


 if (N==1){	//a single episode
	 if (indices) for (int l=0; l<lenmax; l++) (*indices)[0][l]=l;
	 if (sorted) for (int l=0; l<lenmax; l++) (*sorted)[l]=vs[0][l];
	 *len_fin=lens[0];
	 return;
 }


 indbool=i2array(0,N-1,0,lentot-1);
 tempsorted1= (sorted)? *sorted : darray(0,lentot-1);
 tempsorted0= darray(0,lentot-1);
 tempind=i2array(0,N-1,0,lentot-1);
 tempind2=i2array(0,N-1,0,lentot-1);

 templen=0;	//length of temporary array.
 for (int n=0; n<N; n++) {
	 //tempind[n] contains the indices of the elements of tempsorted1 which were already in the previous array.
	 //tempind[n] has length templen+vs[n], i.s the templen value at the end of the loop.
	 merge(tempsorted0, templen, vs[n], lens[n], &tempsorted1, &(tempind[n]), NULL);
	 templen+=lens[n];
	 tempsorted0-=1;	//because of function copy_vector.
	 copy_vector(tempsorted1-1, &tempsorted0, templen);
	 tempsorted0+=1;
 }

 *len_fin=templen;

 //just copy last element into new tempind2:
 for (int j=0; j<*len_fin; j++) tempind2[N-1][j]= tempind[N-1][j];

 //"stretch" tempind arrays to cover final length and copy them into tempind:
 for (int n=N-2; n>=0; n--) {
	 c=0;
	 for (int j=0; j<*len_fin; j++){
		 if (tempind2[n+1][j]==1){
			 tempind2[n][j]= tempind[n][c];
			 c+=1;
		 }
		 else tempind2[n][j]=0;
	 }
 }

 //elements of array nth are those which are 1 in tempind2[n+1] and 0 in tempind2[n].
 if (indices){
	 for (int j=0; j<*len_fin; j++) indbool[N-1][j]= tempind2[N-1][j]==0;
	 for (int n=N-2; n>=0; n--) {
		 for (int j=0; j<*len_fin; j++){
			indbool[n][j]= (tempind2[n+1][j]==1) && (tempind2[n][j]==0);
		 }
	 }

	 for (int n=N-1; n>=0; n--){
		 c=0;
		 for (int j=0; j<lenmax; j++) (*indices)[n][j]=-2;	//later will shift up by 1 and become -1.
		 for (int j=0; j<*len_fin; j++){
			 if (c>=lens[n]) break;
			 if (indbool[n][j]){
				 (*indices)[n][c]=j;
				 c+=1;
			 }
		 }
	 }


 }
 free_darray(tempsorted0,0,lentot-1);
 free_i2array(tempind,0,N-1,0,lentot-1);
 free_i2array(tempind2,0,N-1,0,lentot-1);

}

void merge(double a[], int m, double b[], int n, double **sorted, int **ai, int **bi) {
/* merge 2 arrays: a[0...m-1] b[0...n-1].
 * sorted is the sorted array.
 * ai, bi[0...n+m-1] contain 0s and 1s indicating whether the element in sorted belonged to a or b.
 * *sorted will be allocated if NULL; *ai (*bi) are also allocated if *ai==NULL and ai!=NULL.
 */

  int i, j, k;

  //allocate arrays if needed:
  if (!(*sorted)) *sorted=darray(0,m+n-1);
  if (ai && !(*ai)) *ai=iarray(0,m+n-1);
  if (bi && !(*bi)) *bi=iarray(0,m+n-1);

  j = k = 0;

  for (i = 0; i < m + n; i++) {
	if (ai) (*ai)[i]=0;
	if (bi) (*bi)[i]=0;
  }

  for (i = 0; i < m + n;) {

    if (j < m && k < n) {
      if (a[j] < b[k]) {
        (*sorted)[i] = a[j];
        if (ai) (*ai)[i]=1;
        j++;
      }
      else {
    	(*sorted)[i] = b[k];
        if (bi) (*bi)[i]=1;
        k++;
        }
      i++;
    }
    else if (j == m) {
      for (; i < m + n;) {
    	(*sorted)[i] = b[k];
        if (bi) (*bi)[i]=1;
        k++;
        i++;
      }
    }
    else {
      for (; i < m + n;) {
        (*sorted)[i] = a[j];
        if (ai) (*ai)[i]=1;
        j++;
        i++;
      }
    }
  }
}
