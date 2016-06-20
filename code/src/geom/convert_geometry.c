
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


#include "convert_geometry.h"

#include "../defines.h"
#include "../util/error.h"
#include "../util/moreutil.h"

#include "../util/util1.h"

int convert_geometry(struct crust crst, double *old_v, double **new_v, int sum, int increase_resolution){
/* Convert a vector (old_v) between large/small grid in crst.
 * old_v, new_v have size [1...NP], [1...NPn].
 *
 * Input:
 *  crst: structure containing large grid (crst.nLat_out, crst.nLon_out, crst.nD_out) and small grid (crst.nLat, crst.nLon, crst.nD);
 *  old_v: original vector
 *  increase_resolution: flag. if 1, will use crst.nX_out as starting points (large cells), and calculate new values for crst.nX (small cells);
 *  		otherwise, will convert from larger to smaller cells.
 *  sum: flag. for sum==0, quantity is an intensive variable -> take average. for sum==1, takes sum (only for increase_resolution==0).
 *  		only used if (increase_resolution==0).
 *
 * Output:
 *  new_v: new vector
 *
 * NB: if old_v, new_v are the same (i.e. the refined and output grid are the same), old_v will be copied if (*new_v!=NULL),
 * in which case memory should have already been allocated for the correct size; if *new_v, will set *new_v=old_v (to save memory);
 * in this case, one should be careful not to overwrite old_v later by writing into *new_v. To be safe, allocate memory to *new_v before.
 *
 * points change along lat, then along lon, then along depth.
 * indices start from 1.
 * new_v will be initialized if NULL; otherwise, must have correct no. of elements!
 *
 */

	int P[3], Pn[3], nsub[3];
	int D1, D2, D1D2;	//no of points per dimension (old geometry);
	int D1n, D2n;		//no of points per dimension (new geometry);
	int NP,		//original no. of points
		NPn, 	//new number of points
		nsub_tot;	//ratio.
	int ind;

	NPn= (increase_resolution)? crst.N_allP : crst.nLat_out*crst.nLon_out*crst.nD_out;
	NP= (increase_resolution)? crst.nLat_out*crst.nLon_out*crst.nD_out : crst.N_allP;

	if (!crst.uniform || NP==NPn){
		if (!(*new_v)){
			*new_v=old_v;
		}
		else{
			//if memory has already been allocated for new_v, do not change it (since it may contain extra elements):
			copy_vector(old_v, new_v, NP);
		}
		return 0;
	}

	if (increase_resolution){
		D1=crst.nLat_out;	//large cells
		D2=crst.nLon_out;
		D1D2=D1*D2;
		D1n=crst.nLat;		//small cells
		D2n=crst.nLon;
	}

	else {
		D1=crst.nLat;	//small cells
		D2=crst.nLon;
		D1D2=D1*D2;
		D1n=crst.nLat_out;	//large cells
		D2n=crst.nLon_out;
	}

	//no. of grid points between high resolution geometry and low resolution geometry (per dimension);
	nsub[0]=crst.nLat/crst.nLat_out;
	nsub[1]=crst.nLon/crst.nLon_out;
	nsub[2]=crst.nD/crst.nD_out;
	nsub_tot=nsub[0]*nsub[1]*nsub[2];

	if (crst.nLat%crst.nLat_out!=0 || crst.nLon%crst.nLon_out!=0 || crst.nD%crst.nD_out!=0) {
		print_screen(" ** Error: calculation cells are not a multiple of forecast cells - can not recalculate geometry -> Using old geometry. (convert_geometry)\n");
		print_logfile(" ** Error: calculation cells are not a multiple of forecast cells - can not recalculate geometry -> Using old geometry. (convert_geometry)\n");
		if (!(*new_v)) {
			*new_v=old_v;
		}
		else {
			copy_vector(old_v, new_v, NP);
		}
		return(1);
	}

	if (!(*new_v)) *new_v=darray(1,NPn);
	for (int i=1; i<=NPn; i++) (*new_v)[i]=0.0;

	for (int i=1; i<=NP; i++){
		//reshape linear array into 3x3 array: i -> (P1,P2,P3).
		P[0]=(i-1)%D1+1;
		P[1]=((i-1)%D1D2)/D1 +1;
		P[2]=(i-1)/D1D2 +1;

		if (increase_resolution){
			for (int x=1; x<=nsub[0]; x++){
				Pn[0]=(P[0]-1)*nsub[0]+x;
				for (int y=1; y<=nsub[1]; y++){
					Pn[1]=(P[1]-1)*nsub[1]+y;
					for (int z=1; z<=nsub[2]; z++){
						Pn[2]=(P[2]-1)*nsub[2]+z;
						ind=Pn[0]+D1n*(Pn[1]-1)+D1n*D2n*(Pn[2]-1);
						(*new_v)[ind]= old_v[i];
					}
				}
			}
		}

		else {
			//calculate new indices:
			for (int n=0; n<3; n++) Pn[n]=(P[n]-1)/nsub[n]+1;

			//calculate linear index:
			ind=Pn[0]+D1n*(Pn[1]-1)+D1n*D2n*(Pn[2]-1);
			(*new_v)[ind]= (sum)? (*new_v)[ind]+old_v[i] : (*new_v)[ind]+old_v[i]*(1.0/nsub_tot);
		}
	}

	return(0);
}


int flatten_outgrid(struct crust crst, double *old_v, double **new_v, int *Nfinal){
	/* Flattens values corresponding to a full output grid to a 2D one, summing over depth layers.
	 *
	 * Input:
	 *  crst: structure containing large grid (crst.nLat_out, crst.nLon_out, crst.nD_out).
	 *  old_v: original vector
	 *
	 * Output:
	 *  new_v: new vector
	 *  Nfinal: size of *new_v.
	 *
	 * NB: if old_v, new_v are the same (e.g. if crst.uniform==0), old_v will be copied if (*new_v!=NULL),
	 * in which case memory should have already been allocated for the correct size; if *new_v, will set *new_v=old_v (to save memory);
	 * in this case, one should be careful not to overwrite old_v later by writing into *new_v. To be safe, allocate memory to *new_v before.
	 *
	 * points change along lat, then along lon, then along depth.
	 * indices start from 1.
	 * new_v will be initialized if NULL; otherwise, must have correct no. of elements!
	 *
	 */

	int N0=crst.nLat_out*crst.nLon_out*crst.nD_out;	//initial size.
	int old_i;	//old indices.

	if (!crst.uniform){
		if (!(*new_v)){
			*new_v=old_v;
		}
		else{
			//if memory has already been allocated for new_v, do not change it (since it may contain extra elements):
			copy_vector(old_v, new_v, N0);
		}
		*Nfinal=N0;
		return 0;
	}

	else{
		*Nfinal=N0/crst.nD_out;

		if (!(*new_v)){
			*new_v=darray(1,*Nfinal);
		}

		for (int i=1; i<=*Nfinal; i++){
			(*new_v)[i]=0.0;
			for (int j=1; j<=crst.nD_out; j++){
				old_i=(j-1)*(*Nfinal)+i;
				(*new_v)[i]+=old_v[old_i];
			}
		}
	}

	return 0;

}
