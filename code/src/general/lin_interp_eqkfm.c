
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

#include "../defines.h"
#include "lin_interp_eqkfm.h"
#include "../util/util1.h"

void lin_interp_eqkfm(struct eqkfm **eqkfm_aft, int NF, int L, double *times2, int *ind){

/* Rewrites stressing history in a eqkfm_aft structure based on time new steps given in times2.
 *
 * Input:
 *
 *  eqkfm is the array of eqkfm elements refering to a single event (each element corresponts to a subfault).
 *  NF is the number of subfaults. eqfm[0...NF-1].
 *  NT= no of time steps
 *  times= time steps [0...NT-1]
 *  ind= lists indices of elements in times which where originally in eqkfm[x].ts (NB the 0th element is eqfm[x].t).
 *
 * Output: changes values of eqkfm[x].allslip_xxx arrays, and eqkfm[x].ts, now referring to new time steps.
 *
 */

	double **slipbefore_st, **slipbefore_di, **slip_after_st, **slip_after_di, **slipbefore_op, **slip_after_op;
	double *times1;
	int *NP, NP_tot, Nas;
	int is_str, is_dip, is_open;

	Nas=(*eqkfm_aft)[0].nosnap;
	NP=iarray(0,NF-1);
	NP_tot=0;

	for (int f=0; f<NF; f++) {
		NP[f]= (*eqkfm_aft)[f].np_di*(*eqkfm_aft)[f].np_st;
		NP_tot=MAX(NP_tot,NP[f]);
	}

	//TODO should not allocate memory which may not be needed.
	slipbefore_st=d2array(1,NP_tot,0,Nas);
	slip_after_st=d2array(1,NP_tot,0,L);
	slipbefore_di=d2array(1,NP_tot,0,Nas);
	slip_after_di=d2array(1,NP_tot,0,L);
	slipbefore_op=d2array(1,NP_tot,0,Nas);
	slip_after_op=d2array(1,NP_tot,0,L);
	
	if(!slipbefore_st | !slip_after_st | !slip_after_st | !slip_after_di | !slip_after_op | !slip_after_op)	memory_error_quit;

	times1=darray(0,Nas);
	times1[0]=(*eqkfm_aft)[0].t;	//start time when stress=0.
	for (int t=0; t<Nas; t++) times1[t+1]=(*eqkfm_aft)[0].ts[t];

	for (int f=0; f<NF; f++){
		//check which arrays exist:
		is_str= ((*eqkfm_aft)[f].allslip_str==NULL) ? 0 : 1 ;
		is_dip= ((*eqkfm_aft)[f].allslip_dip==NULL) ? 0 : 1 ;
		is_open= ((*eqkfm_aft)[f].allslip_open==NULL) ? 0 : 1 ;

		for (int pt=1; pt<=NP[f]; pt++){
			if (is_str) slipbefore_st[pt][0]=0.0;
			if (is_dip) slipbefore_di[pt][0]=0.0;
			if (is_open) slipbefore_op[pt][0]=0.0;
			for (int t=0; t<Nas; t++){
				if (is_str) slipbefore_st[pt][t+1]=(*eqkfm_aft)[f].allslip_str[t][pt];
				if (is_dip) slipbefore_di[pt][t+1]=(*eqkfm_aft)[f].allslip_dip[t][pt];
				if (is_open) slipbefore_op[pt][t+1]=(*eqkfm_aft)[f].allslip_open[t][pt];
			}
		}

		if (is_str) fit_lin(times1, times2, Nas+1, L, ind,  NP[f], slipbefore_st, slip_after_st);
		if (is_dip) fit_lin(times1, times2, Nas+1, L, ind,  NP[f], slipbefore_di, slip_after_di);
		if (is_open) fit_lin(times1, times2, Nas+1, L, ind,  NP[f], slipbefore_op, slip_after_op);

		(*eqkfm_aft)[f].ts=times2;
		(*eqkfm_aft)[f].nosnap=L;
		//reallocate memory:
		if (is_str){
			free_d2array((*eqkfm_aft)[f].allslip_str,0,Nas-1, 1,NP[f]);
			(*eqkfm_aft)[f].allslip_str=d2array(0,L-1,1,NP[f]);
		}
		if (is_dip){
			free_d2array((*eqkfm_aft)[f].allslip_dip,0,Nas-1, 1,NP[f]);
			(*eqkfm_aft)[f].allslip_dip=d2array(0,L-1,1,NP[f]);
		}
		if (is_open){
			free_d2array((*eqkfm_aft)[f].allslip_open,0,Nas-1, 1,NP[f]);
			(*eqkfm_aft)[f].allslip_open=d2array(0,L-1,1,NP[f]);
		}

		free_darray((*eqkfm_aft)[f].tot_slip,0,Nas-1);
		(*eqkfm_aft)[f].tot_slip=darray(0,L-1);

		for (int t=0; t<L; t++){
			for (int pt=1; pt<=NP[f]; pt++){
				if (is_str) (*eqkfm_aft)[f].allslip_str[t][pt]=slip_after_st[pt][t];
				if (is_dip) (*eqkfm_aft)[f].allslip_dip[t][pt]=slip_after_di[pt][t];
				if (is_open) (*eqkfm_aft)[f].allslip_open[t][pt]=slip_after_op[pt][t];
			}
		}
	}

	free_iarray(NP,0,NF-1);

	free_d2array(slipbefore_st,1,NP_tot,0,Nas);
	free_d2array(slip_after_st,1,NP_tot,0,L);
	free_d2array(slipbefore_di,1,NP_tot,0,Nas);
	free_d2array(slip_after_di,1,NP_tot,0,L);
	free_d2array(slipbefore_op,1,NP_tot,0,Nas);
	free_d2array(slip_after_op,1,NP_tot,0,L);
	free_darray(times1,0,Nas);

}


void fit_lin(double *times1, double *times2, int Nas, int L, int *ind, int NP, double **slipbefore_st, double **slip_after){

/* Linear interpolation.
 *
 *Input:
 * times1[0...Nas-1]= initial time steps
 * slipbefore_st[1...NP][0...Nas-1]=initial values. If not given, will assume it's 1 for all patches.
 * times2[0...L-1]= final time steps
 * ind= indices of times2 elements which are also in times1.
 *
 * (*slip_after)[1...NP][0...L-1]=final values.
 */


	double delta_slip, delta_t, dt;
	int l;


	for (int p=1; p<=NP; p++){

	  //steps before the first snapshot:
	  l=0;
	  while (l<ind[0]) {
		  slip_after[p][l]=0.0;
		  l++;
	  }

	  //loop over original time steps:
	  for (int j=0; j<Nas-1; j++){
		  //original stress change:
		  delta_slip= (slipbefore_st)? slipbefore_st[p][j+1]-slipbefore_st[p][j] : 1.0;
		  delta_t= times1[j+1]-times1[j];

		  //new stress changes:
		  for (int l=ind[j]; l<ind[j+1]; l++){
			  dt=times2[l+1]-times2[l];
			  slip_after[p][l]=delta_slip*dt/delta_t;
		  }
	  }
	  //steps after the last snapshot:
	  for (int l=ind[Nas-1]; l<L; l++) slip_after[p][l]=0.0;	//stationary after last time step.

	}

}

