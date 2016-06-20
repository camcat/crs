
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


#include <math.h>

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "fit_splines.h"

#include "../util/util1.h"

void splines_eqkfm(struct eqkfm **eqkfm_aft, int Nas, int NF, double *times1, double *times2, int L, long *seed){
	double **slipbefore_st, **slipbefore_di, **slip_after_st, **slip_after_di, **slipbefore_op, **slip_after_op;
	int *NP, NP_tot;
	int is_str, is_dip, is_open;

	NP=iarray(0,NF-1);
	NP_tot=0;

	for (int f=0; f<NF; f++) {
		NP[f]= (*eqkfm_aft)[f].np_di*(*eqkfm_aft)[f].np_st;
		NP_tot=MAX(NP_tot,NP[f]);
	}

	slipbefore_st=d2array(1,NP_tot,1,Nas);
	slip_after_st=d2array(1,NP_tot,1,L);
	slipbefore_di=d2array(1,NP_tot,1,Nas);
	slip_after_di=d2array(1,NP_tot,1,L);
	slipbefore_op=d2array(1,NP_tot,1,Nas);
	slip_after_op=d2array(1,NP_tot,1,L);
	if (!slipbefore_st | !slip_after_st | !slipbefore_di | !slip_after_di | !slipbefore_op | !slip_after_op) memory_error_quit;


	for (int f=0; f<NF; f++){

		//check which arrays exist:
		is_str= ((*eqkfm_aft)[f].allslip_str==NULL) ? 0 : 1 ;
		is_dip= ((*eqkfm_aft)[f].allslip_dip==NULL) ? 0 : 1 ;
		is_open= ((*eqkfm_aft)[f].allslip_open==NULL) ? 0 : 1 ;

		for (int pt=1; pt<=NP[f]; pt++){
			for (int t=0; t<Nas; t++){
				if (is_str) slipbefore_st[pt][t+1]=(*eqkfm_aft)[f].allslip_str[t][pt];
				if (is_dip) slipbefore_di[pt][t+1]=(*eqkfm_aft)[f].allslip_dip[t][pt];
				if (is_open) slipbefore_op[pt][t+1]=(*eqkfm_aft)[f].allslip_open[t][pt];
			}
		}
		if (is_str) fit_splines(times1, times2, Nas, L, NP[f], slipbefore_st, (double *) 0, &slip_after_st, 2, seed);
		if (is_dip) fit_splines(times1, times2, Nas, L, NP[f], slipbefore_di, (double *) 0, &slip_after_di, 2, seed);
		if (is_open) fit_splines(times1, times2, Nas, L, NP[f], slipbefore_op, (double *) 0, &slip_after_op, 2, seed);

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
				if (is_str) (*eqkfm_aft)[f].allslip_str[t][pt]=slip_after_st[pt][t+1];
				if (is_dip) (*eqkfm_aft)[f].allslip_dip[t][pt]=slip_after_di[pt][t+1];
				if (is_open) (*eqkfm_aft)[f].allslip_open[t][pt]=slip_after_op[pt][t+1];
			}
		}
	}

	free_iarray(NP,0,NF-1);

	free_d2array(slipbefore_st,1,NP_tot,1,Nas);
	free_d2array(slip_after_st,1,NP_tot,1,L);
	free_d2array(slipbefore_di,1,NP_tot,1,Nas);
	free_d2array(slip_after_di,1,NP_tot,1,L);
	free_d2array(slipbefore_op,1,NP_tot,1,Nas);
	free_d2array(slip_after_op,1,NP_tot,1,L);

}

