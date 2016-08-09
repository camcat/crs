
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


#include "mem_mgmt.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "../defines.h"
#include "../util/moreutil.h"

#include "../util/util1.h"

#define offset 1

void reduce_eqkfm_memory(struct eqkfm *eqkfm0, int NF){
/* Frees empty arrays in eqkfm0 structure.
 * Input:
 *  eqkfm0= pointer to eqkfm array. (*eqkfm)[0...NF-1]
 *  NF=no of faults, i.e. no. of *eqkfm elements.
 * Output:
 *  frees eqkfm.slip_str, eqkfm.slip_dip, eqkfm.slip_open arrays when needed.
 */

	double 	toll=1e-10;	//tolerance
	int is_str, is_dip, is_open;	//flags used to determine which components of deformations are needed (to save memory).

	//check if all elements are 0, and is so set flag.
	for (int nf=0; nf<NF; nf++){
		check_empty_eqkfm(eqkfm0[nf], toll, &is_str, &is_dip, &is_open);

		//free memory if elements are all 0.
		if (is_str==0){
			//if condition needed since the element may have been freed before:
			if (eqkfm0[nf].slip_str) free_darray(eqkfm0[nf].slip_str,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			eqkfm0[nf].slip_str=NULL;
		}
		if (is_dip==0){
			if (eqkfm0[nf].slip_dip) free_darray(eqkfm0[nf].slip_dip,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			eqkfm0[nf].slip_dip=NULL;
		}
		if (is_open==0){
			if (eqkfm0[nf].open) free_darray(eqkfm0[nf].open,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			eqkfm0[nf].open=NULL;
		}
	}
}

void check_empty_eqkfm(struct eqkfm eqkfm0, double toll, int *is_str, int *is_dip, int *is_open){
/* Checks if eqkfm0.slip_str(dip,open) are empty (i.e. all 0s within tolerance) and returns boolean value.
 *
 * Input:
 *  eqkfm0: single structure eqkfm.
 *  toll: tolerance
 */

	*is_str=0;
	*is_dip=0;
	*is_open=0;

	for (int p=1; p<=eqkfm0.np_st*eqkfm0.np_di; p++){
		if (eqkfm0.slip_str && (*is_str || fabs(eqkfm0.slip_str[p])>toll)) *is_str=1;
		if (eqkfm0.slip_dip && (*is_dip || fabs(eqkfm0.slip_dip[p])>toll)) *is_dip=1;
		if (eqkfm0.open && (*is_open || fabs(eqkfm0.open[p])>toll)) *is_open=1;

		if (*is_str && *is_dip && *is_open) break;
	}
}


void init_crst(struct crust *crst){
/* Initialize variables in crust structure to default values. */


	(*crst).S=NULL;
	(*crst).list_allP=NULL;
	(*crst).lat=NULL;
	(*crst).lon=NULL;
	(*crst).depth=NULL;
	(*crst).lat_out=NULL;
	(*crst).lon_out=NULL;
	(*crst).depth_out=NULL;
	(*crst).dAgrid=NULL;
	(*crst).str0=darray(0,0);
	(*crst).dip0=darray(0,0);
	(*crst).rake0=darray(0,0);
	(*crst).east=NULL;
	(*crst).north=NULL;
	(*crst).rate0=NULL;
	(*crst).mags=NULL;
	(*crst).GRmags=NULL;
	(*crst).nofmzones=1;
	(*crst).fmzone=NULL;
	(*crst).variable_fixmec=0;

	return;
}

void init_cat1(struct catalog *cat, int Zsel){
/* Initialize variables in catalog structure to default values. */

	(*cat).Z=Zsel;
	(*cat).t = darray(1, Zsel);
	(*cat).mag = darray(1, Zsel);
	(*cat).lat0 = darray(1, Zsel);
	(*cat).lon0 = darray(1, Zsel);
	(*cat).east0 = darray(1, Zsel);
	(*cat).north0 = darray(1, Zsel);
	(*cat).depths0 = darray(1, Zsel);
	(*cat).err = darray(1, Zsel);
	(*cat).verr = darray(1, Zsel);
	(*cat).ngrid = iarray(1, Zsel);
	//just allocate first level since subarrays may have different length (and will be initialized later).
	(*cat).ngridpoints=i2array_firstlevel(Zsel);
	(*cat).weights=d2array_firstlevel(Zsel);
	(*cat).b=1.0;

}

struct eqkfm *eqkfm_array(long n1, long n2){
/* Allocate memory to array of eqkfm. */
	struct eqkfm *v;
	v= (struct eqkfm *) malloc((size_t) ((n2-n1+1+offset)*sizeof(struct eqkfm)));

	for (int i=offset; i<=n2-n1+offset; i++){
		v[i].slip_str= NULL;
		v[i].slip_dip= NULL;
		v[i].open= NULL;
		v[i].ts=NULL;
		v[i].nosnap=0;
		v[i].allslip_str= NULL;
		v[i].allslip_dip= NULL;
		v[i].allslip_open= NULL;
		v[i].tevol=NULL;
		v[i].pos_s= NULL;
		v[i].pos_d= NULL;
		v[i].selpoints= NULL;
		v[i].distance= NULL;
		v[i].is_slipmodel=0;
		v[i].np_st=v[i].np_di=v[i].whichfm=v[i].nsel=0;
		v[i].t=0;
		v[i].lat=0;
		v[i].lon=0;
		v[i].depth=0;
		v[i].mag=0;
		v[i].tot_slip=darray(0,0);	//only need one element for earthquake sources, will reallocate for afterslip.
		v[i].L=0;
		v[i].W=0;
		v[i].str1=0;
		v[i].str2=0;
		v[i].dip1=0;
		v[i].dip2=0;
		v[i].rake1=0;
		v[i].rake2=0;
		v[i].index_cat=0;
		v[i].cuts_surf=0;
		v[i].co_aft_pointer=NULL;

	}
	return v-n1+offset;
}

struct pscmp *pscmp_array(long n1, long n2){
/* Allocate memory to array of pscmp. */
	struct pscmp *v;
	v= (struct pscmp *) malloc((size_t) ((n2-n1+1+offset)*sizeof(struct pscmp)));

	for (int i=offset; i<=n2-n1+offset; i++){
		v[i].t=0.0;
		v[i].fdist=NULL;	//distance to fault
		v[i].S=NULL;
		v[i].cmb=NULL;
		v[i].Z=0;
		v[i].st1=NULL;
		v[i].di1=NULL;
		v[i].ra1=NULL;
		v[i].st2=NULL;
		v[i].di2=NULL;
		v[i].ra2=NULL;
		v[i].which_pts=NULL;
		v[i].nsel=0;
	    v[i].index_cat=0;
	    v[i].NF=0;
	}
	return v-n1+offset;
}

struct pscmp *pscmp_arrayinit(struct crust v0, long n1, long n2){
/* Allocate memory to array of pscmp, and initialize variables copied from v0. */

		struct pscmp *v;
		v= (struct pscmp *) malloc((size_t) ((n2-n1+1+offset)*sizeof(struct pscmp)));

		for (int i=offset; i<=n2-n1+offset; i++){
			v[i].t=0.0;
			v[i].fdist=NULL;	//distance to fault
			v[i].S=NULL;
			v[i].cmb=NULL;
			v[i].Z=0;
			v[i].st1=NULL;
			v[i].di1=NULL;
			v[i].ra1=NULL;
			v[i].st2=NULL;
			v[i].di2=NULL;
			v[i].ra2=NULL;
			v[i].which_pts=NULL;
			v[i].nsel=0;
		    v[i].index_cat=0;
		    v[i].NF=0;
		}
		return v-n1+offset;;
}



void free_cat(struct catalog cat){
/* Deallocates memory from variables in catalog structure.
 */

	//assumes that elements have been initialized at position 1.
	free_darray(cat.t,1, 0);
	free_darray(cat.mag,1, 0);
	free_darray(cat.lat0,1, 0);
	free_darray(cat.lon0,1, 0);
	free_darray(cat.east0,1, 0);
	free_darray(cat.north0,1, 0);
	free_darray(cat.depths0,1, 0);
	free_iarray(cat.ngrid,1, 0);
	free_i2array_firstlevel(cat.ngridpoints,1,cat.Z,1,0); //uses 0 for upper index, since it doesn't matter.
	free_d2array_firstlevel(cat.weights,1,cat.Z, 1,0);
}

