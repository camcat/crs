
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

#include "../util/util1.h"
#include "eqkfm_copy.h"

void empty_eqkfm(struct eqkfm *eqkfm0){
/*
 * Sets eqkfm0 to an empty slip model.
 */

	(*eqkfm0).np_st=(*eqkfm0).np_di=0;
	(*eqkfm0).nsel=0;
	(*eqkfm0).nosnap=0;
	(*eqkfm0).tot_slip=darray(0,0);	//only need one element for earthquake slip models, will reallocate for afterslip.
	(*eqkfm0).L=0;
	(*eqkfm0).W=0;
	(*eqkfm0).ts=NULL;
	(*eqkfm0).tevol=NULL;
	(*eqkfm0).slip_str=NULL;
	(*eqkfm0).slip_dip=NULL;
	(*eqkfm0).open=NULL;
	(*eqkfm0).allslip_str=NULL;
	(*eqkfm0).allslip_dip=NULL;
	(*eqkfm0).allslip_open=NULL;
	(*eqkfm0).pos_s=NULL;
	(*eqkfm0).pos_d=NULL;
	(*eqkfm0).distance=NULL;
	(*eqkfm0).selpoints=NULL;		//indices of cell points affected by this event.
}

void copy_eqkfm_nolocation_noindex_notime(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
/*
 * Copies eqkfm1 into eqkfm2, but preserves time, location and index_cat.
 */


	double lat, lon, depth, t;
	int i;

	i=(*eqkfm2).index_cat;
	lat=(*eqkfm2).lat;
	lon=(*eqkfm2).lon;
	depth=(*eqkfm2).depth;
	t=(*eqkfm2).t;

	copy_eqkfm_all(eqkfm1, eqkfm2);

	(*eqkfm2).t=t;
	(*eqkfm2).index_cat=i;
	(*eqkfm2).lat=lat;
	(*eqkfm2).lon=lon;
	(*eqkfm2).depth=depth;

}

void copy_eqkfm_noindex_notime(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
/*
 * Copies eqkfm1 into eqkfm2, but preserves time and index_cat.
 */
	double lat, lon, depth, t;
	int i;

	i=(*eqkfm2).index_cat;
	t=(*eqkfm2).t;

	copy_eqkfm_all(eqkfm1, eqkfm2);

	(*eqkfm2).t=t;
	(*eqkfm2).index_cat=i;
}

void copy_eqkfm_attributes(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
/*
 * Copies various properties of eqkfm1 into eqkfm2 (not focal mechanisms and slip models).
 */

	//general properties:
	(*eqkfm2).is_slipmodel=eqkfm1.is_slipmodel;
    (*eqkfm2).index_cat=eqkfm1.index_cat;	//index of event in catalog (only for the catalog used for LL calculation).
	(*eqkfm2).nsel=eqkfm1.nsel;
	(*eqkfm2).selpoints= eqkfm1.selpoints;
	(*eqkfm2).distance= eqkfm1.distance;
	(*eqkfm2).cuts_surf=eqkfm1.cuts_surf;
	(*eqkfm2).top=eqkfm1.top;
	(*eqkfm2).nosnap=eqkfm1.nosnap;


	//earthquake/afterslip properties:
	(*eqkfm2).t=eqkfm1.t;
	(*eqkfm2).ts=eqkfm1.ts;
	(*eqkfm2).tevol=eqkfm1.tevol;
	(*eqkfm2).mag=eqkfm1.mag;
	(*eqkfm2).lat=eqkfm1.lat;
	(*eqkfm2).lon=eqkfm1.lon;
	(*eqkfm2).east=eqkfm1.east;
	(*eqkfm2).north=eqkfm1.north;
	(*eqkfm2).depth=eqkfm1.depth;

	return;

}

void copy_eqkfm_focmec(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
/*
 * Copies focal mechanisms of eqkfm1 into eqkfm2.
 */

	//focal mechanism:
	(*eqkfm2).whichfm=eqkfm1.whichfm;
	(*eqkfm2).str1=eqkfm1.str1;
	(*eqkfm2).dip1=eqkfm1.dip1;
	(*eqkfm2).rake1=eqkfm1.rake1;
	(*eqkfm2).str2=eqkfm1.str2;
	(*eqkfm2).dip2=eqkfm1.dip2;
	(*eqkfm2).rake2=eqkfm1.rake2;

	return;

}

void copy_eqkfm_slipmodel(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
/*
 * Copies slip model of eqkfm1 into eqkfm2.
 */
	//slipmodel properties:
	(*eqkfm2).tot_slip=eqkfm1.tot_slip;
	(*eqkfm2).L=eqkfm1.L;
	(*eqkfm2).W=eqkfm1.W;
	(*eqkfm2).np_st=eqkfm1.np_st;
	(*eqkfm2).np_di=eqkfm1.np_di;
	(*eqkfm2).pos_s=eqkfm1.pos_s;
	(*eqkfm2).pos_d=eqkfm1.pos_d;
	(*eqkfm2).slip_str=eqkfm1.slip_str;
	(*eqkfm2).slip_dip=eqkfm1.slip_dip;
	(*eqkfm2).open=eqkfm1.open;
	(*eqkfm2).allslip_str=eqkfm1.allslip_str;
	(*eqkfm2).allslip_dip=eqkfm1.allslip_dip;
	(*eqkfm2).allslip_open=eqkfm1.allslip_open;

	return;

}

void copy_eqkfm_all(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
/*
 * Copies all properties of eqkfm1 into eqkfm2.
 */

	copy_eqkfm_attributes(eqkfm1, eqkfm2);
	copy_eqkfm_focmec(eqkfm1, eqkfm2);
	copy_eqkfm_slipmodel(eqkfm1, eqkfm2);

	return;
}


void copy_eqkfm_noslipmodel(struct eqkfm eqkfm1, struct eqkfm *eqkfm2){
/*
 * Copies eqkfm1 into eqkfm2, preserving slip model.
 */

	copy_eqkfm_attributes(eqkfm1, eqkfm2);
	copy_eqkfm_focmec(eqkfm1, eqkfm2);

	return;

}
