
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

#include "background_rate.h"

int background_rate(char *catfile, struct crust *crst_in, struct tm reftime,
		double Mmain, double *minmag, double *rate, double **rate_grid, double dR, double dZ, double min_smoothing, int ord) {
	/*
	 * Calculates background based on a catalog. It will decluster it using the Knopoff-Gardner 74 criterium and calculate the smoothed rate using the
	 * algorithm from Helmstetter et al (2007).
	 *
	 * Input:
	 *  catfile: file containing catalog in ZMAP format.
	 *  crst_in: crust structure containing calculation geometry.
	 *  reftime: reference time (IssueTime)
	 *  Mmain: minimum magnitude to be used for declustering
	 *  dR, dZ: extra horizontal, vertical distance to be used when selecting earthquakes from catalog.
	 *  ord: order for Helmstetter algorithm (e.g. if ord=2, 2nd nearest neighbour distance will be used as st. dev. for smoothing).
	 *
	 *
	 * Output:
	 *  minmag: completeness magnitude (seismicity rate if for events Mw>=minmag)
	 *  *rate: scalar value of seismicity rate (no. of events per day)
	 *  **rate_grid: array containing the spatial distribution of seismicity (by convention adds up to crst.N_allP). size [1...crst.N_allP]
	 */

	struct catalog cat;
	struct crust crst=*crst_in;
	double *zlist=NULL;
	double t0;
	int NP= crst.uniform? (crst.nLat*crst.nLon) : crst.N_allP;
	double T, *rate_h=NULL, *rate_v=NULL;
	int h_ind, v_ind, *sel, sel_no=0, err;
	double *weights=NULL;
	cat.Mc=20.0;	//by convention, this value will make readZMAP find completeness magnitude.

	read_firstlineZMAP(catfile, reftime, &t0);	//just to know starting time
	err=readZMAP(&cat, (struct eqkfm **) 0, NULL, catfile, crst, reftime, t0, 0.0, t0, 0.0, 10, 0.0, dR, dZ, 0.0, 0);

	if(err || cat.Z==0) {
		return 1;
	}

	// Decluster catalog:
	// sel contains flag (1,0) for selected/excluded events, weights is a weight accounting for the fact that declustering algorithm also removes some background events.
	// (see comments in decluster_catalog).
	sel=decluster_catalog(cat, Mmain, &weights, 0);
	T=cat.tend-cat.tstart;

	for (int i=1; i<=cat.Z; i++) {
		cat.err[i]=fmax(cat.err[i], min_smoothing);	//values used later (in smoothed_rate_Helmstetter.c).
		sel_no+=sel[i];
	}

	// If events are too few, background rate estimation is less reliable:
	if (sel_no<=1000){
		if (sel_no<=100 || T<1.0)	{
			print_screen("**Warning: too few events, or too short time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
			print_logfile("**Warning: too few events, or too short a valid time period found in catalog [%ld events, %.3lf days]: can not estimate background rate.** \n", cat.Z, T);
			return 1;
		}
		else{
			print_screen("**Warning: only %ld events used to calculate background seismicity rate.** \n", cat.Z);
			print_logfile("**Warning: only %ld events used to calculate background seismicity rate.** \n", cat.Z);
		}
	}
	else {
		print_logfile("**Background seismicity rate calculated from catalog %s (%d events between t=[%.5lf, %.5lf]).** \n", catfile, sel_no, cat.tstart, cat.tend);
	}

	// find minmag and background rate;
	if (minmag) *minmag= cat.Mc;
	if (rate) {
		*rate=0;
		for (int i=1; i<=cat.Z; i++) {
			*rate+=weights[i];	//instead of adding 1, use weights[i] which accounts for removal of background events by declustering algorithm.
		}
		(*rate)*=1.0/T;
	}


	if (crst.uniform){
		zlist=darray(1,crst.nD);
		for (int i=1; i<=crst.nD; i++) zlist[i]=crst.depth[1+(i-1)*NP];

		// Depth-average rate on horizontal grid calculated from Helmstetter algorithm;
		// depth distribution obtained separately.
		rate_h=smoothed_rate_Helmstetter_cat(cat, crst, weights, ord);
		rate_v=fit_depth(zlist, zlist[2]-zlist[1], crst.nD, cat.depths0, cat.verr, weights, cat.Z);

		//normalize rate vectors:
		normv(rate_h, NP);
		normv(rate_v, crst.nD);

		*rate_grid=darray(1,crst.N_allP);
		for (int p=1; p<=crst.N_allP; p++) {
			//calculate linear index:
			h_ind=(p-1)%NP+1;
			v_ind=(p-1)/NP+1;
			// by convention, (*rate_grid)[p] add up to crst.N_allP (since for simplicity they are all assumed to be 1 if *rate_grid==NULL)
			(*rate_grid)[p]=rate_h[h_ind]*rate_v[v_ind]*crst.N_allP;
		}
	}

	else{
		// Assume vertically homogeneous rate, since can not use depth info (non uniform grid does not necessarily have layers).
		*rate_grid=smoothed_rate_Helmstetter_nonuni(crst.east, crst.north, crst.N_allP, cat.east0, cat.north0, cat.err, weights, cat.Z, ord);
		normv(*rate_grid, crst.N_allP);
	}

	free_cat(cat);
	if (zlist) free_darray(zlist,1,crst.nD);
	if (rate_h)	free_darray(rate_h, 1, NP);
	if (rate_v) free_darray(rate_v, 1, crst.nD);

	return 0;

}

