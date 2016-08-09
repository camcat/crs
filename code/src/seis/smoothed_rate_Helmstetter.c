
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

#include "smoothed_rate_Helmstetter.h"

double *smoothed_rate_Helmstetter(double *xgrid, double *ygrid, double dx, double dy, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord){
/* Calculates background rate from a catalog, using the algorithm from:
 * Helmstetter et al (2007) "High-resolution Time-independent Grid-based Forecast for M ≥ 5 Earthquakes in California"
 *
 * Input:
 *  xgrid, ygrid: x,y grid for which rate should be calculated. size [1...Ngrid]
 *  dx, dy: spacing between grid points.
 *  xs, ys: x,y coordinates of events. size [1...N]
 *  err= location error associated with each event. size [1...N]
 *  weights: weight assigned to each event. Can be used as flag for declustering (0/1 for excluded/selected events), or to weight some events more than others.
 *  		 if weights==NULL, all events are selected with weight=1. size [1...N]
 *  ord= 1,2: indicates if first o second nearest neighbour should be used.
 *
 * Output:
 *  vector of size [1...Ngrid] containing number of events in each grid point for the total time period.
 */

	double d, w;
	double *dist=NULL;
	int *ind, no_ind;
	double *rate, *rate_tot;

	ind=iarray(1,Ngrid);
	rate=darray(0,Ngrid);
	rate_tot=darray(1,Ngrid);
	for (int i=1; i<=Ngrid; i++) rate_tot[i]=0.0;

	switch (ord){
	//calculate first or second nearest neighbor distance for each earthquake:
		case 1:
			all_nearestneighbours(xs, ys, N, NULL, &dist);
			break;
		case 2:
			all_2ndnearestneighbours(xs, ys, N, NULL, &dist);
			break;
		default:
			print_screen("** Error:  illegal value for variable 'ord' in smoothed_rate_Helmstetter.c. \n**");
			print_logfile("** Error:  illegal value for variable 'ord' in smoothed_rate_Helmstetter.c. \n**");
			return NULL;
	}

	//todo could parallelize (omp)
	for (int eq=1; eq<=N; eq++){
		if (!weights || (weights[eq]>0.0)){	//skip if weight=0.0
			d=fmax(dist[eq],err[eq]);
			find_gridpoints_exact(ygrid, xgrid, NULL, dx, dy, 0.0, Ngrid, Ngrid, ys[eq], xs[eq], d, 0.0, 0.0, 10000, &no_ind, &ind, &rate, 1, 0);
			w= (weights)? weights[eq] : 1.0;
			for (int i=1; i<=no_ind; i++) rate_tot[ind[i]]+=w*rate[i];
		}
	}

	free_iarray(ind,1,Ngrid);
	free_darray(rate,0,Ngrid);

	return rate_tot;
}

double *smoothed_rate_Helmstetter_nonuni(double *xgrid, double *ygrid, int Ngrid, double *xs, double *ys, double *err, double *weights, int N, int ord){
	/* Very similar to Helmstetter function, but for non uniform grid.
	 * Since cannot use exact integral of Gaussian function in each grid point, will only use value at the grid cell center.
	 *
	 * Input:
	 *  xgrid, ygrid: x,y grid for which rate should be calculated. size [1...Ngrid]
	 *  dx, dy: spacing between grid points.
	 *  xs, ys: x,y coordinates of events. size [1...N]
	 *  err= location error associated with each event. size [1...N]
	 *  weights: weight assigned to each event. Can be used as flag for declustering (0/1 for excluded/selected events), or to weight some events more than others.
	 *  		 if weights==NULL, all events are selected with weight=1. size [1...N]
	 *  ord= 1,2: indicates if first o second nearest neighbour should be used.
	 *
	 * Output:
	 *  vector of size [1...Ngrid] containing number of events in each grid point for the total time period.
	 */

	double d, w;
	double *dist=NULL;
	int *ind, no_ind;
	double *rate, *rate_tot;

	ind=iarray(1,Ngrid);
	rate=darray(1,Ngrid);
	rate_tot=darray(1,Ngrid);
	for (int i=1; i<=Ngrid; i++) rate_tot[i]=0.0;

	switch (ord){
		case 1:
			all_nearestneighbours(xs, ys, N, NULL, &dist);
			break;
		case 2:
			all_2ndnearestneighbours(xs, ys, N, NULL, &dist);
			break;
		default:
			print_screen("** Error:  illegal value for variable 'ord' in smoothed_rate_Helmstetter_nonuni. \n**");
			print_logfile("** Error:  illegal value for variable 'ord' in smoothed_rate_Helmstetter_nonuni. \n**");
			return NULL;
	}

	//todo could parallelize (omp)
	for (int eq=1; eq<=N; eq++){
		if (!weights || (weights[eq]>0.0)){
			d=fmax(dist[eq],err[eq]);
			// smooth event based on coordinates of the center of each grid point (not exact calculation as in Helmstetter function).
			find_gridpoints(ygrid, xgrid, NULL, NULL, Ngrid, ys[eq], xs[eq], d, 0.0, 0.0001, 1, &no_ind, &ind, &rate, 1, 0);
			w= (weights)? weights[eq] : 1.0;
			for (int i=1; i<=no_ind; i++) rate_tot[ind[i]]+=w*rate[i];
		}
	}

	return rate_tot;
}

double *smoothed_rate_Helmstetter_cat(struct catalog cat, struct crust crst, double *weights, int ord) {
	/* Wrapper for applying function Helmstetter to variables in catalog structure.
	 * Returns background rate on a 2D grid (only uses horizontal distance).
	 *
	 * Input:
	 *  cat: structure containing earthquakes to be used to estimate background rate.
	 *  crst: contains grid geometry.
	 *  weights: each event in cat will be scaled by a factor weights (see decluster_catalog function for details).
	 *  ord: order (e.g. if ord=2, 2nd nearest neighbour distance will be used for smoothing).
	 *
	 * Output: background rate in a vector the same size as crst.nLon*crst.nLat (i.e. horizontal plane).
	 */

	double dx, dy;

	dx=(DEG2RAD*crst.dlon)*Re*cos(DEG2RAD*crst.lat0);
	dy=(DEG2RAD*crst.dlat)*Re;

	return smoothed_rate_Helmstetter(crst.east, crst.north, dx, dy, crst.nLon*crst.nLat, cat.east0, cat.north0, cat.err, weights, cat.Z, ord);

}

double *fit_depth(double *zgrid, double dz, int Ngrid, double *zs, double *err, double *weights, int N){
/* Returns the number of events in depth bins.
 *
 * Input:
 *
 *  zgrid, ygrid: vertical grid for which rate should be calculated. size [1...Ngrid]
 *  dz: spacing between grid points.
 *  zs: depth of events. size [1...N]
 *  err= location error associated with each event. size [1...N]. Used as st. dev. for smoothing.
 *  weights: weight assigned to each event. Can be used as flag for declustering (0/1 for excluded/selected events), or to weight some events more than others.
 *  		 if weights==NULL, all events are selected with weight=1. size [1...N]
 *
 * Output:
 *  vector of size [1...Ngrid] containing number of events in each grid point for the total time period.
 *
 */

	double *prob0, *prob;
	double probCum;
	double rz, w;

	prob0=darray(1,Ngrid);
	prob=darray(1,Ngrid);

	for (int n=1; n<=Ngrid; n++) prob[n]=0;

	//todo could parallelize(omp)
	for (int eq=1; eq<=N; eq++){
		probCum=0;
		for (int n=1; n<=Ngrid; n++){
			rz=fabs(zgrid[n]-zs[eq]);
			//smooth event across grid points:
			prob0[n]=exact_prob_1d(rz, dz, err[eq]);
		}
		w=(weights)? weights[eq] : 1.0;
		for (int n=1; n<=Ngrid; n++) prob[n]+=w*prob0[n];
	}

	free_darray(prob0,1,Ngrid);

	return prob;
}


