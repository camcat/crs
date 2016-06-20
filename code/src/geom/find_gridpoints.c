
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


#include "find_gridpoints.h"
#include <math.h>
#include <stddef.h>

#include "../defines.h"
#include "../util/error.h"
#include "../util/moreutil.h"

#include "../util/util1.h"

int find_closest_point(double *ys, double *xs, double *depths, int N, double y, double x, double Depth){

	/* Selects the closest point associated with a given source.
	 *
	 * Input:
	 *  ys, xs, depths: arrays with grid points coordinates. Range [1...N]
	 *  y,x: earthquake coordinates
	 *  Depth: earthquake depth
	 *
	 * Returns:
	 *  index of selected points
	 *
	 */

	double r, rmin=1e30, dz, probCum, *prob;
	int p,p2;
	int closestp;
	int ngridj;
	double y1, y2, x1, x2, D1, D2, A, Vfrac=1;

	for (p=1; p<=N; p++){
		r= sqrt(pow(ys[p]-y,2)+pow(xs[p]-x,2)+pow(depths[p]-Depth,2));
		if (r<rmin){
			rmin=r;
			closestp=p;
		}
	}
	return closestp;

}

	

int find_gridpoints(double *ys, double *xs, double *dAs, double *depths, int N, double y, double x, double SD, double Depth, double SDd,
		int cut_sd, int *ngridj0, int **ngridpointj, double **weightsj, int inside, int d3){

	/* Selects the grid points associated with a given source.
	 * Events are weighted according to a radial Gaussian distribution with a cutoff at SD*cut_sd (SDd*cut_sd vertically).
	 * If the source events is close to a boundary, the sum of the weights is <1,since part of the area affected by the earthquake is outside the domain.
	 * In this case the sum of the weights equal the ratio between the total grid point area (sum(dAs)) and the area comprised in the cutoff radius.
	 *
	 * Input:
	 *  ys, xs, depths: arrays with grid points coordinates. Range [1...N]
	 *  inside: flag indicating whether the sum of the weights should be set to 1.
	 *  d3: flag indicating if the 3D distance should be used.
	 *  y,x: earthquake coordinates
	 *  SD: standard deviation (horizontal)
	 *  Depth: earthquake depth
	 *  SDd: standard deviation (depth)
	 *
	 * Output:
	 *  ngridj0: no. of selected points
	 *  ngridpointj: list of indices of selected points
	 *  weightsj: weights of each point
	 *
	 */

	double r, rmin=1e30, dz, probCum, *prob;
	int *ngridpointj_temp;
	int p,p2;
	int closestp;
	int ngridj;
	double y1, y2, x1, x2, D1, D2, A, Vfrac=1;

	prob=darray(0,N+1);
	ngridpointj_temp=iarray(0,N+1);	//temporary array to store indices of selected points.
	if (!dAs) inside=1; 	//can't calculate total area, so assume it's all inside.


	ngridj=0;
	probCum=0;

	y1=y-cut_sd*SD;
	y2=y+cut_sd*SD;
	x1=x-cut_sd*SD;
	x2=x+cut_sd*SD;
	D1=Depth-cut_sd*SDd;
	D2=Depth+cut_sd*SDd;

	A=0;	//total area occupied by selected grid points

	for (p=1; p<=N; p++){
		r= (d3)? sqrt(pow(ys[p]-y,2)+pow(xs[p]-x,2)+pow(depths[p]-Depth,2)) : sqrt(pow(ys[p]-y,2)+pow(xs[p]-x,2));
		if (r<rmin){
			rmin=r;
			closestp=p;
		}
		if (ys[p]>=y1 && ys[p]<=y2 && xs[p]>=x1 && xs[p]<=x2 && (!d3 || (depths[p]>=D1 && depths[p]<=D2))){
			dz= (d3) ? depths[p]-Depth : 0.0;

			if (r<=cut_sd*SD && dz<=cut_sd*SDd)
			{
				ngridj+=1;
				ngridpointj_temp[ngridj]=p;
				prob[ngridj]= (d3)? exp(-pow(r,2)/(2*pow(SD,2)))*exp(-pow(dz,2)/(2*pow(SDd,2))) : exp(-pow(r,2)/(2*pow(SD,2)));
				probCum+=prob[ngridj];
				if (dAs) A+= dAs[p];
			}
		}
	}

	// allocate memory to *ngridpointj, unless it was allocated before:
	if (!(*ngridpointj)) {
		*ngridpointj=iarray(1,MAX(ngridj,1)+1);
	}
	if (!(*weightsj)) *weightsj=darray(0,MAX(ngridj,1)+1);
	for (p=1; p<=ngridj; p++){
		(*ngridpointj)[p]=ngridpointj_temp[p];
	}

// K*SD is cutoff radius (gaussian would imply inf points), so that total area considered is pi*(K*SD). Vfrac is the fraction of area inside grid (if event is located outside grid).
// Vfrac is fraction of area inside selected area.
	switch (inside) {
	case 0:
	  Vfrac=fabs(A/(pi*(cut_sd*SD)*(cut_sd*SD)));
	  break;
	case 1:
	  Vfrac=1;
	  break;
	default:
	  print_screen("Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
	  print_logfile( "Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
	  return 1;
	  break;
	}

	if (Vfrac>1.1) {
		print_screen("Warning: Vfrac>1 (%lf)\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
		print_logfile("Warning: Vfrac>1 (%lf) \n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
	}
	if (Vfrac>1) Vfrac=1;

	for (p2=1; p2<=ngridj; p2++) (*weightsj)[p2]=Vfrac*prob[p2]/probCum;
	(*weightsj)[0]=1-Vfrac;

	//if no point is selected, select nearest point:
	if (ngridj==0){
		ngridj=1;
		(*ngridpointj)[1]=closestp;
		(*weightsj)[1]=1;
	}

	if (ngridj0) *ngridj0=ngridj;
	free_darray(prob, 0, N+1);
	free_iarray(ngridpointj_temp, 0, N+1);
	return(0);
}

int find_gridpoints_d(double *ys, double *xs, double *depths, int *already_selected, int Nsel0, int N, double y_eq, double x_eq, double Depth, double m,
		double dDCFS, int *ngridj, int **ngridpointj){

/* Similar to find_gridpoints, with two differences: (1) it selects all points in the array already_selected, and it does not calculate weights.
 * The selection radius is magnitude dependent: R~pow(M0/dDCFS), with M0=seismic moment. Exact form was found empirically.
 *
 * Input:
 * 	 ys, xs, depths: arrays with grid points coordinates. Range [1...N]
 * 	 already_selected: sorted indices of points that must be selected. Range [1...Nsel0]
 * 	 y_eq,x_eq, Depth: earthquake coordinates
 *   m: earthquake magnitude.
 *   dDCFS: minimum stress chance that should be resolved (smaller values -> more points selected).
 *
 * Output:
 * 	 ngridj: no of grid points selected.
 * 	 ngridpointj: list of grid points selected.
 */

	int *points_temp;
	double r, rmin=1e30;
	double y1, y2, x1, x2, D1, D2;
	double k=0.85, R;	//k found empirically looking at distrib. of DCFS(r,m).
	int counter=1, closestp;

	points_temp=iarray(1,N);
	*ngridj=0;

	R=pow(k*pow(10.0,3.0*m/2.0)/dDCFS,1.0/3.0);
	y1=y_eq-R;
	y2=y_eq+R;
	x1=x_eq-R;
	x2=x_eq+R;
	D1=fmax(Depth-R, 0.0);
	D2=Depth+R;

	for (int p=1; p<=N; p++){
		if (counter<=Nsel0 && already_selected[counter]==p){
			counter+=1;
			*ngridj+=1;
			points_temp[*ngridj]=p;
		}
		else {
			r=sqrt(pow(ys[p]-y_eq,2)+pow(xs[p]-x_eq,2)+pow(depths[p]-Depth,2));
			if (r<rmin){
				rmin=r;
				closestp=p;
			}
			if (ys[p]>=y1 && ys[p]<=y2 && xs[p]>=x1 && xs[p]<=x2 && depths[p]>=D1 && depths[p]<=D2){
				if (r<=R){
					*ngridj+=1;
					points_temp[*ngridj]=p;
				}
			}
		}
	}

	if (*ngridj==0){
		*ngridj=1;
		*ngridpointj=iarray(1,(*ngridj));
		(*ngridpointj)[1]=closestp;
	}

	else{
		*ngridpointj=iarray(1,(*ngridj));
		for (int p2=1; p2<=*ngridj; p2++) (*ngridpointj)[p2]= points_temp[p2];
	}
	free_iarray(points_temp,1,N);

	return(0);
}

int find_gridpoints_exact(double *ys, double *xs, double *depths, double dx, double dy, double dz, int N, int Nselmax, double y, double x,
		double SD, double Depth, double SDd, int cut_sd, int *ngridj, int **ngridpointj, double **weightsj, int inside, int d3){

/* Similar to find_gridpoints, but instead of simply using center point, integrates over each cell to calculate exact weights.
 * Events are weighted according to a radial Gaussian distribution with a cutoff at SD*cut_sd (SDd*cut_sd vertically).
 * If the source events is close to a boundary, the sum of the weights is <1,since part of the area affected by the earthquake is outside the domain.
 * In this case the sum of the weights equal the ratio between the total grid point area (sum(dAs)) and the area comprised in the cutoff radius.
 *
 *
 * Input:
 *  ys, xs, depths: arrays with grid points coordinates. Range [1...N]
 *  dx, dy, dz: size of cells.
 *  Nselmax: size of ngridpointj, if pre-allocated.
 *  Can also set to *ngridpointj=NULL, so it will be allocated here with right size (in that case set Nselmax>=N to avoid error below).
 *
 *  inside: flag indicating whether the sum of the weights should be set to 1.
 *  d3: flag indicating if the 3D distance should be used.
 *  y,x: earthquake coordinates
 *  SD: standard deviation (horizontal)
 *  Depth: earthquake depth
 *  SDd: standard deviation (depth)
 *
 * Output:
 *  ngridj: no. of selected points
 *  ngridpointj: list of indices of selected points
 *  weightsj: weights of each point
 *
 */

	double r, rmin=1e30, probCum, prob[N+1];
	double rx, ry, rz;
	int p,p2;
	int K, Kd;
	int *ngridpointj_temp;
	int closestp, ngridj_int=0;
	double y1, y2, x1, x2, D1, D2, A, Vfrac=1;

	ngridpointj_temp=iarray(0,N+1);	//temporary array to store indices of selected points.

	if (ngridj) *ngridj=0;
	probCum=0;

// K, Kd determine cutoff radius/depth.
	K=cut_sd;
	Kd=cut_sd;

	y1=y-K*SD;
	y2=y+K*SD;
	x1=x-K*SD;
	x2=x+K*SD;
	D1=Depth-Kd*SDd;
	D2=Depth+Kd*SDd;

	for (p=1; p<=N; p++){
		rx=xs[p]-x;
		ry=ys[p]-y;
		rz= (d3)? depths[p]-Depth : 0.0;
		r= sqrt(pow(ry,2)+pow(rx,2)+pow(rz,2));
		if (r<rmin){
			rmin=r;
			closestp=p;
		}
		if (ys[p]>=y1 && ys[p]<=y2 && xs[p]>=x1 && xs[p]<=x2 && (!d3 || (depths[p]>=D1 && depths[p]<=D2))){
			if (r<=K*SD && (!d3 || rz<=Kd*SDd)){
				ngridj_int+=1;
				if (ngridj_int>Nselmax){
					print_screen("*Error: *ngridj>Nselmax in find_gridpoints.c - need to choose larger value for Nselmax. Exiting. **\n");
					print_logfile("*Error: *ngridj>Nselmax in find_gridpoints.c - need to choose larger value for Nselmax. Exiting. **\n");
					return(1);
				}
				ngridpointj_temp[ngridj_int]=p;
				prob[ngridj_int]= exact_prob(rx,ry,rz,dx, dy,dz,SD, SD, SDd, d3);
				probCum+=prob[ngridj_int];
			}
		}
	}

	A=ngridj_int*dx*dy;

	// allocate memory to *ngridpointj, unless it was allocated before:
	if (!(*ngridpointj)) {
		*ngridpointj=iarray(1,MAX(ngridj_int,1));
	}
	if (!(*weightsj)) {
		*weightsj=darray(0,MAX(ngridj_int,1));
	}
	for (p=1; p<=ngridj_int; p++){
		(*ngridpointj)[p]=ngridpointj_temp[p];
	}


// K*SD is cutoff radius (gaussian would imply inf points), so that total area considered is pi*(K*SD). Vfrac is the fraction of area inside grid (if event is located outside grid).
// Vfrac is fraction of area inside selected area.

	switch (inside) {
		case 0:
		  Vfrac=fabs(A/(pi*(K*SD)*(K*SD)));
		  break;
		case 1:
		  Vfrac=1;
		  break;
		default:
		print_screen("Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
		print_logfile("Error: variable 'inside' has illegal value (%d)! (find_gridpoints).\n");
	  return 1;
		  break;
	}

	if (Vfrac>1.1) {
		print_screen("Warning: Vfrac>1 (%lf)\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
		print_logfile("Warning: Vfrac>1 (%lf) \n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t\n [%lf,%lf]\t%lf+/-%lf]t (find_gridpoints).\n", Vfrac,xs[1],xs[N],x,SD,ys[1],ys[N],y,SD,depths[1],depths[N],Depth,SDd);
	}
	if (Vfrac>1) Vfrac=1;

	for (p2=1; p2<=ngridj_int; p2++) (*weightsj)[p2]=Vfrac*prob[p2]/probCum;
	(*weightsj)[0]=1-Vfrac;

	if (ngridj_int==0){
		ngridj_int=1;
		(*ngridpointj)[1]=closestp;
		(*weightsj)[1]=1;
	}

	if (ngridj) *ngridj=ngridj_int;

	return(0);
}

double exact_prob_1d(double r, double dr, double sd){

/* Integrates a Gaussian distribution between r-0.5dr and r+0.5dr.
 *
 * Input:
 *  r= distance of cell center to central point of gaussian;
 *  dr= cell size;
 *  sd= st.dev. of gaussian.
 *
 * Output:
 *  returns value of integral.
 */

	double 	rmin = r-0.5*dr,	\
			rmax = r+0.5*dr;

	return (erf(rmax/(sqrt(2.0)*sd))-erf(rmin/(sqrt(2.0)*sd)))/2.0;

}

double exact_prob(double rx, double ry, double rz, double dx, double dy, double dz, double sdx, double sdy, double sdz, int d3){
/*
 * Integrates a 3D Gaussian distribution within a cell.
 *
 * Input:
 *  for each dimension:
 *  rX= distance of cell center to central point of gaussian; \
 *  dX= cell size;
 *  sdX= st.dev. of gaussian.
 *  d3= flag to indicate is 3D distance should be used (otherwise, ignore vertical coordinate).
 *
 * Output:
 *   returns value of integral.
 */

	double 	xmin = rx-0.5*dx,	\
			xmax = rx+0.5*dx,	\
			ymin = ry-0.5*dy,	\
			ymax = ry+0.5*dy,	\
			zmin = rz-0.5*dz,	\
			zmax = rz+0.5*dz;

	double Ix, Iy, Iz, I;

	Ix=erf(xmax/(sqrt(2.0)*sdx))-erf(xmin/(sqrt(2.0)*sdx));
	Iy=erf(ymax/(sqrt(2.0)*sdy))-erf(ymin/(sqrt(2.0)*sdy));
	if (d3) Iz=erf(zmax/(sqrt(2.0)*sdz))-erf(zmin/(sqrt(2.0)*sdz));

	I= (d3)? Ix*Iy*Iz/(8*dx*dy*dz) : Ix*Iy/(4*dx*dy);

	return I;

}

int all_nearestneighbours(double *x, double *y, int N, int **pts, double **dist){
	/*
	 * Given a list of points in 2D space, it finds all nearest neighbors.
	 *
	 * Input:
	 *  x, y: coordinates (indices: [1...N]);
	 *  pts: index of nearest neighbour;
	 *
	 * Output:
	 *  dist: distance to nearest neighbour;
	 *  pts, dist are pointers to 1D arrays. If NULL, ignored; if they point to NULL, memory will be allocated. Otherwise, arrays of the correct size should be passed.
	 *
	 * Returns:
	 *  no. of operations (used to test efficiency).
	 */

	double *xs=NULL;	//sorted copies of x,y, (will be sorted).
	int *ind=NULL;
	int *x_ind=NULL;
	int *x_order;
	double d, *dmin= NULL;
	double y_indp;
	int x_toofar, x_toofarahead, x_toofarbehind;
	int indp0, indp;
	int n_op=0;	//no. of times distance is computed.

	if (pts) ind=*pts;
	if (dist) dmin=*dist;
	if (!dmin) {
		dmin=darray(1,N);
		if (dist) *dist=dmin;
	}
	if (!ind) {
		ind=iarray(1,N);
		if (pts) *pts=ind;
	}
	for (int i=1; i<=N; i++) dmin[i]=1e30;

	// sort element by x;

	mysort(N, x, &x_ind, &xs);
	x_order=iarray(1,N);
	for (int i=1; i<=N; i++) x_order[x_ind[i]]=i;

	for (int i=1; i<=N; i++){
		indp0=indp=x_order[i];
		x_toofar=x_toofarahead=x_toofarbehind=0;
		//search in x direction.
		while (!x_toofar){
			indp=(indp>indp0)? 2*indp0-indp : 2*indp0-indp+1;
			if (indp>indp0){
				if (x_toofarahead) continue;
				if (indp>N || fabs(xs[indp]-x[i])>dmin[i]) {
					x_toofarahead=1;
					if (x_toofarbehind) x_toofar=1;
					continue;
				}
			}
			else{
				if (x_toofarbehind) continue;
				if (indp<=0 || fabs(xs[indp]-x[i])>dmin[i]) {
					x_toofarbehind=1;
					if (x_toofarahead) x_toofar=1;
					continue;
				}
			}
			y_indp= y[x_ind[indp]];
			if (fabs(y_indp-y[i])> dmin[i]) continue;
			d=sqrt(pow(xs[indp]-x[i],2)+pow(y_indp-y[i],2));
			n_op+=1;
			if (d<dmin[i]){
				dmin[i]=d;
				ind[i]=x_ind[indp];
			}
			if (d<dmin[x_ind[indp]]){
				dmin[x_ind[indp]]=d;
				ind[x_ind[indp]]=i;
			}
		}
	}

	return n_op;
}

int all_2ndnearestneighbours(double *x, double *y, int N, int **pts, double **dist){
	/*
	 * Given a list of points in 2D space, it finds all 2nd nearest neighbors.
	 *
	 * Input:
	 *  x, y: coordinates (indices: [1...N]);
	 *  pts: index of nearest neighbour;
	 *
	 * Output:
	 *  dist: distance to nearest neighbour;
	 *  pts, dist are pointers to 1D arrays. If NULL, ignored; if they point to NULL, memory will be allocated. Otherwise, arrays of the correct size should be passed.
	 *
	 * Returns:
	 *  no. of operations (used to test efficiency).
	 */

	double *xs=NULL;	//sorted copies of x,y, (will be sorted).
	int **ind=NULL;
	int *x_ind=NULL;
	int *x_order;
	double d, **dmin;
	double y_indp;
	int x_toofar, x_toofarahead, x_toofarbehind;
	int indp0, indp;
	int n_op=0;	//no. of times distance is computed.

	dmin=d2array(1,2,1,N);
	ind=i2array(1,2,1,N);

	for (int i=1; i<=N; i++) dmin[1][i]=dmin[2][i]=1e30;

	// sort element by x;
	mysort(N, x, &x_ind, &xs);
	x_order=iarray(1,N);
	for (int i=1; i<=N; i++) x_order[x_ind[i]]=i;

	for (int i=1; i<=N; i++){
		indp0=indp=x_order[i];
		x_toofar=x_toofarahead=x_toofarbehind=0;
		//search in x direction.
		while (!x_toofar){
			indp=(indp>indp0)? 2*indp0-indp : 2*indp0-indp+1;
			if (indp>indp0){
				if (x_toofarahead) continue;
				if (indp>N || fabs(xs[indp]-x[i])>dmin[2][i]) {
					x_toofarahead=1;
					if (x_toofarbehind) x_toofar=1;
					continue;
				}
			}
			else{
				if (x_toofarbehind) continue;
				if (indp<=0 || fabs(xs[indp]-x[i])>dmin[2][i]) {
					x_toofarbehind=1;
					if (x_toofarahead) x_toofar=1;
					continue;
				}
			}
			y_indp= y[x_ind[indp]];
			if (fabs(y_indp-y[i])> dmin[2][i]) continue;
			d=sqrt(pow(xs[indp]-x[i],2)+pow(y_indp-y[i],2));
			n_op+=1;
			if (d<dmin[2][i] && x_ind[indp]!=ind[1][i]){
				if (d<dmin[1][i]){
					dmin[2][i]=dmin[1][i];
					ind[2][i]=ind[1][i];
					dmin[1][i]=d;
					ind[1][i]=x_ind[indp];
				}
				else {
					dmin[2][i]=d;
					ind[2][i]=x_ind[indp];
				}
			}

			if (d<dmin[2][x_ind[indp]] && ind[1][x_ind[indp]]!=i){
				if (d<dmin[1][x_ind[indp]]){
					dmin[2][x_ind[indp]]=dmin[1][x_ind[indp]];
					ind[2][x_ind[indp]]=ind[1][x_ind[indp]];
					dmin[1][x_ind[indp]]=d;
					ind[1][x_ind[indp]]=i;
				}
				else {
					dmin[2][x_ind[indp]]=d;
					ind[2][x_ind[indp]]=i;
				}
			}
		}
	}

	if (pts) *pts=ind[2];
	if (dist) *dist=dmin[2];

	return n_op;
}
