
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


#include "struct_conversions.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

//------------------ combining -------------------//

int *combine_eqkfm(struct eqkfm *eqkfm1, struct eqkfm *eqkfm2, int N1, int N2,
				   double dt, double dM, double dR, int overwrite) {

	/* Finds common elements of two eqkfm structures, where each member corresponds to an earthquake (i.e. no multifault events).
	 *
	 * Input:
	 *  eqkfm1, eqkfm2:	structures to be combined. range [0,N1-1] and [0,N2-1].
	 *  N1, N2: length of eqkfm1, eqkfm2
	 *  dt, dM, dR: ranges within which earthquakes are considered to be the same.
	 *  overwrite: flag indicating if elements of eqkfm1 should be overwritten with values from eqkfm2.
	 *
	 * Output:
	 *  return array of indices of eqkfm1 corresponding to elements in eqkfm2. Value -1 if event is not found in eqkfm1.
	 *  if (overwrite==1), also changes values in eqkfm1.
	 */

	// Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int n1=0, n10=0, n12=0; //indices of next and previous event and closest in time.
	double dist20, dist2, dlon;
	int selected, *sel, *sel1;
	double dx, dy, dz, r;
	int outfile=0, not_selected=0;	//outfile can be reactivated for printing out selected events (e.g. for debugging).
	char fname[120];
	FILE *fout;

	if (N1==0 | N2==0) {
		print_screen("Warning - one of the eqkfm structures is empty! (combine_eqkfm) \n");
		print_logfile("Warning - one of the eqkfm structures is empty! (combine_eqkfm) \n");
		return NULL;
	}

	sel=iarray(0,N2-1);
	sel1=iarray(0,N1-1);
	for (int n=0; n<N1; n++) sel1[n]=0;

	for (int n2=0; n2<N2; n2++){
		selected=0;

		n10=n12;
		while (n10<N1-1 && eqkfm1[n10].t<=eqkfm2[n2].t-dt) n10++;
		n1=n10;
		while (n1<N1-1 && eqkfm1[n1].t<=eqkfm2[n2].t+dt) n1++;

		dist2=1e30;
		for (int n=n10; n<=n1; n++){
			if (sel1[n]!=0) continue;
			dist20=pow((eqkfm2[n2].t-eqkfm1[n].t)/dt,2)+pow((eqkfm2[n2].mag-eqkfm1[n].mag)/dM,2);
			if (dist2>dist20){
				dist2=dist20;
				n12=n;
			}
		}

		dx=Re*(eqkfm2[n2].lat-eqkfm1[n12].lat)*DEG2RAD;
		dlon=eqkfm2[n2].lon-eqkfm1[n12].lon;
		if (fabs(dlon)>180) dlon=(dlon>0) ? dlon-360 : dlon+360;
		dy=Re*(eqkfm2[n2].lon-eqkfm1[n12].lon)*DEG2RAD*cos(eqkfm1[n12].lat*DEG2RAD);
		dz=eqkfm2[n2].depth-eqkfm1[n12].depth;
		r=sqrt(dx*dx+dy*dy);

		if (fabs(eqkfm2[n2].mag-eqkfm1[n12].mag) <=(dM+0.001) && fabs(eqkfm2[n2].t-eqkfm1[n12].t) <=dt && r<=dR){
			if(procId == 0) {
				if (outfile) fprintf(fout,"%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",eqkfm2[n2].t, eqkfm2[n2].lat, eqkfm2[n2].lon, eqkfm2[n2].depth, eqkfm2[n2].mag, 1, eqkfm1[n12].t, eqkfm1[n12].lat, eqkfm1[n12].lon, eqkfm1[n12].depth, eqkfm1[n12].mag,r);
			}
			if (overwrite==1){
				//copy_eqkfm_nolocation_noindex_notime(eqkfm2[n2], eqkfm1+n12);	//in this case the location is taken from eqkfm2
				copy_eqkfm_noindex_notime(eqkfm2[n2], eqkfm1+n12);	//in this case the location is taken from eqkfm1
			}
			selected+=1;
			sel[n2]=n12;
			sel1[n12]=1;
		}

		else {
			if(procId == 0) {
				if (outfile) fprintf(fout,"%.8lf\t%lf\t%lf\t%lf\t%lf\t%d\t%.8lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",eqkfm2[n2].t, eqkfm2[n2].lat, eqkfm2[n2].lon, eqkfm2[n2].depth, eqkfm2[n2].mag, 0, 0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0);
			}
		}

		if (selected!=1) {
			sel[n2]=-1;
			if (!selected) {
				not_selected+=1;
				if (extra_verbose) {
					print_screen("Warning: element %d [t=%lf, Mw=%lf, d=%.3lf] from eqkfm2 missing in eqkfm1 (function: combine_eqkfm)!!\n",n2,eqkfm2[n2].t,eqkfm2[n2].mag, eqkfm2[n2].depth);
					print_logfile("Warning: element %d [t=%lf, Mw=%lf, d=%.3lf] from eqkfm2 missing in eqkfm1 (function: combine_eqkfm)!!\n",n2,eqkfm2[n2].t,eqkfm2[n2].mag, eqkfm2[n2].depth);
				}
			}
		}
	}
			
	if (not_selected){
		print_screen("Warning: %d elements from focal mechanism catalog missing from earthquake catalog.\n", not_selected);
		print_logfile("Warning: %d elements from focal mechanism catalog missing from earthquake catalog.\n", not_selected);
	}

	if(procId == 0) {
		if (outfile) fclose(fout);
	}

	return sel;
}

int *combine_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM){
	/*
	 * Finds common elements of two earthquake catalogs based on time and magnitudes
	 *
	 * Input:
	 *  t1, m1: event times, magnitudes in catalog 1; range [0...N1-1]
	 *  t2, m2: event times, magnitudes in catalog 2; range [0...N2-1]
	 *  N1, N2: size of tx, mx.
	 *  dt, dM: tolerance
	 *
	 * Output:
	 *  returns array of indices of [t1,m1] corresponding to elements in [t2, m2]. Value -1 if event is not found in [t1, m1].	 *
	 */


	int n1=0, n10=0, n12=0; //indices of next and previous event and closest in time.
	int selected;
	int *sel, *sel1;
	double dist2, dist20;

	if (!N2) return NULL;
	if (!N1) {
		sel= iarray(0,N2-1);
		for (int n=0; n<N2; n++) sel[n]=-1;
		return sel;
	}

	sel= iarray(0,N2-1);
	sel1= iarray(0,N1-1);
	for (int n=0; n<N1; n++) sel1[n]=0;

	for (int n2=0; n2<N2; n2++){
		selected=0;

		n10=n12;
		while (n10<N1-1 && t1[n10]<=t2[n2]-dt) n10++;
		n1=n10;
		while (n1<N1-1 && t1[n1]<=t2[n2]+dt) n1++;

		dist2=1e30;
		for (int n=n10; n<=n1; n++){
			if (sel1[n]!=0) continue;
			dist20=	pow((t2[n2]-t1[n])/dt,2)+pow((m2[n2]-m1[n])/dM,2);
			if (dist2>dist20){
				dist2=dist20;
				n12=n;
			}
		}

		if (fabs(m2[n2]-m1[n12]) <=(dM+0.001) && fabs(t2[n2]-t1[n12]) <=dt){
			selected+=1;
			sel[n2]=n12;
			sel1[n12]=1;
		}
		else sel[n2]=-1;

	}
	return sel;
}

double **union_cats(double *t1, double *t2, double *m1, double *m2, int N1, int N2, double dt, double dM, int ***ind, int *tot){
	/*
	 * Combines two earthquake catalogs based on time and magnitudes.
	 * gives times and magnitude from both (also non common elements).
	 *
	 * Input:
	 *  t1, m1: event times, magnitudes in catalog 1; range [0...N1-1]
	 *  t2, m2: event times, magnitudes in catalog 2; range [0...N2-1]
	 *  N1, N2: size of tx, mx.
	 *  dt, dM: tolerance
	 *
	 * Output:
	 *  ind: contains original indices of the elements in the returned vector (res). range [1-2,0...tot-1].
	 *  value of -1 in ind[x][y] means that element y was not found in one [tx mx], otherwise index is given.
	 *  	 e.g. if (*ind)[1][n]=m, then res[1][m]=t1[n],res[2][m]=m1[n]
	 *  	      if (*ind)[2][n]=m, then res[2][m]=t2[n],res[2][m]=m2[n]
	 *  	      if (*ind)[X][n]=-1, then the element is not found in tX, mX.
	 *
	 *
	 * Returns:
	 *  res: time, magnitude of combined catalog. range [0...tot-1].
	 */

	int n1=0, n10=0, n12=0, n120; //indices of next and previous event, closest and closest to previous element.
	int selected, count=0;
	int *sel, *sel1;
	double dist2, dist20;
	double **res;

	if (ind) *ind=i2array(1,2,0,N1+N2);
	res=d2array(1,2,0,N1+N2);

	if (!N1){
		if (tot) *tot=N2;
		for (int i=0; i<N2; i++){
			if (ind){
				(*ind)[1][i]=-1;
				(*ind)[2][i]=i;
			}
			res[1][i]=t2[i];
			res[2][i]=m2[i];
		}
		return res;
	}

	if (!N2){
		if (tot) *tot=N1;
		for (int i=0; i<N1; i++){
			if (ind){
				(*ind)[1][i]=i;
				(*ind)[2][i]=-1;
			}
			res[1][i]=t1[i];
			res[2][i]=m1[i];
		}
		return res;
	}


	sel=iarray(0,N2-1);
	sel1=iarray(0,N1-1);

	for (int n=0; n<N1; n++) sel1[n]=0;

	for (int n2=0; n2<N2; n2++){
		selected=0;

		n10=n120=n12;
		while (n10<N1-1 && t1[n10]<=t2[n2]-dt) n10++;
		n1=n10;
		while (n1<N1-1 && t1[n1]<=t2[n2]+dt) n1++;

		dist2=1e30;
		for (int n=n10; n<=n1; n++){
			if (sel1[n]!=0) continue;
			dist20=	pow((t2[n2]-t1[n])/dt,2)+pow((m2[n2]-m1[n])/dM,2);
			if (dist2>dist20){
				dist2=dist20;
				n12=n;
			}
		}

		if (fabs(m2[n2]-m1[n12]) <=(dM+0.001) && fabs(t2[n2]-t1[n12]) <=dt){
			selected+=1;
			sel[n2]=n12;
			sel1[n12]=1;
		}
		else sel[n2]=-1;

		int n0=n120;
		while (n0<n12 && t1[n0]<t2[n2]){
			if(sel1[n0]==0) {
				res[1][count]=t1[n0];
				res[2][count]=m1[n0];
				if (ind!= (int ***) 0){
					(*ind)[1][count]=n0;
					(*ind)[2][count]=-1;
				}
				count+=1;
			}
			n0+=1;
		}
		res[1][count]=t2[n2];
		res[2][count]=m2[n2];
		if (ind!= (int ***) 0){
			(*ind)[1][count]=sel[n2];
			(*ind)[2][count]=n2;
		}
		count+=1;
		for	(int n=n0; n<n12; n++){
			if (sel1[n]==0){
				res[1][count]=t1[n];
				res[2][count]=m1[n];
				if (ind!= (int ***) 0){
					(*ind)[1][count]=n;
					(*ind)[2][count]=-1;
				}
				count+=1;
			}
		}
	}

	if (sel1[n12]==1) n120=n12+1;
	else n120=n12;

	for (int n=n120; n<N1; n++){
		res[1][count]=t1[n];
		res[2][count]=m1[n];
		if (ind!= (int ***) 0){
			(*ind)[1][count]=n;
			(*ind)[2][count]=-1;
		}
		count+=1;
	}

	*tot=count;
	return res;
}


double **union_cats2(struct catalog cat, struct pscmp *DCFS, int N2, int ***ind, int *tot){
	/*
	 * Combines two strutures (cat, pscmp) based on the catalog indices contained in pscmp structure.
	 * gives times and magnitude from both (also non common elements).
	 *
	 * Input:
	 *  cat: cat.t, cat.m contain event times, magnitudes in catalog; range [1...cat.Z];
	 *  DCFS: DCFS[x].t, DCFS[x].mag, range x=[0...N2-1]
	 *  N2: size of DCFS.
	 *
	 * Output:
	 *  ind: contains original indices of the elements in the returned vector (res). range [1...2,0...tot-1].
	 *  value of -1 in ind[x][y] means that element y was not found in one [tx mx], otherwise index is given.
	 *  	 e.g. if (*ind)[1][m]=n, then res[1][m]=cat.t[n],res[2][m]=cat.mag[n]
	 *  	      if (*ind)[2][m]=n, then res[1][m]=DCFS[n].t,res[2][m]=DCFS[n].m
	 *  	      if (*ind)[X][m]=-1, then the element is not found in cat, DCFS.
	 *
	 * Returns:
	 *  res: time, magnitude of combined catalog. range [0...tot-1].
	 */

	double t2;
	double **res;
	int N1=cat.Z;
	int c1, c2, k;

	if (ind) *ind=i2array(1,2,0,N1+N2);
	res=d2array(1,2,0,N1+N2);

	if (!N1){
		if (tot) *tot=N2;
		for (int i=0; i<N2; i++){
			if (ind){
				(*ind)[1][i]=-1;
				(*ind)[2][i]=i;
			}
			res[1][i]=DCFS[i].t;
			res[2][i]=DCFS[i].m;
		}
		return res;
	}

	if (!N2){
		if (tot) *tot=N1;
		for (int i=0; i<N1; i++){
			if (ind){
				(*ind)[1][i]=i;
				(*ind)[2][i]=-1;
			}
			res[1][i]=cat.t[i];
			res[2][i]=cat.mag[i];
		}
		return res;
	}

	c1=1;
	c2=k=0;

	while (c1<=cat.Z){
		t2=cat.t[c1];
		//loop over elements of DCFS preceding cat.t[c1]:
		while (c2<N2 && DCFS[c2].t<t2 && DCFS[c2].index_cat<c1){
			res[1][k]=DCFS[c2].t;
			res[2][k]=DCFS[c2].m;
			(*ind)[1][k]=-1;	//element is not in cat;
			(*ind)[2][k]=c2;	//element is in DCFS;
			c2++;
			k++;
		}

		//add elements c1 (which may also be found in DCFS):
		res[1][k]=cat.t[c1];
		res[2][k]=cat.mag[c1];
		(*ind)[1][k]=c1;
		if (c2<N2 && DCFS[c2].index_cat==c1){	//current DCFS element corresponds to this cat element.
			(*ind)[2][k]=c2;
			c2++;
		}
		else (*ind)[2][k]=-1; //element not found in DCFS
		c1++;
		k++;
	}

	//add remaining elements from DCFS:
	while (c2<N2){
		res[1][k]=DCFS[c2].t;
		res[2][k]=DCFS[c2].m;
		(*ind)[1][k]=-1;	//element is not in cat;
		(*ind)[2][k]=c2;	//element is in DCFS;
		c2++;
		k++;
	}

	*tot=k;
	return res;
}


//------------------ filtering -------------------//

void eqk_filter(struct eqkfm **eqkfm1, int *Ntot, double Mag, double Depth){
	/*
	 * Filters eqkfm structure base magnitude and depth.
	 *
	 * Input:
	 * 	eqkfm1 [0...Ntot-1] structure
	 * 	Ntot: size of eqkfm1
	 * 	Mag: selection magnitude (m>=Mag)
	 * 	Depth: depth magnitude (depth<=Depth)
	 *
	 * Output:
	 * 	eqkfm1 is substituted with filtered version.
	 */

	struct eqkfm *eqkfm0;
	int j=0;
	int Ntot_new=0;

	for (int i=0; i<(*Ntot); i++){
		if ((*eqkfm1)[i].mag>=Mag && (*eqkfm1)[i].depth<=Depth) Ntot_new+=1;
	}

	eqkfm0 = eqkfm_array(0, Ntot_new-1);

	for (int i=0; i<(*Ntot); i++){
		if ((*eqkfm1)[i].mag>=Mag && (*eqkfm1)[i].depth<=Depth){
			copy_eqkfm_all((*eqkfm1)[i],eqkfm0+j);
			j+=1;
		}
	}

	*eqkfm1 = eqkfm_array(0, Ntot_new-1);

	for (int i=0; i<Ntot_new; i++) copy_eqkfm_all(eqkfm0[i],(*eqkfm1)+i);
	*Ntot=Ntot_new;
	print_screen("%d events with Mw>=%.3lf, z<=%.3lf selected from eqkfm (eqk_filter).\n",Ntot_new, Mag, Depth);
	print_logfile("%d events with Mw>=%.3lf, z<=%.3lf selected from eqkfm (eqk_filter).\n",Ntot_new, Mag, Depth);
	return;
}

//--------------------extracting 1d arrays------------------------//

double *timesfromeqkfm(struct eqkfm *eqkfm1, int N, int *NF){
	/* Copies times from eqkfm to double vector.
	 *
	 * Input:
	 *  eqkfm1: range [0...nel-1], where nel is the sum of the elements in NF
	 *  N: number of events in eqkfm1 (may be <nel is events are multiple fault events).
	 *  NF: number of faults for each event; range [0, N-1]. If NULL, assume single fault events.
	 *
	 * Output
	 *  returns array of events times; range [0...N-1].
	 */


	double *times=darray(0,N-1);
	int counter=0;

	for (int i=0; i<N; i++) {
		times[i]=eqkfm1[counter].t;
		counter= (NF==NULL)? counter+1 : counter+NF[i];
	}
	return times;

}

double *magssfromeqkfm(struct eqkfm *eqkfm1, int N, int *NF){
	/* Copies magnitudes from eqkfm to double vector.
	 *
	 * Input:
	 *  eqkfm1: range [0...nel-1], where nel is the sum of the elements in NF
	 *  N: number of events in eqkfm1 (may be <nel is events are multiple fault events).
	 *  NF: number of faults for each event; range [0, N-1]. If NULL, assume single fault events.
	 *
	 * Output
	 *  returns array of events magnitudes; range [0...N-1].
	 */


	double *mags=darray(0,N-1);
	int counter=0, NF_i;
	double M0;
	for (int i=0; i<N; i++) {

		M0=0.0;
		NF_i= (NF==NULL)? 1 : NF[i];
		for (int k=0; k<NF_i; k++)	M0+=pow(10,1.5*(eqkfm1[counter+k].mag+6.0));
		mags[i]=(2.0/3.0)*log10(M0)-6.0;
		counter= counter+NF_i;
	}
	return mags;

}

double *timesfrompscmp(struct pscmp *DCFS, int N){
	/* Copies times from DCFS to double vector.
	 *
	 * Input:
	 *  DCFS: range [0...N-1].
	 *  N: size of DCFS
	 *
	 * Output
	 *  returns array of events times; range [0...N-1].
	 */

double *times=darray(0,N-1);
for (int i=0; i<N; i++) times[i]=DCFS[i].t;
return times;

}

double *magsfrompscmp(struct pscmp *DCFS, int N){
	/* Copies magnitudes from DCFS to double vector.
	 *
	 * Input:
	 *  DCFS: range [0...N].
	 *  N: size of DCFS
	 *
	 * Output
	 *  returns array of events magnitudes; range [0...N-1].
	 */

double *mags=darray(0,N-1);
for (int i=0; i<N; i++) mags[i]=DCFS[i].m;
return mags;

}

void eqkfm2dist(struct eqkfm *eqkfm1, double *lats, double *lons, double *depths, int N, int Ntot, int all){
	/*
	 * Calculates distances between sources in eqjfm1 and grid points, and writes them into eqkfm1[i].dist.
	 *
	 * Input:
	 *  eqkfm1:	structure containing sources. Range [0...Ntot-1]
	 *  lats, lons, depths: coordinates of grid points. range [1...N]
	 *  N, Ntot: sizes of grid (lat, lon, depths) and eqkfm structure.
	 *  all: flag indicating if distance should calculated for all elements of eqkfm, or only when eqkfm[i].is_slipmodel==1.
	 *
	 */

	double x,y, *xs, *ys, Depth;
	double lat0, lon0;
	int nsel, pt;

	ys=darray(1,N);
	xs=darray(1,N);
	lat0=0.5*(lats[N]+lats[1]);
	lon0=0.5*(lons[N]+lons[1]);
	for (int k0=1; k0<=N;k0++) latlon2localcartesian(lats[k0], lons[k0], lat0, lon0, ys+k0, xs+k0);

	for (int i=0; i<Ntot; i++){
		if (all==1 || eqkfm1[i].is_slipmodel==0){
			nsel=eqkfm1[i].nsel;
			if (nsel==0) continue;

			latlon2localcartesian(eqkfm1[i].lat, eqkfm1[i].lon, lat0, lon0, &y, &x);
			Depth=eqkfm1[i].depth;
			eqkfm1[i].distance=darray(1,nsel);

			for (int p=1; p<=nsel; p++) {
				pt=eqkfm1[i].selpoints[p];
				eqkfm1[i].distance[p]= sqrt(pow(ys[pt]-y,2)+pow(xs[pt]-x,2)+pow(depths[pt]-Depth,2));
			}
		}
	}
}
