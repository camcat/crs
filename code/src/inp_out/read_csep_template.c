
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


#include "read_csep_template.h"

#include <math.h>
#include <stddef.h>

#include "../defines.h"
#include "../geom/convert_geometry.h"

#include "../util/util1.h"
#include "read_matrix.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int read_rate(struct crust crst, char *fname, double **bg_rate, double *r0, double *minmag){
	/* Read 9th column of a forecast-like file (same structure as output grid file), and converts it to internal geometry.
	 *
	 * Input:
	 * 	fname: file name with structure of output file (or forecast template)
	 * 	crst: crst structure used to convert geometry
	 *
	 * Ouput:
	 * 	bg_rate: rate (col. 9 of file), reshaped with internal geometry.
	 * 	r0:	backgroud rate (over entire volume).
	 * 	minmag: min. mag bin found in file.
	 */

	double *dum_rate, dmag;
	double r0int;
	int err;

	err=read_csep_template(fname, 0,0,0,0,0,0,0,0,&dmag,0,0,0,&dum_rate,0,0,0,0,0,0,minmag,0, NULL);
	err+=convert_geometry(crst, dum_rate, bg_rate, 1, 1);
	if (minmag) *minmag-=0.5*dmag;

	r0int=0;
	for (int i=1; i<=crst.N_allP; i++) r0int+= (*bg_rate)[i];	//calculate background rate from file;
	for (int i=1; i<=crst.N_allP; i++) (*bg_rate)[i]*=crst.N_allP/r0int;	//by convention;

	if (r0) *r0=r0int;

	return err;
}

int read_csep_template(char *fname, int *no_magbins, int *nlat, int *nlon,
					   int *ndep, int *ng, double *dlat, double *dlon,
					   double *ddep, double *dmag, double **lats, double **lons,
					   double **deps, double **rate, double *minlat, double *maxlat,
					   double *minlon, double *maxlon, double *mindep, double *maxdep,
					   double *minmag, double *maxmag, int *uni) {

/* Reads a template file in csep format.
 *
 * Input:
 * 	fname, name of txt file.
 *
 * Output:
 * 	no_magbins: no. of magnitude bins;
 * 	nlat, nlon ndep: no of gridpoints with different lat, lon, dep (may not work for non homogeneous grid).
 * 	ng: tot no of grid points.
 * 	dlat, dlon, ddep, dmag: shortst distance between grid points (for lat, lon, this is the "defaultCellDimension" value; for depth, it is calculated).
 * 	lats, lons, deps: 1D arrays containing grid points coordinates: [1...ng].
 * 	minX, maxX: edges of domain (lat, lon, depth);
 * 	minmag, maxmag: smallest, largest *centers* of magnitude bin.
 * 	rate: rate in each cell (summed over magnitude bins).
 *
 * All output variables can be passed as null pointers, and will be ignored.
 *
 * */

	// Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int NC, NL, NH=0, NP;
	int n1, nmag, err=0;
	int miss_col=0;	//used to switch between cell/point (9/6 columns) format.
	int ncolclass=0;
	long NR;
	double **data;
	double lat0=1e30, lat1=-1e30, lon0=1e30, lon1=-1e30, dep0=1e30, dep1=-1e30, mag0, mag1;
	double dlati, dloni, ddepi, dmagi;
	double toll=1e-6;

	if(procId == 0) {
		NC = countcol(fname);
		NL = countline(fname);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NC, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NL, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&miss_col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (NC>=8 && NC<=10) ncolclass=1;	//cells
	else if (NC>=5 && NC<=7) ncolclass=2;	//points


        switch (ncolclass){
                case 1:
                        miss_col=0;
                        if (uni)  *uni=1;
                        break;
		case 2:
                        miss_col=3;
                        if (uni) *uni=0;
                        break;
                
		default:
                        print_screen("Wrong number of columns in template file %s. Exiting.\n", fname);
                        print_logfile("Wrong number of columns in template file %s. Exiting.\n", fname);
                        return 1;
        }


	if (NL>0 && NC>0) {
		data = d2array(1,NC, 1, NL+1);
		if (!data) memory_error_quit;	

		if(procId == 0) {
			print_screen("Reading template file %s\n", fname);
			print_logfile("Reading template file %s\n", fname);
			err = read_matrix(fname, NC, NH, data, &NR);
		}

		#ifdef _CRS_MPI
			MPI_Bcast(&err, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		#endif

		if(err) {
			return 1;
		}

		#ifdef _CRS_MPI
			MPI_Bcast(&NR, 1, MPI_LONG, 0, MPI_COMM_WORLD);

			long nrl=1, nrh=NC, ncl=1, nch=NL+1;
			long nrow=nrh-nrl+1, ncol=nch-ncl+1;
			MPI_Bcast(data[nrl], nrow*ncol+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		#endif
	}
	else {
		return 1;
	}

	nmag=2;
	mag0=data[7-miss_col][1];
	mag1=data[8-miss_col][1];
	while (fabs(data[8-miss_col][1]-data[8-miss_col][nmag])>toll) {
		mag0=fmin(mag0, data[7-miss_col][nmag]);
		mag1=fmax(mag0, data[8-miss_col][nmag]);
		nmag++;
	}
	nmag-=1;

	NP=NR/nmag;
	//assume cells have all the same size.
	
	if (NC>=8){
	  dlati=data[4][1]-data[3][1];
	  dloni=data[2][1]-data[1][1];
	  ddepi=data[6][1]-data[5][1];
	}
	else{
	  dlati=dloni=ddepi=0.0;
	}

	dmagi=data[8-miss_col][1]-data[7-miss_col][1];

	if (lats) *lats=darray(1,NP);
	if (lons) *lons=darray(1,NP);
	if (deps) *deps=darray(1,NP);
	if (rate) *rate=darray(1,NP);

	//closest_lat=1e30;
	//closest_lon=1e30;
	//closest_dep=1e30;

	for (int n=1; n<=NP; n++){
		n1=nmag*n;
		if (NC>=8){
		 	if (lats) (*lats)[n]= 0.5*(data[3][n1]+data[4][n1]);
			if (lons) (*lons)[n]= 0.5*(data[1][n1]+data[2][n1]);
		   	if (deps) (*deps)[n]= 0.5*(data[5][n1]+data[6][n1]);
	                lat0=fmin(lat0, data[3][n1]);
	                lat1=fmax(lat1, data[4][n1]);
	                lon0=fmin(lon0, data[1][n1]);
	                lon1=fmax(lon1, data[2][n1]);
	                dep0=fmin(dep0, data[5][n1]);
	                dep1=fmax(dep1, data[6][n1]);
		}
		else{
		        if (lats) (*lats)[n]=data[2][n1];
                        if (lons) (*lons)[n]=data[1][n1];
                        if (deps) (*deps)[n]=data[3][n1];
                        lat0=fmin(lat0, data[2][n1]);
                        lat1=fmax(lat1, data[2][n1]);
                        lon0=fmin(lon0, data[1][n1]);
                        lon1=fmax(lon1, data[1][n1]);
                        dep0=fmin(dep0, data[3][n1]);
                        dep1=fmax(dep1, data[3][n1]);
		}
		
		if (rate) {
			if (NC==8 | NC==5){
			   //rate column is missing. This is ok for template file, but not for background seismicity file.
			    if (n==1){
				print_screen("Warning: last column not given in file: %s.\n", fname);
				print_logfile("Warning: last column not given in file: %s.\n", fname);
				(*rate)[n]=1;
			   }
			}
			else{
				(*rate)[n]=1;
				for (int i=n1-nmag+1; i<=n1; i++)(*rate)[n]+= data[9-miss_col][i];
			}
		}

	}

        if (nlat) *nlat= (int) (toll+(lat1-lat0)/dlati);
        if (nlon) *nlon= (int) (toll+(lon1-lon0)/dloni);
        if (ndep) *ndep= (int) (toll+(dep1-dep0)/ddepi);

	if (no_magbins) *no_magbins = nmag;
	if (ng) *ng= NP;

	if (dlat) *dlat=dlati;
	if (dlon) *dlon=dloni;
	if (ddep) *ddep=ddepi;
	if (dmag) *dmag=dmagi;

	if (minlat) *minlat=lat0;
	if (maxlat) *maxlat=lat1;
	if (minlon) *minlon=lon0;
	if (maxlon) *maxlon=lon1;
	if (mindep) *mindep=dep0;
	if (maxdep) *maxdep=dep1;
	if (minmag) *minmag=mag0+0.5*dmagi;
	if (maxmag) *maxmag=mag1-0.5*dmagi;

	if(procId == 0) {
		NC = countcol(fname);
		NL = countline(fname);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NC, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NL, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	free_d2array(data,1,NC, 1, NL+1);

	return (0);
}
