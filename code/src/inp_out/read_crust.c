
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


#include "read_crust.h"
#include "read_matrix.h"
#include "read_csv.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int read_crust(char *fnametemplate, char *focmecgridfile, struct crust *crst, double resxy, double resz, int vary_focmec){
/*
 * Read crust master file into crst structure.
 *
 * input:
 * 	fname= file containing info about the crust (pscmp or farfalle format).
 * 	fnametemplate= grid file (CSEP format).
 * 	focmecgridfile= grid file containing foc. planes for each grid point.
 *  resxy, resz= resired grid resolution (for calculations);
 *
 * output:
 * 	crst= structure containing info about the domain;
 *
 */

	int err=0, err1=0;
	int NG, ind, NLat, NLon, Nd;
	int no_subpointsx, no_subpointsy, no_subpointsz;
	int no_magbins;
	int is_refined=0;
	double mag1, mag2;
	double dx, dy, dAeq;
	double lat0, lat1,lon0, lon1, d0, d1;
	double *strtmp=0, *diptmp=0, *raketmp=0;
	double str0tmp, dip0tmp;
	double *dumzoneindex=NULL, *zoneindex=NULL;

	print_screen("Loading model setup...");
	print_logfile("\nEntering read_crust...\n");

	double *olats, *olons, *odeps;

	//--------------read grid file:------------------------------//


	//The last column of the file contains indices with indices of focal mechanism zones; it should only be read if (vary_focmec)==1.
	
	err = read_csep_template(fnametemplate, &no_magbins, &((*crst).nLat_out), &((*crst).nLon_out),
							 &((*crst).nD_out), &((*crst).N_allP), &((*crst).dlat_out), &((*crst).dlon_out),
							 &((*crst).ddepth_out), &((*crst).dmags), &olats, &olons, &odeps, (vary_focmec) ? &dumzoneindex : NULL,
							 &((*crst).latmin), &((*crst).latmax), &((*crst).lonmin), &((*crst).lonmax),
							 &((*crst).depmin), &((*crst).depmax), &mag1, &mag2, &((*crst).uniform));

	if(err) {
		print_logfile("Error while reading grid template (%s). Exiting.\n", fnametemplate);
		error_quit(" ** Error while reading grid template. Exiting. **\n");
	}



	lat0=(*crst).latmin;
	lat1=(*crst).latmax;
	lon0=(*crst).lonmin;
	lon1=(*crst).lonmax;
	d0=(*crst).depmin;
	d1=(*crst).depmax;

	(*crst).lat_out=olats;
	(*crst).lon_out=olons;
	(*crst).depth_out=odeps;
	(*crst).lat0=0.5*(lat1+lat0);
	(*crst).lon0=0.5*(lon1+lon0);

	if ((*crst).uniform){
		print_logfile( "Model domain: \n lat=[%.2lf, %.2lf], %d points; \n lon=[%.2lf, %.2lf], %d points; \n dep=[%.2lf, %.2lf], %d points; \n",
					lat0, lat1, (*crst).nLat_out, lon0, lon1, (*crst).nLon_out, d0, d1, (*crst).nD_out);
	}
	else{
		print_logfile( "Model domain: \n lat=[%.2lf, %.2lf]; \n lon=[%.2lf, %.2lf]; \n dep=[%.2lf, %.2lf]; \n",
					lat0, lat1, lon0, lon1, d0, d1);

	}
	print_logfile(" %s grid found.\n", ((*crst).uniform)? "Uniform" : "Non uniform");

	//--------------calculate magnitude bins:-------------//

	(*crst).nmags=no_magbins;
	(*crst).mags=darray(1,no_magbins);

	for (int i=1; i<=no_magbins; i++) (*crst).mags[i]=mag1+(i-1)*(*crst).dmags;

	print_logfile(" mag=[%.2lf, %.2lf], %d bins.\n", (*crst).mags[1], (*crst).mags[(*crst).nmags], (*crst).nmags);

	//--------------calculate refined geometry:-------------//

	dy=Re*DEG2RAD*(*crst).dlat_out;
	dx=Re*DEG2RAD*cos(DEG2RAD*(0.5*(lat1+lat0)))*(*crst).dlon_out;
	no_subpointsx= (int) (0.01+ceil(dx/resxy));
	no_subpointsy= (int) (0.01+ceil(dy/resxy));
	no_subpointsz= (int) (0.01+ceil((*crst).ddepth_out/resz));

	is_refined= (no_subpointsx*no_subpointsy*no_subpointsz>1) & (*crst).uniform;

	if ((*crst).uniform){
		NLat=(*crst).nLat_out;
		NLon=(*crst).nLon_out;
		Nd=(*crst).nD_out;

		NLat*=no_subpointsy;
		NLon*=no_subpointsx;
		Nd*=no_subpointsz;
		(*crst).nLat=NLat;
		(*crst).nLon=NLon;
		(*crst).nD=Nd;
		(*crst).N_allP=NG=(*crst).nLat*(*crst).nLon*(*crst).nD;

		//assume that lat0, lat1 are the boundaries of the domain (*not* the coordinates of the outermost cell centers).
		(*crst).dlat=(lat1-lat0)/(*crst).nLat;
		(*crst).dlon=(lon1-lon0)/(*crst).nLon;
		(*crst).ddepth=(d1-d0)/(*crst).nD;
		(*crst).lat=darray(1,NG);
		(*crst).lon=darray(1,NG);
		(*crst).depth=darray(1,NG);
		(*crst).east=darray(1,NG);
		(*crst).north=darray(1,NG);
		(*crst).dAgrid=darray(1,NG);
		(*crst).list_allP=iarray(1,NG);
		for (int i=1; i<=NG; i++) (*crst).list_allP[i]=i;

		for (int d=1; d<=Nd; d++){
			for (int lo=1; lo<=NLon; lo++){
				for (int la=1; la<=NLat; la++){
					ind=(d-1)*NLat*NLon+(lo-1)*NLat+la;
					(*crst).lat[ind]=lat0+(la-0.5)*(*crst).dlat;
					(*crst).lon[ind]=lon0+(lo-0.5)*(*crst).dlon;
					(*crst).depth[ind]=d0+(d-0.5)*(*crst).ddepth;
				}
			}
		}

		print_logfile("Forecast resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km;\n", dy, dx, (*crst).ddepth_out);
		if (is_refined) {
			print_logfile( "Internal resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km -> %d x %d x %d = %d grid points.\n", resxy, resxy, resz, NLat, NLon, Nd, NG);
		}
		else {
			print_logfile("Internal resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km-> %d x %d x %d = %d grid points.\n", dy, dx, (*crst).ddepth_out, (*crst).nLat_out, (*crst).nLon_out, (*crst).nD_out, (*crst).N_allP);
			print_logfile("Real int.resolution: dlat=%.2lf km, dlon=%.2lf km, ddep=%.2lf km.\n", dy, dx, (*crst).ddepth_out);
		}
	}

	else {
		if (no_subpointsx!=1 || no_subpointsy!=1 || no_subpointsz!=1) {
			print_screen("** Warning: non uniform grid in file %s, can not refine geometry (read_crust.c).**\n",fnametemplate);
			print_logfile("** Warning: non uniform grid in file %s, can not refine geometry (read_crust.c).**\n",fnametemplate);
		}
		(*crst).nLat=0;
		(*crst).nLon=0;
		(*crst).nD=0;
		NG=(*crst).N_allP;

		//assume that lat0, lat1 are the boundaries of the domain (*not* the coordinates of the outermost cell centers).
		(*crst).dlat=(*crst).dlat_out;
		(*crst).dlon=(*crst).dlon_out;
		(*crst).ddepth=(*crst).ddepth_out;
		(*crst).lat=(*crst).lat_out;
		(*crst).lon=(*crst).lon_out;
		(*crst).depth=(*crst).depth_out;
		(*crst).east=darray(1,NG);
		(*crst).north=darray(1,NG);
		(*crst).dAgrid=darray(1,NG);
		(*crst).list_allP=iarray(1,NG);
		for (int i=1; i<=NG; i++) (*crst).list_allP[i]=i;
	}

	//--------------calculate area of each grid cell, and local coordinates:-------------//

	dAeq=pow(Re*PI/180,2)*(*crst).dlon*(*crst).dlat;
	for (int k=1; k<=NG;k++){
		(*crst).dAgrid[k]= dAeq*cos((*crst).lat[k]*PI/180);
		latlon2localcartesian((*crst).lat[k], (*crst).lon[k], (*crst).lat0, (*crst).lon0, (*crst).north+k, (*crst).east+k);
	}

	//--------------read value of focal mechanism grid from file:-------------//
	//----------(this applies is a single foc mec is given per grid point)----//

	//if focmecgridfile is given, each point has its own focal mechanism:
	if (focmecgridfile && strcmp(focmecgridfile,"")!=0){

		//save regional mechanism str0, dip0 (if OOPs are used, these values are use to establish which OOPs should be selected):
		str0tmp=(*crst).str0[0];
		dip0tmp=(*crst).dip0[0];

		//read file of focal mechanisms:
		err1 = read_focmecgridfile(focmecgridfile, crst);

		if(err1){
			print_screen("*Warning: errors occurred while reading focmecgridfile (%s)*\n", focmecgridfile);
			print_logfile("*Warning: errors occurred while reading focmecgridfile*\n");
		}

		// if a refined grid is being used, focal mechanism values should me mapped to refined geometry:
		if (is_refined){
			//allocate strtmp and diptmp since 0th element is also needed:
			strtmp=darray(0,(*crst).N_allP);
			diptmp=darray(0,(*crst).N_allP);

			err1+=convert_geometry((*crst),(*crst).str0, &strtmp, 0, 1);
			(*crst).str0=strtmp;
			err1+=convert_geometry((*crst),(*crst).dip0, &diptmp, 0, 1);
			(*crst).dip0=strtmp;
			//err1+=convert_geometry((*crst),(*crst).rake0, &raketmp, 0, 1);
			//(*crst).rake0=strtmp;
		}

		//copy regional mechanism into 0th element:
		(*crst).str0[0]=str0tmp;
		(*crst).dip0[0]=dip0tmp;

		if(err1){
			print_screen("*Warning: errors occurred while reading focmecgridfile (%s)*\n", focmecgridfile);
			print_logfile("*Warning: errors occurred while reading focmecgridfile*\n");
		}

	}

	//--------------if (vary_focmec==1), read indices and no. of foc. mec zones into crst-----------//

	//NB: vary_focmec and focmecgridfile are mutually exclusive (complains in read_inputfile.c if both files are given).

	else{
		if(vary_focmec) {
			//allocate vector containing indices of zones:
			(*crst).fmzone=iarray(1,(*crst).N_allP);
			//convert indices to refined grid:
			err+=convert_geometry(*crst, dumzoneindex, &zoneindex, 1, 1);
			//indices are in range [0...nofmzones-1]:
			for (int n=1; n<=(*crst).N_allP; n++) 
				{
				(*crst).fmzone[n]= (zoneindex[n]>0) ?  zoneindex[n]-1 : 1;	//max ... needed value is set to 0.0 is the column is missing.
				}
			//the number of zones is equal to the largest index:
			(*crst).nofmzones=(int) max_v(zoneindex+1,(*crst).N_allP);
			//free arrays:
			if (dumzoneindex!=zoneindex) free_darray(dumzoneindex, 1, 1);
			free_darray(zoneindex, 1, 1);
		}
	}

	print_screen("done\n");

	return(err!=0 || err1!=0);
}

int read_focmecgridfile(char *fname, struct crust *crst) {
	/* Reads file containing focal mechanisms for individual grid points into crst structure. The file should use the output geometry.
	 *
	 * Input:
	 * 	fname: 	file name
	 *
	 * Output:
	 *  crst.str0, dip0:	spatially variable receiver faults;
	 *  crst.nofmzones, crst.fmzone:	indices of rec.fault zones (one per grid point; each contains one receiver fault).
	 */


	// Variables used for MPI
	int procId = 0;
	int fileError = 0;

	print_logfile("Reading gridded focal mechanisms from file: %s\n", fname);

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double **data;
	long NL;
	int err=0,
		NP= (*crst).uniform? (*crst).nLat_out*(*crst).nLon_out*(*crst).nD_out : (*crst).N_allP;	//expected no. of lines in file (no. of forecast points).

	data=d2array(1,3,1,NP);
	if (!data) memory_error_quit;	

	if(procId == 0) {
		err = read_matrix(fname, 3, 0, data, &NL);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NL, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&err, 1, MPI_LONG, 0, MPI_COMM_WORLD);

		long nrl=1, nrh=2, ncl=1, nch=(*crst).N_allP;
		long nrow=nrh-nrl+1, ncol=nch-ncl+1;
		MPI_Bcast(data[nrl], nrow*ncol+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	if (err || ((int)NL)!=NP){
		if (err) {
			print_screen("Error: can not open file %s (read_focmecgridfile), exiting.\n", fname);
			print_logfile("Error: can not open file %s (read_focmecgridfile), exiting.\n", fname);
		}
		else {
			print_screen("Error: wrong number of lines in file %s (%d lines found; %d expected). (read_focmecgridfile).\n", fname, NL, NP);
			print_logfile("Error: wrong number of lines in file %s (%d lines found; %d expected). (read_focmecgridfile).\n", fname, NL, NP);
		}

		free_d2array(data, 1,3,1,NP);

		return 1;
	}

	(*crst).variable_fixmec=1;

	//NB zeroth element assigned because it will contain regional mechanism.
	(*crst).str0=darray(0,(*crst).N_allP);
	(*crst).dip0=darray(0,(*crst).N_allP);
	(*crst).rake0=darray(0,(*crst).N_allP);

	for (int i=1; i<=crst->N_allP; i++){
		(*crst).str0[i]=data[1][i];
		(*crst).dip0[i]=data[2][i];
		(*crst).rake0[i]=data[3][i];
	}

	//assign a different zone to each grid point:
	(*crst).nofmzones=(*crst).N_allP;
	(*crst).fmzone=iarray(1,(*crst).N_allP);
	for (int i=1; i<=(*crst).N_allP; i++) (*crst).fmzone[i]=i-1;

	free_d2array(data, 1,2,1,NP);


	return 0;
}

