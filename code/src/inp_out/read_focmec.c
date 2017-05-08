
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


#include "read_focmec.h"
#include "read_matrix.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int readmultiplefocmec(char **focmecfiles, int nofiles,
					  struct crust crst, double border, double dz, double dDCFS,
					  struct tm reftime, double t0, double t1, double tfocmec,
					  double mag, double ***focmec,	int **firstelements, int *NFM,
					  int *NFM_timesel, struct eqkfm **eqkfm,int sel, int fm2) {
/*
 * Reads a set of focal mechanisms catalogs and fills in a focmec and a eqkfm structure.
 *
 * Input:
 *  focmecfiles: list of focal mechanisms catalog files [0...nofiles-1]
 *  nofiles: number of focal mechanisms catalogs
 *  crst: structure containing domain information
 *  border, dz: extra (horizontal/vertical) distance to be considered for spatial selection outside of grid in crst.
 *  dDCFS:	min. value for which grid points should be selected for calculating stress changes from source events
 *  reftime: IssueTime (times will be calcualted with reference to this time)
 *  t0,t1: start and end time for selection of focal mechanisms for sources (eqkfm)
 *  tfocmec: end time for selection of focal mechanisms for receiver faults (focmec). Start time is not bounded.
 *  mag: minimum magnitude for selection of focal mechanisms for sources (eqkfm). no lower bound for focmec.
 *  sel: flag indicating is spatial selection should be done.
 *  fm2: flag indicating is both foc. planes should be selected (if fm2=0, only selects first plane).
 *
 * Output:
 *  focmec:	array containing focal mechanisms [1...NC][1...NFM], where NC=no. of columns in file (set below).
 *  first_elements:	indices of focmec elements which correspond to the first element of a new focal mechanism area (i.e. a new focal mechanisms catalog)
 *  	NB firstelements should not be initialized (is done inside this function).
 *  NFM: length of focmec
 *  eqkfm: structure containing sources [0...NFM_timesel-1]
 *  NFM_timesel: length of eqkfm
 *
 *  if focmec=NULL or eqkfm=NULL, they will be ignored.
 *
 */

	double **focmectemp;
	int nfmtemp, nfmtemp2, cl;
	int ntotmax=0, err=0, nfm_sofar=0, nfm_sofar2=0;
	struct eqkfm *eqkfmtemp;

	// Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	if(procId == 0) {
		for(int n=0; n<nofiles; n++) {
			cl=countline(focmecfiles[n]);

			if (cl>=0) {
				ntotmax+=countline(focmecfiles[n]);
			}
			else {
				fileError = 1;
				break;
			}
		}
	}


	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&ntotmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (fm2) ntotmax*=2;

	if(focmec) {
		*focmec = d2array(1,4,1,ntotmax);
		if (!(*focmec)) memory_error_quit;	
	}
	if (firstelements) *firstelements=iarray(0,nofiles);
	if (eqkfm) *eqkfm=eqkfm_array(0,ntotmax-1);

	for (int n=0; n<nofiles; n++) {
		err += readfocmec(focmecfiles[n], crst, border, dz, dDCFS, reftime, t0,
						  t1, tfocmec, mag, (focmec)? &focmectemp : NULL, &nfmtemp, &nfmtemp2,
						  (eqkfm)? &eqkfmtemp : NULL, sel, fm2);
		if (err) break;

		if (focmec && !nfmtemp) {
			print_screen("**Warning: no focal mechanisms selected from file %s. (readmultiplefocmec).**\n", focmecfiles[n]);
			print_logfile("**Warning: no focal mechanisms selected from file %s. (readmultiplefocmec).**\n", focmecfiles[n]);
		}

		if (firstelements) (*firstelements)[n]=nfm_sofar+1;
		if (focmec) {
			for (int ns=1; ns<=4; ns++){
				for (int n=1; n<=nfmtemp; n++) (*focmec)[ns][n+nfm_sofar]=focmectemp[ns][n];
			}
			nfm_sofar+=nfmtemp;
		}
		if (eqkfm){
			for (int i=0; i<nfm_sofar2+nfmtemp2; i++) copy_eqkfm_all(eqkfmtemp[i], (*eqkfm)+nfm_sofar2+i);
			nfm_sofar2+=nfmtemp2;
		}
	}

	if (firstelements) (*firstelements)[nofiles]=nfm_sofar+1;
	if (NFM) *NFM=nfm_sofar;
	if (NFM_timesel) *NFM_timesel=nfm_sofar2;

	return (err>0);
}

int readfocmec(char *focmecfile, struct crust crst,
			   double border, double dz, double dDCFS, struct tm reftime,
			   double t0, double t1, double tfocmec, double mag,
			   double ***focmec, int *NFM, int *NFM_timesel,
			   struct eqkfm **eqkfm,int sel, int fm2) {

/*
 * Reads a focal mechanisms catalog and fills in a focmec and a eqkfm structure.
 *
 * Input:
 *  focmecfile: foc. mec. file
 *  crst: structure containing domain information
 *  border, dz: extra (horizontal/vertical) distance to be considered for spatial selection outside of grid in crst.
 *  dDCFS:	min. value for which grid points should be selected for calculating stress changes from source events
 *  reftime: IssueTime (times will be calcualted with reference to this time)
 *  t0,t1: start and end time for selection of focal mechanisms for sources (eqkfm)
 *  tfocmec: end time for selection of focal mechanisms for receiver faults (focmec). Start time is not bounded.
 *  mag: minimum magnitude for selection of focal mechanisms for sources (eqkfm). no lower bound for focmec.
 *  sel: flag indicating is spatial selection should be done.
 *  fm2: flag indicating is both foc. planes should be selected (if fm2=0, only selects first plane).
 *
 * Output:
 *  focmec:	array containing focal mechanisms [1...NC][1...NFM], where NC=no. of columns in file (set below).
 *  NFM: length of focmec
 *  eqkfm: structure containing sources [0...NFM_timesel-1]
 *  NFM_timesel: length of eqkfm
 *
 *  if focmec=NULL or eqkfm=NULL, they will be ignored.
 */

	// Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	print_logfile("\nReading focal mechanisms from file %s.\n", focmecfile);
	if (fm2) print_logfile("Using both focal mechanisms.\n");
	else print_logfile("Using only first focal mechanism.\n");

	FILE *fin;

	//check if file exists and can be opened.
	if(procId == 0) {
		fin = fopen(focmecfile,"r");
		if(fin == NULL) {
			print_screen("** Error: could not open focal mechanisms catalog %s (readfocmec.c). **\n", focmecfile);
			print_logfile("Error: could not open focal mechanisms catalog (readfocmec.c).\n");
			fileError = 1;
		}
		else {
			fclose(fin);
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return (1);
	}

	int NFMmax;
	double **focmec0;
	int *selected, *selectedsources, p;
	int NC;
	int NFM2=0, NFM2sources=0;
	long NFM0;
	int ignore_time=0;
	int err=0;
	int H, lat_col, lon_col, dep_col, mag_col, str_col, dip_col, rake_col, time_col;	//H=no. of header lines.
	struct tm ev;
	char timestr[20];
	double slip;
	double *times;
	double 	latmin=crst.latmin-180*border/(Re*PI), \
			latmax=crst.latmax+180*border/(Re*PI), \
			lonmin=crst.lonmin-180*border/(Re*PI*cos(crst.lat0*PI/180)), \
			lonmax=crst.lonmax+180*border/(Re*PI*cos(crst.lat0*PI/180)), \
			depthmin=fmin(0.0, crst.depmin-dz), \
			depthmax=crst.depmax+dz;

	int adjust_fault_depth=0;	//flag: if set to 1, synthetic slip models whose top edge is above the surface are shifted down so the top is ad depth=0.


	if(procId == 0) {
		NFMmax = (fm2==1)? 2*countline(focmecfile) : countline(focmecfile);
		NC = countcol_header(focmecfile,1);	//assume one header line (safer than counting columns in header, which may have unwanted spaces);
		if (NC==-1) countcol(focmecfile); //previous function will return -1 is the file has a single line (and no header); in this case, count columns of this line.
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NFMmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NC, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(NFMmax < 0) {
		return 1;
	}

	print_screen("Reading catalog of focal mechanisms...\n");

	lat_col=3;
	lon_col=4;
	dep_col=12;
	mag_col=11;
	str_col=5;	//str, dip, rake have also columns 8, 9, 10.
	dip_col=6;
	rake_col=7;
	time_col=2;
	H=1;

	focmec0=d2array(1,NC,1,NFMmax);
	if (!focmec0) memory_error_quit;	
	selected=iarray(1,NFMmax);
	selectedsources=iarray(1,NFMmax);
	times=darray(1,NFMmax);

	if(procId == 0) {
		err = read_matrix(focmecfile, NC, H, focmec0, &NFM0);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NFM0, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&err, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		long nrl=1, nrh=NC, ncl=1, nch=NFMmax;
		long nrow=nrh-nrl+1, ncol=nch-ncl+1;
		MPI_Bcast(focmec0[nrl], nrow*ncol+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	if (err) return err;


	if(procId == 0) {
		for (int p=1; p<=NFM0; p++){
			if ((sel==0) | ((focmec0[lat_col][p]>=latmin && focmec0[lat_col][p]<=latmax && focmec0[lon_col][p]>=lonmin && focmec0[lon_col][p]<=lonmax && focmec0[dep_col][p]>=(depthmin-1) && focmec0[dep_col][p]<=(depthmax+1)))){
				if (!ignore_time){
					sprintf(timestr, "%.0lf", focmec0[time_col][p]);
					sscanf(timestr, "%4d%2d%2d%2d%2d%2d", &(ev.tm_year), &(ev.tm_mon), &(ev.tm_mday), &(ev.tm_hour), &(ev.tm_min), &(ev.tm_sec));
					ev.tm_year-=1900;
					ev.tm_mon-=1;
					times[p]=difftime(mktime(&ev),mktime(&reftime))*SEC2DAY;
				}		

				if (focmec){
					if (ignore_time || (times[p]<=tfocmec)){
						NFM2+=1;
						selected[NFM2]=p;
					}
				}

				if (eqkfm){
					if (ignore_time || (times[p]>=t0 && times[p]<=t1 && focmec0[mag_col][p]>=mag)){
						NFM2sources+=1;
						selectedsources[NFM2sources]=p;
					}
				}
			}
		}

		print_logfile("%d events selected as sample of focal planes, %d selected as sources.\n", NFM2, NFM2sources);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NFM2, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&NFM2sources, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(times, NFMmax+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(selected, NFMmax+1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(selectedsources, NFMmax+1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	//-------------fill in matrix of focal mechanisms-----------------//

	if (NFM) *NFM= (focmec)? (fm2 ? 2*NFM2 : NFM2) : 0;
	if (focmec){
		*focmec=d2array(1,4,1,*NFM);
		if (!(*focmec)) memory_error_quit;	
		for (int p0=1; p0<=NFM2; p0++){
			p=selected[p0];
			(*focmec)[1][p0]=focmec0[str_col][p];
			(*focmec)[2][p0]=focmec0[dip_col][p];
			(*focmec)[3][p0]=focmec0[rake_col][p];
			(*focmec)[4][p0]= (ignore_time) ? -1e30 : times[p];

			if (fm2==1){
				(*focmec)[1][p0+NFM2]=focmec0[str_col+3][p];
				(*focmec)[2][p0+NFM2]=focmec0[dip_col+3][p];
				(*focmec)[3][p0+NFM2]=focmec0[rake_col+3][p];
				(*focmec)[4][p0+NFM2]=(ignore_time) ? -1e30 : times[p];
			}
		}
	}

	//-------------fill in catalog of sources with foc mec---------------//

	if (NFM_timesel) *NFM_timesel= (eqkfm)? NFM2sources : 0;
	if (eqkfm){
		*eqkfm=eqkfm_array(0,NFM2sources-1);

		#pragma omp parallel for private(p) reduction(+:err)
		for (int p0=0; p0<NFM2sources; p0++){
			p=selectedsources[p0+1];
			(*eqkfm)[p0].t=times[p];
			(*eqkfm)[p0].lat=focmec0[lat_col][p];
			(*eqkfm)[p0].lon=focmec0[lon_col][p];
			(*eqkfm)[p0].depth=focmec0[dep_col][p];
			(*eqkfm)[p0].mag=focmec0[mag_col][p];
			(*eqkfm)[p0].index_cat= 0;
			(*eqkfm)[p0].str1=focmec0[str_col][p];
			(*eqkfm)[p0].dip1=focmec0[dip_col][p];
			(*eqkfm)[p0].rake1=focmec0[rake_col][p];
			if ((*eqkfm)[p0].rake1<0) (*eqkfm)[p0].rake1+=360;
			if (fm2){
				(*eqkfm)[p0].str2=focmec0[str_col][p+3];
				(*eqkfm)[p0].dip2=focmec0[dip_col][p+3];
				(*eqkfm)[p0].rake2=focmec0[rake_col][p+3];
				if ((*eqkfm)[p0].rake2<0) (*eqkfm)[p0].rake2+=360;
				(*eqkfm)[p0].whichfm=0;
			}
			else (*eqkfm)[p0].whichfm=1;

			(*eqkfm)[p0].is_slipmodel=1;
			(*eqkfm)[p0].np_st=1;
			(*eqkfm)[p0].np_di=1;
			(*eqkfm)[p0].pos_s=darray(1,1);	//location of patches within fault; [0], [0] for single patch events.
			(*eqkfm)[p0].pos_d=darray(1,1);
			(*eqkfm)[p0].pos_s[1]=0;	//location of patches within fault; [0], [0] for single patch events.
			(*eqkfm)[p0].pos_d[1]=0;
			if ((*eqkfm)[p0].whichfm){
				(*eqkfm)[p0].slip_str=darray(1,1);
				(*eqkfm)[p0].slip_dip=darray(1,1);
				(*eqkfm)[p0].open=NULL;
			}
			else {
				(*eqkfm)[p0].slip_str=darray(1,2);
				(*eqkfm)[p0].slip_dip=darray(1,2);
				(*eqkfm)[p0].open=NULL;
			}
			err+=find_gridpoints_d(crst.north, crst.east, crst.depth, (int *) 0, 0, crst.N_allP, (*eqkfm)[p0].north, (*eqkfm)[p0].east, (*eqkfm)[p0].depth,  (*eqkfm)[p0].mag, dDCFS,  &((*eqkfm)[p0].nsel), &((*eqkfm)[p0].selpoints));
			WellsCoppersmith((*eqkfm)[p0].mag, (*eqkfm)[p0].rake1, &((*eqkfm)[p0].L), &((*eqkfm)[p0].W), &slip);
			slip=(*eqkfm)[p0].tot_slip[0]=pow(10,(1.5*((*eqkfm)[p0].mag+6)))*(1.0/(crst.mu*pow(10,12)*(*eqkfm)[p0].W*(*eqkfm)[p0].L));
			//shift depth to make sure slip model is not outside domain:
			if ((*eqkfm)[p0].depth<0.5*(*eqkfm)[p0].W*sin(DEG2RAD*(*eqkfm)[p0].dip1)) {
				print_screen("** Warning: synthetic slip model for source event at t=%.3e, mag=%.2lf is partially above ground.\n", (*eqkfm)[p0].t, (*eqkfm)[p0].mag);
				print_logfile("** Warning: synthetic slip model for source event at t=%.3e, mag=%.2lf is partially above ground.\n", (*eqkfm)[p0].t, (*eqkfm)[p0].mag);
			}
			(*eqkfm)[p0].slip_str[1]=slip*cos(DEG2RAD*(*eqkfm)[p0].rake1);
			(*eqkfm)[p0].slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfm)[p0].rake1);

			//by convention, slip_xxx[2] contains the slip for second foc mech (for single patch events only!).
			if (!(*eqkfm)[p0].whichfm){
				if (adjust_fault_depth & (*eqkfm)[p0].depth<0.5*(*eqkfm)[p0].W*sin(DEG2RAD*(*eqkfm)[p0].dip2)) {
					(*eqkfm)[p0].depth=0.5*(*eqkfm)[p0].W*sin(DEG2RAD*(*eqkfm)[p0].dip2);
				}
				(*eqkfm)[p0].slip_str[2]=slip*cos(DEG2RAD*(*eqkfm)[p0].rake2);
				(*eqkfm)[p0].slip_dip[2]=-slip*sin(DEG2RAD*(*eqkfm)[p0].rake2);
			}
	    	latlon2localcartesian((*eqkfm)[p0].lat, (*eqkfm)[p0].lon, crst.lat0, crst.lon0, &((*eqkfm)[p0].north), &((*eqkfm)[p0].east));

		}
	}

	free_d2array(focmec0,1,NC,1,NFMmax);
	free_iarray(selected,1,NFMmax);
	free_iarray(selectedsources,1,NFMmax);
	free_darray(times,1,NFMmax);

	print_screen("done.\n");

	return (err!=0);

}

void select_fm_time(double **focmec, int *NFM, double Tstart){
	/* Selects focal mechanisms occuring before Tstart.
	 *
	 * Input:
	 *  focmec[1...*NFM]: array with focal mechanisms.
	 *
	 * Output:
	 *  overwrites focmec with selected elements. Final Range: [1...*NFM]
	 *  *NFM also overwritten.
	 */

	int ntot=0;

	for (int n=1; n<=*NFM; n++){
		if (focmec[4][n]<Tstart) {
			ntot+=1;
			for (int i=1; i<=4; i++) focmec[i][ntot]=focmec[i][n];
		}
	}

	*NFM=ntot;
}
