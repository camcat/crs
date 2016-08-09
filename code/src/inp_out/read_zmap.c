
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


/*
 *
 *  Created on: Sep 15, 2013
 *      Author: camilla
 */

#include "read_zmap.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int readZMAP (struct catalog *cat, struct eqkfm **eqfm, int *Ntot, char *file,
			  struct crust crst, struct tm reftime, double t0s, double t1s,
			  double t0c, double t1c, double Mmain, double tw, double border,
			  double extra_d, double dDCFS, int findgridpoints) {
/*
 * Reads catalog in zmap format into structure cat (used for LL calculation) and structure eqfm (for stress sources).
 *
 * Input:
 * 	file: zmap catalog.
 * 	crst: structure containing model domain information (used to distribute events across grid points);
 * 	reftime: reference time (IssueTime)
 * 	t0s,t1s: start/end time for considering sources.
 * 	t0c,t1c: start/end for catalog (LL period).
 * 	To avoid using incomplete parts of the catalog, a time window of length tw will be discarded after each event of magnitude >=Mmain in calculating Mc, b values.
 * 		(these events should not be used in the LL calculations).
 * 	border (in km): extra horizontal distance outside of model domain to be considered for sources.
 * 	extra_d (in km): extra vertical distance outside of model domain to be considered for sources.
 * 	findgridpoints: flag indicating if grid points corresponding to events should be calculated (expensive).
 * 	dDCFS: minimum stress value for which grid points should be considered.
 *
 * Output:
 * 	cat: structure containing events to be used for LL calculation (both retrospective and for forecast)
 * 	eqkfm: contains events to be used as sources. It can be given as NULL, and will be skipped.
 * 	Otherwise, everything in these structures is initialized.
 *
 */

	// Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	double std_merr=0.1, std_verr=5.0, std_herr=1.0;	//standard error (km), used if not given in the catalog.
	int line_length=2000;
	int hh, lines=0, valid=0, empty=0, missing_values=0;
	int * old2new;
	FILE *fin;
	char line[line_length], *st;
	int Z;
	double Mc_offset=0.3, dM=0.9;	//Mc_offset since max curvature method tends to underestimate it, and we should not use an incomplete catalog.
	int cut_sd=3.0;	//no. of s.dev. for cutting off gaussian.
	double t_last_large;
	int k;
	struct tm ev;
	int lon_out_of_range, lat_out_of_range, date_out_of_range, time_out_of_range, mag_out_of_range, dep_out_of_range;

	int mon, day, hour, min, sec;
	double fyear, fmon, fday, fhour, flmin, fsec;
	double year, this_year;
	time_t t;

	int eq1=0, eq2=0, eq;
	int errP=0, file_err=0;
	int *seleq1, *seleq2, *catindex;	//catindex: indices of events in catalog, to be copied into eqkfm.
	double SD, SDd, SDlat, SDlon, f=1.0;
	double lat0l, lat1l, lon0l, lon1l, dep0l, dep1l;
	double east, north;
	double Mc;
	double 	*eastgrid=crst.east, \
			*northgrid=crst.north, \
			*depgrid=crst.depth, \
			*dAgrid=crst.dAgrid;
	int N=crst.N_allP;


	//If there is no file given:
	if (strcmp(file,"")==0){
		if (cat){
			(*cat).exists=0;
			(*cat).Z = 0;
			//by convention, Mc>20 means that it should be calculated from the catalog. Bu there is no catalog, hence assume it is complete:
			if ((*cat).Mc>20) (*cat).Mc=-100;
			(*cat).b=1.0;
		}

		if (eqfm) *eqfm=NULL;
		if (Ntot) *Ntot=0;
		return 0;
	}
	else{
		if (cat) (*cat).exists=1;
	}



	if(procId == 0) {
		Z = countline(file)+1;
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&Z, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (Z <= 0) {
		return 1;
	}

	t=time(NULL);
	this_year=(*localtime(&t)).tm_year+1900;

	double 	lon0=crst.lonmin, \
			lon1=crst.lonmax, \
			lat0=crst.latmin, \
			lat1=crst.latmax, \
			dep0=crst.depmin, \
			dep1=crst.depmax;

	double 	*lat=darray(1,Z), \
			*lon=darray(1,Z), \
			*mag=darray(1,Z), \
			*mag2=darray(1,Z), \
			*dep=darray(1,Z), \
			*herr=darray(1,Z), \
			*verr=darray(1,Z), \
			*merr=darray(1,Z), \
			*times=darray(1,Z);  //days from start of catalog (t1c)?

	setenv("TZ", "UTC", 1);
	
	if(procId == 0) {	
		fin = fopen(file,"r");

		if(fin == NULL) {
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		print_screen("**Error: unable to open input file %s (readZMAP).**\n", file);
		print_logfile("**Error: unable to open input file %s (readZMAP).**\n", file);
		return (1);
	}


	//------------------------------read entire catalog:-------------------------//

	int exp_not=0;
	char ex[]="e";

	if(procId == 0) {
		while (!feof(fin)) {
			lines+=1;
			st=fgets(line,line_length,fin);

			//determine if notation is exponential:
			if (lines==1) {
				for (int i=0; i<strlen(line) && !exp_not; i++) {
					if (line[i]==ex[0]) exp_not=1;
				}
			}

			//initialize to out of range values (so will notice if columns are missing):
			lon[valid+1]=999;
			lat[valid+1]=999;
			year=0;
			mon=12;
			day=32;
			mag[valid+1]=999;
			dep[valid+1]=-999;
			hour=25;
			min=62;
			sec=62.0;

			//initialize to standard values (in case they are not given in catalog).
			merr[valid+1]=std_merr;
			verr[valid+1]=std_verr;
			herr[valid+1]=std_herr;

			if (exp_not) {
				hh=sscanf(line, "%25le %25le %25le %25le %25le %25le %25le %25le %25le %25le %25le %25le %25le",
							lon+valid+1, lat+valid+1, &fyear, &fmon, &fday, mag+valid+1, dep+valid+1, &fhour, &flmin, &fsec, herr+valid+1, verr+valid+1, merr+valid+1);

				year=(int) fyear;
				mon=(int) fmon;
				day=(int) fday;
				hour=(int) fhour;
				min= (int) flmin;
				sec= (int) fsec;
			}
			else {
				hh=sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",	lon+valid+1, lat+valid+1, &fyear, &fmon, &fday, mag+valid+1, dep+valid+1, &fhour, &flmin, &fsec, herr+valid+1, verr+valid+1, merr+valid+1);
				year=(int) fyear;
				mon=(int) fmon;
				day=(int) fday;
				hour=(int) fhour;
				min= (int) flmin;
				sec= (int) fsec;
			}

			if (fabs(herr[valid+1])<tol0) herr[valid+1]=std_herr;
			if (fabs(verr[valid+1])<tol0) verr[valid+1]=std_herr;
			if (fabs(merr[valid+1])<tol0) merr[valid+1]=std_herr;

			ev.tm_year=floor(year)-1900;
			ev.tm_mon=mon-1;
			ev.tm_mday=day;
			ev.tm_isdst=0;
			ev.tm_hour=hour;
			ev.tm_min=min;
			ev.tm_sec=(int) sec;

			times[valid+1]=difftime(mktime(&ev),mktime(&reftime))*SEC2DAY;

			//check if catalog is chronological:
			if (valid>=1 && times[valid+1]<times[valid]){
				print_logfile("Error: catalog %s not chronological (ln.%d). Exiting.\n", file, lines);
                                print_screen("Error: catalog %s not chronological (ln.%d). Exiting.\n", file, lines);
				file_err=1;
			}

			if (st==NULL) empty+=1;
			else {
				lon_out_of_range = lat_out_of_range = date_out_of_range = time_out_of_range = mag_out_of_range = dep_out_of_range = 0;

				if (lon[valid+1]<-180 || lon[valid+1]>360) lon_out_of_range=1;
				if (lat[valid+1]<-90 || lat[valid+1]>90) 	 lat_out_of_range=1;
				if (year<1000 || year>this_year ||	mon<1 || mon>12 || day<1 || day>31) date_out_of_range=1;
				if (hour<0 || hour>24 ||	min<0 || min>59 || sec<0 || sec>61) time_out_of_range=1;
				if (mag[valid+1]<-10 || mag[valid+1]>12) mag_out_of_range=1;
				if (dep[valid+1]<0 || dep[valid+1]>2000) dep_out_of_range=1;

				if (lon_out_of_range || lat_out_of_range || date_out_of_range || time_out_of_range || mag_out_of_range || dep_out_of_range) {
					missing_values+=1;

					print_screen("** Warning: line %d has following columns out of range:", lines);
					if (lon_out_of_range) print_screen("lon, ");
					if (lat_out_of_range) print_screen("lat, ");
					if (dep_out_of_range) print_screen("dep, ");
					if (mag_out_of_range) print_screen("mag, ");
					if (date_out_of_range) print_screen("date, ");
					if (time_out_of_range) print_screen("time, ");
					print_screen(" and will be skipped. ** \n");

					print_logfile("Warning: line %d has following columns out of range:", lines);
					if (lon_out_of_range) print_logfile("lon (%.2lf), ", lon[valid+1]);
					if (lat_out_of_range) print_logfile("lat (%.2lf), ", lat[valid+1]);
					if (dep_out_of_range) print_logfile("dep (%.2lf), ", dep[valid+1]);
					if (mag_out_of_range) print_logfile("mag (%.2lf), ", mag[valid+1]);
					if (date_out_of_range) print_logfile("date (%2d-%2d-%4.lf), ", day, mon, year);
					if (time_out_of_range) print_logfile("time (%2d:%2d:%2d). ", hour, min, sec);
					print_logfile(" and will be skipped.\n");
				}
				else valid+=1;
			}
		}
		fclose(fin);

		print_screen("%d of %d lines are valid; %d have missing values; %d are empty\n", valid, lines, missing_values, empty);
		print_logfile("%d of %d lines are valid; %d have missing values; %d are empty\n", valid, lines, missing_values, empty);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&file_err, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(valid == 0){
		(*cat).Z = 0;

		return(1);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(lat, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(lon, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(mag, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(mag2,  Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(dep, 	 Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(herr,  Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(verr,  Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(times, Z+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	//------------------------------select events:-------------------------//

	//by convention, Mc>20 means that Mc should be calculated by program (later).
	if (!cat || (*cat).Mc>=20) Mc=-100;	//select all events.
	else Mc=(*cat).Mc;		//select events>=(*cat).Mc.

	//define large boundaries (for sources):
	lat0l=lat0-180*border/(Re*PI);
	lat1l=lat1+180*border/(Re*PI);
	lon0l=lon0-180*border/(Re*PI*cos(crst.lat0*PI/180));
	lon1l=lon1+180*border/(Re*PI*cos(crst.lat0*PI/180));
	dep0l=fmin(0.0,dep0-extra_d);
	dep1l=dep1+extra_d;

	//count points in catalog, and sources:
	seleq1=iarray(0,valid);
	catindex=iarray(0,valid);
	seleq2=iarray(0,valid);
	old2new=iarray(0,valid);
	t_last_large=-1e30;
	for (int i=1; i<=valid; i++){
		SD=f*herr[i];
		SDd=f*verr[i];
		SDlat=180*SD/(Re*PI);
		SDlon=180*SD/(Re*PI*cos(0.5*(lat0+lat1)*PI/180));

		if (mag[i]>= Mc-dM && lat[i]+SDlat>=lat0l && lat[i]-SDlat<=lat1l && lon[i]+SDlon>=lon0l && lon[i]-SDlon<=lon1l && dep[i]+SDd>=dep0l && dep[i]-SDd<=dep1l){
			if (times[i]>=t0s && times[i]<=t1s) {
				seleq1[eq1]=i;
				eq1+=1;
			}
			if (mag[i]>= Mc && times[i]>=t0c && times[i]<=t1c && lat[i]+SDlat>=lat0 && lat[i]-SDlat<=lat1 && lon[i]+SDlon>=lon0 && lon[i]-SDlon<=lon1 && dep[i]+SDd>=dep0 && dep[i]-SDd<=dep1) {
				if (times[i]>=t_last_large+tw){
					if (times[i]>=t0s && times[i]<=t1s) catindex[eq1-1]=eq2+1;	//+1 since cat.XX[1...cat.Z].
					seleq2[eq2]=i;
					eq2+=1;
				}
				else if (times[i]>=t0s && times[i]<=t1s) catindex[eq1-1]=0;
			}
			else if (times[i]>=t0s && times[i]<=t1s) catindex[eq1-1]=0;
		}
		if (mag[i]>=Mmain) t_last_large=times[i];
	}

	if (!eq1 && !eq2) {
		if (cat) (*cat).Z=0;
		if (Ntot) *Ntot=0;
		print_logfile("No events found in catalog, either as sources or for LL inversion. Exiting.\n");
		print_screen("No events found in catalog, either as sources or for LL inversion. Exiting.\n");
		return 1;
	}

	print_logfile("%d events selected for LL inversion (spatio-temporal selection). \n", eq2);
	print_logfile("%d events selected as sources (spatio-temporal selection). \n", eq1);


	//----------------------------find completeness magnitude and-------------------------//
	//------------------------------select only events with M>Mc--------------------------//


	if (!cat || (*cat).Mc>=20){

		//extract vector with magnitude of selected events:
		for (int i=1; i<=eq2; i++){
			eq=seleq2[i-1];
			mag2[i]=mag[eq];
		}

		//find completeness magnitude using maximum curvature method (Zhuang et al, 2011, Techniques for Analyzing Seismicity Basic models of seismicity: Temporal models");
		(*cat).Mc=Mc_maxcurv(mag2+1, eq2)+Mc_offset;

		//select events about completeness magnitude (for catalog):
		k=0;
		old2new[0]=0;	//events not in seleq2 before are of course still not there.
		for (int i=0; i<eq2; i++){
			if (mag2[i+1]>=(*cat).Mc){
				old2new[i+1]=k+1;		//+1 since cat.XX[1...cat.Z].
				seleq2[k]=seleq2[i];	//doesn't overwrite since k<=i.
				k++;
			}
		}
		eq2=k;

		//select events about completeness magnitude (for eqkfm):
		k=0;
		for (int i=0; i<eq1; i++){
			eq=seleq1[i];
			if (mag[eq]>=(*cat).Mc){
				seleq1[k]=seleq1[i];	//doesn't overwrite since k<=i.
				catindex[k]=old2new[catindex[i]];	//updated indices
				k++;
			}
		}
		eq1=k;

		print_logfile("Calculated completeness magnitude (using maximum curvature): Mc=%.2lf\n", (*cat).Mc);
		print_logfile("%d events selected for LL inversion. \n", eq2);
		print_logfile("%d events selected as sources. \n", eq1);
	}

	//------------------------------fill in catalog:-------------------------//

	if (cat){
		init_cat1(cat, eq2);

		#pragma omp parallel for private(eq, SD, SDd, east, north) reduction(+:errP)
		for (int i=1; i<=eq2; i++){
		if (errP) continue;
		eq=seleq2[i-1];
			SD=f*herr[eq];
			SDd=f*verr[eq];
			(*cat).t[i]=times[eq];
			(*cat).mag[i]=mag[eq];
			(*cat).lat0[i]=lat[eq];
			(*cat).lon0[i]=lon[eq];
			(*cat).depths0[i]=dep[eq];
			(*cat).err[i]=herr[eq];
			(*cat).verr[i]=verr[eq];
			*((*cat).ngrid + i)=0;
			latlon2localcartesian(lat[eq], lon[eq], crst.lat0, crst.lon0, &north, &east);
			(*cat).east0[i]=east;
			(*cat).north0[i]=north;
	    	if (findgridpoints){
				errP+=find_gridpoints(northgrid, eastgrid, dAgrid, depgrid, N, north, east, SD, dep[eq], SDd, cut_sd, (*cat).ngrid + i, &((*cat).ngridpoints[i]), &((*cat).weights[i]), 1, 1);
			}
		}
		(*cat).tstart= ((*cat).Z==0)? t0c : fmax(t0c, (*cat).t[1]);
		(*cat).tend= ((*cat).Z==0)? t1c : fmin(t1c, (*cat).t[(*cat).Z]);

		if (!errP) {
			if ((*cat).Mc>=20) (*cat).Mc=Mc_maxcurv((*cat).mag+1, (*cat).Z)+Mc_offset;
			if ((*cat).Z<100){
				print_logfile("Warning: catalog has fewer than 100 events -> will set b-value to 1.\n");
				print_screen("Warning: catalog has fewer than 100 events -> will set b-value to 1.\n");
				(*cat).b=1.0;
			}
			else{
				(*cat).b=calculatebvalue((*cat).mag+1, (*cat).Z, (*cat).Mc);
			}

			print_screen("Estimated GR values for catalog: Mc=%.2lf, b=%.3lf\n", (*cat).Mc, (*cat).b);
			print_logfile("Estimated GR values for catalog: Mc=%.2lf, b=%.3lf\n", (*cat).Mc, (*cat).b);
		}
	}

	//------------------------------fill in eqkfm:-------------------------//


	if (eqfm){

		if (Ntot) *Ntot=eq1;
		*eqfm=eqkfm_array(0,eq1-1);

		#pragma omp parallel for private(eq) reduction(+:errP)
		for (int i=0; i<eq1; i++){
			if (errP) continue;
			eq=seleq1[i];
			(*eqfm)[i].t=times[eq];
			(*eqfm)[i].lat=lat[eq];
			(*eqfm)[i].lon=lon[eq];
			(*eqfm)[i].depth=dep[eq];
			(*eqfm)[i].mag=mag[eq];
			(*eqfm)[i].index_cat= catindex[i];
	    	latlon2localcartesian((*eqfm)[i].lat, (*eqfm)[i].lon, crst.lat0, crst.lon0, &((*eqfm)[i].north), &((*eqfm)[i].east));
	    	if (findgridpoints){
				if (catindex[i]!=0) errP+=find_gridpoints_d(northgrid, eastgrid, depgrid, (*cat).ngridpoints[catindex[i]], (*cat).ngrid[catindex[i]], N, (*eqfm)[i].north, (*eqfm)[i].east, (*eqfm)[i].depth,  (*eqfm)[i].mag, dDCFS,  &((*eqfm)[i].nsel), &((*eqfm)[i].selpoints));
				else errP+=find_gridpoints_d(northgrid, eastgrid, depgrid, (int *) 0, 0, N, (*eqfm)[i].north, (*eqfm)[i].east, (*eqfm)[i].depth,  (*eqfm)[i].mag, dDCFS,  &((*eqfm)[i].nsel), &((*eqfm)[i].selpoints));
	    	}

		}
	}

	else {
		if (Ntot) *Ntot=0;
	}

	return (errP!=0);
}

int read_firstlineZMAP(char *file, struct tm reftime, double *stime){
/*
 * Reads first line of ZMAP file and return time of the event. Used to determine start time of a catalog (assumed to be chronological).
 *
 * Input:
 * 	file: ZMAP catalog file
 * 	reftime: reference time
 *
 * Output:
 *  stime: time of first event, in days, relative to reftime.
 */


	// Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int line_length=2000;
	FILE *fin;
	char line[line_length], *st;
	int mon, day, hour, min, sec;
	double fyear, fmon, fday, fhour, flmin, fsec;
	double year, this_year;
	struct tm ev;
	time_t t;

	t=time(NULL);
	this_year=(*localtime(&t)).tm_year+1900;

	setenv("TZ", "UTC", 1);

	if(procId == 0) {
		fin = fopen(file,"r");

		if(fin == NULL) {
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		print_screen("**Error: unable to open input file %s (readZMAP).**\n", file);
		print_logfile("**Error: unable to open input file %s (readZMAP).**\n", file);
		return (1);
	}


	//------------------------------read first line if catalog:-------------------------//

	int exp_not=0;
	char ex[]="e";

	if(procId == 0) {
		st=fgets(line,line_length,fin);

		//determine if notation is exponential:
		for (int i=0; i<strlen(line) && !exp_not; i++) {
			if (line[i]==ex[0]) exp_not=1;
		}

		if (exp_not) {
			sscanf(line, "%*25le %*25le %25le %25le %25le %*25le %*25le %25le %25le %25le %*25le %*25le %*25le",
						&fyear, &fmon, &fday, &fhour, &flmin, &fsec);

			year=(int) fyear;
			mon=(int) fmon;
			day=(int) fday;
			hour=(int) fhour;
			min= (int) flmin;
			sec= (int) fsec;
		}
		else {
			sscanf(line, "%*lf %*lf %lf %lf %lf %*lf %*lf %lf %lf %lf %*lf %*lf %*lf\n",	&fyear, &fmon, &fday, &fhour, &flmin, &fsec);
			year=(int) fyear;
			mon=(int) fmon;
			day=(int) fday;
			hour=(int) fhour;
			min= (int) flmin;
			sec= (int) fsec;
		}

		ev.tm_year=floor(year)-1900;
		ev.tm_mon=mon-1;
		ev.tm_mday=day;
		ev.tm_isdst=0;
		ev.tm_hour=hour;
		ev.tm_min=min;
		ev.tm_sec=(int) sec;

		*stime=difftime(mktime(&ev),mktime(&reftime))*SEC2DAY;
	}

	#ifdef _CRS_MPI
		MPI_Bcast(stime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	return 0;
}

