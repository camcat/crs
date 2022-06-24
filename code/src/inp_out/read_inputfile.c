
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


#include "read_inputfile.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../defines.h"
#include "../util/error.h"

#include "../util/util1.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int read_inputfile(char *input_fname, char *outname, char *fore_template,
		char *catname, char ***focmeccat, char *background_rate_grid, char *background_rate_cat, char *fixedmecfile, char *slipmodelfile, char *aseismicmodelfile,
		char *model_parameters_file, char *Logfile, struct tm *reftime,
		double *Tstart, double *Tend, double *tstartLL, double *tendLL, long *seed, int *num_fm){

	// todo: could check string length is enough
	/* Read master input file.
	 *
	 * input: file name input_fname
	 *
	 * output:
	 * 		outname: output file name (without extension)
	 * 		reftime: structure containing the reference time (issue time)
	 * 		Tstart, Tend: time in days from reftime
	 * 		fore_template: forecast template
	 * 		catname: catalog (ZMAP)
	 * 		focmeccat: catalog of focal mechanisms
	 * 		background_rate_grid: file containing background seismicity model
	 * 		slipmodefile, aseismicmodelfile: files containing a list of slip models/aseismicmodel snapshots
	 * 		model_parameters_file: file containing model parameters
	 * 		Logfile: log file to be written
	 * 		seed: seed for random number generator
	 * 		num_fm: number of focal mechanisms areas
	 *
	 * NB: all pointers will be ignored if NULL is passed; otherwise, char* should already be initialized.
	 */

	setenv("TZ", "UTC", 1);

	// Variables used for MPI
	int procId = 0;
	int fileError = 0;
	int bCastFocmeccat = 0;
	int focmeccatIsNull = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	FILE *fin;
	int Nchar=1000;
	char line[Nchar], listfocmeccat[Nchar];
	char *key, *value;
	int NP=19, i, err=0;
	struct tm times0, times1, times2, times3, times4;
	int value_found[NP];
	int listfm=0, nofm=0;	//listfm is a flag indicating whether multiple focal mechanism catalogs are given. nofm is the number of such catalogs.
	char comment[]="#", comm=comment[0];

	for (int n=0; n<NP; n++) value_found[n]=0;
	if (tendLL) *tendLL=0.0;	//Default value: use time span up to IssueTime.

	char *keys[]={
	/*0*/	"IssueDate",
	/*1*/	"ForecastStartDate", \
	/*2*/	"ForecastEndDate",	\
	/*3*/	"OutputForecastFile", \
	/*4*/	"InputCatalogFile",	\
	/*5*/	"InputCatalogFocMecFile",	\
	/*6*/	"InputListCatalogFocMecFile",	\
	/*7*/	"ForecastTemplate", \
	/*8*/	"", \
	/*9*/	"InputListSlipModels", \
	/*10*/	"InputListAseismicModels", \
	/*11*/	"InputBackgroundRateGrid",\
	/*12*/	"InputModelParametersFile",\
	/*13*/	"RandomSeedValue",\
	/*14*/	"Logfile",\
	/*15*/	"FixedMecFile", \
	/*16*/	"InputBackgroundRateCatalog", \
	/*17*/	"InversionStartDate", \
	/*18*/	"InversionEndDate"
	};


	// NB: arguments 5,6 are alternative (one or single receiver fault catalog).

	// If there is a file error, only root will know about it.
	// So it is important that the error is broadcast to all
	// processes, and they all return from the function with
	// the same error code.
	if(procId == 0) {
		fin = fopen(input_fname, "r");
		if(fin == NULL) {
			print_screen("Error read_input: unable to open input file %s.\n", input_fname);
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	if(procId == 0) {
		while(!feof(fin)) {
			fgets(line,Nchar,fin);
			if (line[0]==comm) continue;
			if (ferror(fin)) {
				error_quit("Error reading input data using fgets!\n");
				fileError = 1;
				break;
			}
			key=strtok(line,"=");
			value=strtok(NULL,"=");
			if (!value) continue;
			i=0;
			while (i<NP && strcmp(key,keys[i])) i++;
			if (i>=NP){
				print_screen("Error read_inputfile: parameter \" %s\" in file \"%s\" not recognized.\n", key, input_fname);
				print_logfile("Error read_inputfile: parameter \" %s\" in file \"%s\" not recognized.\n", key, input_fname);
				fileError = 1;
				break;
			}

			value_found[i]=1;

			// Fill in output variables with values from file:
			switch(i){
				case 0:
					sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times0.tm_year), &(times0.tm_mon), &(times0.tm_mday), &(times0.tm_hour), &(times0.tm_min), &(times0.tm_sec));
					times0.tm_year-=1900;
					times0.tm_mon-=1;
					times0.tm_isdst=0;
					if (reftime) *reftime=times0;
					if (Tstart && value_found[1]) *Tstart=difftime(mktime(&times1),mktime(&times0))*SEC2DAY;
					if (Tend && value_found[2]) *Tend=difftime(mktime(&times2),mktime(&times0))*SEC2DAY;
					if (tstartLL && value_found[17]) *tstartLL=difftime(mktime(&times3),mktime(&times0))*SEC2DAY;
					if (tendLL && value_found[18]) *tendLL=difftime(mktime(&times4),mktime(&times0))*SEC2DAY;
					break;

				case 1:
					sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times1.tm_year), &(times1.tm_mon), &(times1.tm_mday), &(times1.tm_hour), &(times1.tm_min), &(times1.tm_sec));
					times1.tm_year-=1900;
					times1.tm_mon-=1;
					times1.tm_isdst=0;
					if (Tstart && value_found[0]) *Tstart=difftime(mktime(&times1),mktime(&times0))*SEC2DAY;
					break;
				case 2:
					sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times2.tm_year), &(times2.tm_mon), &(times2.tm_mday), &(times2.tm_hour), &(times2.tm_min), &(times2.tm_sec));
					times2.tm_year-=1900;
					times2.tm_mon-=1;
					times2.tm_isdst=0;
					if (Tend && value_found[0]) *Tend=difftime(mktime(&times2),mktime(&times0))*SEC2DAY;
					break;
				case 3:
					if (outname) sscanf(value,"%s",outname);
					break;
				case 4:
					if (catname) sscanf(value,"%s",catname);
					break;
				case 5:
					if (focmeccat) {
						// In this case a single focal mechanism catalog is given:
						// the file name can be read directly into the 0th element of focmeccat.
						*focmeccat= malloc(sizeof(char*));
						(*focmeccat)[0]=malloc(120*sizeof(char));
						sscanf(value,"%s",(*focmeccat)[0]);

						bCastFocmeccat = 1;	// Flag, indicating that 'focmeccat' should be broadcast.
					}
					if (num_fm) *num_fm=1;
					break;
				case 6:
					// In this case several catalog of focal mechanisms are given:
					// the names of the catalogs, listed in the file "listfocmeccat" will be read later on.
					if (focmeccat) sscanf(value,"%s",listfocmeccat);
					listfm=1;
					break;
				case 7:
					if (fore_template) sscanf(value,"%s",fore_template);
					break;
				case 8:
					// just empty line
					break;
				case 9:
					if (slipmodelfile) sscanf(value,"%s",slipmodelfile);
					break;
				case 10:
					if (aseismicmodelfile) sscanf(value,"%s",aseismicmodelfile);
					break;
				case 11:
					if (background_rate_grid) sscanf(value,"%s",background_rate_grid);
					break;
				case 12:
					if (model_parameters_file) sscanf(value,"%s",model_parameters_file);
					break;
				case 13:
					if (seed) sscanf(value,"%ld",seed);
					break;
				case 14:
					if (Logfile) sscanf(value,"%s",Logfile);
					break;
				case 15:
					if (fixedmecfile) sscanf(value,"%s",fixedmecfile);
					break;
				case 16:
					if (background_rate_cat) sscanf(value,"%s",background_rate_cat);
					break;
				case 17:
					sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times3.tm_year), &(times3.tm_mon), &(times3.tm_mday), &(times3.tm_hour), &(times3.tm_min), &(times3.tm_sec));
					times3.tm_year-=1900;
					times3.tm_mon-=1;
					times3.tm_isdst=0;
					if (tstartLL && value_found[0]) *tstartLL=difftime(mktime(&times3),mktime(&times0))*SEC2DAY;
					break;
				case 18:
					sscanf(value, "%d-%d-%dT%d:%d:%dZ", &(times4.tm_year), &(times4.tm_mon), &(times4.tm_mday), &(times4.tm_hour), &(times4.tm_min), &(times4.tm_sec));
					times4.tm_year-=1900;
					times4.tm_mon-=1;
					times4.tm_isdst=0;
					if (tendLL && value_found[0]) *tendLL=difftime(mktime(&times4),mktime(&times0))*SEC2DAY;
					break;
			}
		}

		fclose(fin);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}


	#ifdef _CRS_MPI
		// The file names are used in conditions in main.c for
		// setting certain flags. catname is used in setup.c.
		MPI_Bcast(Tstart, 				 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(Tend, 				 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(tstartLL, 			 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(tendLL, 			 	 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(seed, 				 1,   MPI_LONG,   0, MPI_COMM_WORLD);
		MPI_Bcast(num_fm, 				 1,   MPI_INT, 	  0, MPI_COMM_WORLD);
		MPI_Bcast(&bCastFocmeccat, 		 1,   MPI_INT, 	  0, MPI_COMM_WORLD);
		MPI_Bcast(&listfm, 				 1,   MPI_INT, 	  0, MPI_COMM_WORLD);

		if(bCastFocmeccat) {
			if(procId != 0) {
				*focmeccat = malloc(sizeof(char*));
				(*focmeccat)[0] = malloc(120*sizeof(char));
			}

			MPI_Bcast((*focmeccat)[0], 120, MPI_CHAR, 0, MPI_COMM_WORLD);
		}
	#endif

	if (focmeccat && listfm) {
		err = read_focmecfiles(listfocmeccat, focmeccat, num_fm);
		if (err) {
			print_screen("Error: could not read file %s.\n", listfocmeccat);
			print_logfile("Error: could not read file %s.\n", listfocmeccat);
		}
	}

	// Print out warnings or errors for missing or redundant parameters:
	if(procId == 0) {

		if ((value_found[6] + value_found[5])>1){
			print_screen("Error: parameters %s, %s are alternative to each other. Exit.\n", keys[5], keys[6]);
			print_logfile("Error: parameters %s, %s are alternative to each other. Exit.\n", keys[5], keys[6]);
			err=1;

		}
		if ((value_found[11] + value_found[16])>1){
			print_screen("Error: parameters %s, %s are alternative to each other. Exit.\n", keys[11], keys[16]);
			print_logfile("Error: parameters %s, %s are alternative to each other. Exit.\n", keys[11], keys[16]);
			err=1;
		}
		nofm=0;
		for (int n=0; n<NP; n++) {
			if (!value_found[n]) {
				switch (n){
				case 3:
					if (extra_verbose) {
						print_screen("Warning: parameter %s not given in %s -> will use default value 'output'.\n", keys[n], input_fname);
						print_logfile("Warning: parameter %s not given in %s -> will use default value 'output'.\n", keys[n], input_fname);
					}
					if (outname) strcpy(outname,"output");
					break;

				case 4:
					if (extra_verbose) {
						print_screen("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
						print_logfile("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
					}
					if (catname) strcpy(catname,"");
					break;

				case 5:
					if (nofm){
						if (extra_verbose) {
							print_screen("Warning: parameters %s, %s not given in %s.\n", keys[5], keys[6], input_fname);
							print_logfile("Warning: parameters %s, %s not given in %s.\n", keys[5], keys[6], input_fname);
						}
						if (focmeccat) {
							*focmeccat=NULL;

							focmeccatIsNull = 1;
						}
						if (num_fm) *num_fm=1;
					}
					else nofm=1;
					break;
				case 6:
					if (nofm){
						if (extra_verbose) {
							print_screen("Warning: parameters %s, %s not given in %s.\n", keys[5], keys[6], input_fname);
							print_logfile("Warning: parameters %s, %s not given in %s.\n", keys[5], keys[6], input_fname);
						}
						if (focmeccat) {
							*focmeccat=NULL;

							focmeccatIsNull = 1;
						}
						if (num_fm) *num_fm=1;
					}
					else nofm=1;
					break;
				case 8:
					//just empty line
					break;
				case 9:
					if (extra_verbose) {
						print_screen("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
						print_logfile("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
					}
					if (slipmodelfile) strcpy(slipmodelfile,"");
					break;
				case 10:
					if (extra_verbose) {
						print_screen("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
						print_logfile("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
					}
					if (aseismicmodelfile) strcpy(aseismicmodelfile,"");
					break;
				case 11:
					if (extra_verbose) {
						print_screen("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
						print_logfile("Warning: parameter %s not given in %s.\n", keys[n], input_fname);
					}
					if (background_rate_grid) strcpy(background_rate_grid,"");
					break;
				case 13: 
					print_screen("Warning: parameter %s not given in %s, will use default value.\n", keys[n], input_fname);
					print_logfile("Warning: parameter %s not given in %s, will use default value.\n", keys[n], input_fname);
					*seed=-37284630;
					break;
				case 14:
					if (Logfile) strcpy(Logfile,"");
					break;

				case 15:
					if (fixedmecfile) strcpy(fixedmecfile,"");
					break;
				case 16:
					if (background_rate_cat) strcpy(background_rate_cat,"");
					break;
				case 17:
					if (tstartLL) *tstartLL= NAN;	//this is used later (will give error if LLinversion=1 in parameter file).
					break;
				case 18:
					print_logfile("InversionEndDate not given in %s, will use IssueTime as end of LL inversion period.\n", input_fname);
					if (tendLL) *tendLL=0.0;
					break;

				default:
					print_screen("Error: parameter %s not given in %s.\n", keys[n], input_fname);
					print_logfile("Error: parameter %s not given in %s.\n", keys[n], input_fname);
					err=1;
					break;
				}
			}
		}
	}

	#ifdef _CRS_MPI

		// Broadcast file names:
		// even though the files will only be used by procId=0, all ranks may use them to set flags (e.g. if they are not given).
		MPI_Bcast(catname,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(outname,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(slipmodelfile,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(aseismicmodelfile,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(background_rate_grid,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(background_rate_cat,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);
		MPI_Bcast(fixedmecfile,  			 120, MPI_CHAR,   0, MPI_COMM_WORLD);

		MPI_Bcast(num_fm, 				1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nofm, 				1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&focmeccatIsNull, 	1, MPI_INT, 0, MPI_COMM_WORLD);
		// Broadcase these again since they may have changed.
		MPI_Bcast(tstartLL, 			 	 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(tendLL, 			 	 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(seed, 				 1,   MPI_LONG,   0, MPI_COMM_WORLD);

		if(procId != 0) {
			if(focmeccatIsNull) {
				*focmeccat=NULL;
			}
		}
	#endif

	return (err || fileError);
}

int read_focmecfiles(char *inputfile, char ***listfiles, int *nfiles) {

/* Read a file containing a list of focal mechanism catalogs files (corresponding to different areas).
 *
 * input:
 *  inputfile
 *
 * output:
 *  listfiles, a list of filenames;
 *  nfiles, the number of files.
 *
 */

	int Nchar=1000;
	char line[Nchar];
	char comment[]="#", comm=comment[0];
	FILE *fin;

	// Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	if(procId == 0) {
		if (!(fin=fopen(inputfile,"r"))) {
			print_screen("Error: can not open file %s (read_focmecfiles), Exiting.\n", inputfile);
			print_logfile("Error: can not open file %s (read_focmecfiles), Exiting.\n", inputfile);
			fileError = 1;
		}
	}
	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		*listfiles = NULL;

		return 1;
	}

	if(procId == 0) {
		line[0]=comm;
		while (line[0]==comm)fgets(line,Nchar,fin);
		if (ferror(fin)) print_screen("Error reading input data (file: %s) using fgets!\n", inputfile);
		sscanf(line,"%d", nfiles);

		*listfiles = malloc((*nfiles)*sizeof(char*));
		for (int nn=0; nn<(*nfiles); nn++) {
			(*listfiles)[nn] = malloc(120 * sizeof(char));
			line[0]=comm;
			while (line[0]==comm)fgets(line,Nchar,fin);
			if (ferror(fin)) {
				print_screen("Error reading input data (file: %s) using fgets!\n", inputfile);

				fileError = 1;

				break;
			}
			sscanf(line,"%s", (*listfiles)[nn]);
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(nfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// The following will generate an error if ferror(fin)
		// was true in any of the iterations of the while loop above.
		// The print statements above should help with debugging
		if(procId != 0) {
			*listfiles = malloc((*nfiles)*sizeof(char*));
			for(int i = 0; i < (*nfiles); ++i) {
				(*listfiles)[i] = malloc(120 * sizeof(char));
			}
		}

		for(int i = 0; i < (*nfiles); ++i) {
			MPI_Bcast((*listfiles)[i], 120, MPI_CHAR, 0, MPI_COMM_WORLD);
		}
	#endif

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	return 0;
}

// TODO: MPI code this this function requires optimization ...
int read_listslipmodel(char *input_fname, struct tm reftime, struct slipmodels_list *allslipmodels,
					   double res, int is_aseismic, int *aseismic_linear, double *t0log, int *flag_multisnap) {
	/*
	 * Read a file containing a list of slip models.
	 *
	 * input:
	 * 	input_fname (file name)
	 * 	reftime: structure containing the value of corresponding to reference time (IssueTime).
	 * 	is_aseismic: flag indicating whether the file contains coseismic slip models or aseismic slip (files are organized differently).
	 * 	res: desired slip model resolution
	 *
	 * output:
	 * 	allslipmodels: list of slip models;
	 * 	aseismic_log: flag indicating if logarithmic time steps should be used for aseismic processes (splines or single snapshot).
	 *
	 */

	// Variables used for MPI
	int procId = 0;
	int fileError = 0, size_slipmodels = 0;
	int aseismic_log=0, aseismic_splines=0;	//flags.

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	FILE *fin;
	int Nchar=1000;
	char line[Nchar];
	char time_str[50], log_evol[10];
	char comment[]="#", comm=comment[0];
	struct tm times;
	int Nm0, nsm, no_slipmod;

	//if input_fname=="", this was not given in the input file. Variables should be set to NULL.
	if (strcmp(input_fname,"")==0){
		print_screen("Warning: no slip model file given input file.\n");
		print_logfile("\nWarning no slip model file given input file.\n");
		(*allslipmodels).NSM=0;
		(*allslipmodels).is_aseismic=is_aseismic;
		(*allslipmodels).tmain= NULL;
		(*allslipmodels).mmain= NULL;
		(*allslipmodels).disc= NULL;
		(*allslipmodels).Nfaults=NULL;
		(*allslipmodels).no_slipmodels=NULL;
		(*allslipmodels).is_used=NULL;
		(*allslipmodels).cut_surf=NULL;
		return 0;
	}

	if(procId == 0) {
		fin = fopen(input_fname, "r");
		if(fin == NULL) {
			fileError = 1;
		}
	}

	if (is_aseismic) *flag_multisnap=-1;	//initial value, meaning not assigned yet. (will change later)


	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		print_screen("Warning: slip model file %s not found (read_listslipmodel).\n", input_fname);
		print_logfile("\nWarning: slip model file %s not found (read_listslipmodel).\n", input_fname);
		(*allslipmodels).NSM=0;
		(*allslipmodels).is_aseismic=is_aseismic;
		(*allslipmodels).tmain= NULL;
		(*allslipmodels).mmain= NULL;
		(*allslipmodels).disc= NULL;
		(*allslipmodels).Nfaults=NULL;
		(*allslipmodels).no_slipmodels=NULL;
		(*allslipmodels).is_used=NULL;
		(*allslipmodels).cut_surf=NULL;
		return 1;
	}

	else {
		line[0]=comm;
		if(procId == 0) {
			while (line[0]==comm)fgets(line,Nchar,fin);
			if (ferror(fin)) fprintf(stderr, "ERROR reading input data (file: %s) using fgets!\n", input_fname);
			sscanf(line,"%d", &Nm0);
			fgets(line,Nchar,fin);
//			(*allslipmodels).cmb_format=(char *)malloc(120*sizeof(char));
			if (is_aseismic) {
				sscanf(line,"%s %s %lf", &((*allslipmodels).cmb_format), log_evol, t0log);
				if (!strcmp(log_evol, "log")) {
					*aseismic_linear=0;
					aseismic_log=1;
				}

				else {
					if (!strcmp(log_evol, "splines")) {
						*aseismic_linear=0;
						aseismic_splines=1;
					}
					else {
						if (!strcmp(log_evol, "lin")) {

							*aseismic_linear=1;
						}
						else{
							print_screen("Illegal value for type of temporal evolution (should be one of: lin/log).\n", input_fname);
							print_logfile("Illegal value for type of temporal evolution (should be one of: lin/log.\n", input_fname);
							fileError=1;
						}
					}
				}
			}
			else {
				sscanf(line,"%s %d", &((*allslipmodels).cmb_format), &((*allslipmodels).constant_geometry));
			}
		}

		#ifdef _CRS_MPI
			MPI_Bcast(&Nm0, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (is_aseismic){
				MPI_Bcast(&aseismic_log, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&aseismic_splines, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(aseismic_linear, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(t0log, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
			MPI_Bcast(&((*allslipmodels).constant_geometry), 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&((*allslipmodels).cmb_format), 120, MPI_CHAR, 0, MPI_COMM_WORLD);
		#endif

		#ifdef _no_numerical_recipes
		    if (aseismic_splines){
			print_logfile(" Error: can not use splines if Numerical Recipes source code is not available.");
			print_screen(" Error: can not use splines if Numerical Recipes source code is not available.");
			return(1);
		    }
		#endif

		(*allslipmodels).NSM=Nm0;
		(*allslipmodels).is_aseismic=is_aseismic;
		(*allslipmodels).tmain= darray(0,Nm0-1);	//-1 element to store mainshock time (when aseismic slip starts).
		(*allslipmodels).mmain= (is_aseismic)? NULL : darray(0,Nm0-1);
		(*allslipmodels).cut_surf=malloc(Nm0*sizeof(int));// iarray(0,Nm0-1);
		(*allslipmodels).disc=darray(0,Nm0-1);
		(*allslipmodels).Nfaults=iarray(0,Nm0-1);
		(*allslipmodels).no_slipmodels=iarray(0,Nm0-1);
		(*allslipmodels).is_used=iarray(0,Nm0-1);
		(*allslipmodels).slipmodels = malloc(Nm0*sizeof(char*));
		(*allslipmodels).tsnap= NULL;	// will allocate if needed later on.
		nsm=0;
		if(procId == 0) {
			for (int nn=0; nn<Nm0; nn++) {
				(*allslipmodels).disc[nn] = res;
				fgets(line,Nchar,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				if (is_aseismic){
					sscanf(line,"%s %d", time_str, &no_slipmod);	//NB: no_slipmod changes at each nn iteration.

					//check if the number of snapshots is in agreement with flags:
					if (aseismic_splines==1 & no_slipmod<2){
						print_screen("Error: more than 1 snapshot required if time step mode set to splines (file %s).\n", input_fname);
						print_logfile("Error: more than 1 snapshot required if time step mode set to splines (file %s).\n", input_fname);
						fileError=1;
					}

					if (aseismic_log==1 & no_slipmod!=1){
						print_screen("Error: multiple snapshots not allowed if time step mode set to log (file %s).\n", input_fname);
						print_logfile("Error: multiple snapshots not allowed if time step mode set to log (file %s).\n", input_fname);
						fileError=1;
					}

					//check if no. of snapshots (single or multiple) is always the same:
					if ((*flag_multisnap==1 & no_slipmod<2) | (*flag_multisnap==0 & no_slipmod!=1)){
						print_screen("Error: aseismic events should either have a single snapshot each, or multiple snapshots each (file %s).\n", input_fname);
						print_logfile("Error: aseismic events should either have a single snapshot each, or multiple snapshots each (file %s).\n", input_fname);
						fileError=1;
					}

					//update value of flag_multisnap:
					*flag_multisnap= (no_slipmod>1) ? 1 : 0;
				}	

				else{
                    sscanf(line,"%s %lf %d", time_str, (*allslipmodels).mmain+nn, &no_slipmod);     //NB: no_slipmod changes at each nn iteration.
				}						
				sscanf(time_str, "%d-%d-%dT%d:%d:%dZ", &(times.tm_year), &(times.tm_mon), &(times.tm_mday), &(times.tm_hour), &(times.tm_min), &(times.tm_sec));
				times.tm_year-=1900;
				times.tm_mon-=1;
				times.tm_isdst=0;
				(*allslipmodels).tmain[nn]=difftime(mktime(&times),mktime(&reftime))*SEC2DAY;

				//check if catalog is chronological:
 			    if (nn>=1 && (*allslipmodels).tmain[nn]<(*allslipmodels).tmain[nn-1]){
					print_logfile("Error: slip model list in file %s no chronological. Exiting.\n", input_fname);
					print_screen("Error: slip model list in file %s not chronological. Exiting.\n", input_fname);
					fileError=1;
				}


				 (*allslipmodels).no_slipmodels[nn]=no_slipmod;
				 if (is_aseismic) (*allslipmodels).tsnap= (double *) realloc((*allslipmodels).tsnap, (nsm+1+no_slipmod) * sizeof(double));
				 if (nsm+1+no_slipmod>Nm0) {
					 (*allslipmodels).slipmodels=realloc((*allslipmodels).slipmodels, (nsm+1+no_slipmod) * sizeof(char*));
					(*allslipmodels).cut_surf=realloc((*allslipmodels).cut_surf, (nsm+1+no_slipmod) * sizeof(int)); 
					size_slipmodels = nsm+1+no_slipmod; // [Fahad] Used to bcast the final size
				 }

				 for (int n=1; n<=no_slipmod; n++){
					(*allslipmodels).slipmodels[nsm] = malloc(120 * sizeof(char));
					fgets(line,Nchar,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
					if (is_aseismic){
//						(*allslipmodels).Nfaults[nsm]=1; //actual value found later.
						sscanf(line,"%lf %s", (*allslipmodels).tsnap+nsm, (*allslipmodels).slipmodels[nsm]);
						(*allslipmodels).tsnap[nsm]+=(*allslipmodels).tmain[nn];
						(*allslipmodels).cut_surf[nsm]=0;
					}
			
					else{
						sscanf(line,"%d %s", (*allslipmodels).cut_surf+nsm, (*allslipmodels).slipmodels[nsm]);
					}
					nsm++;

				}
				
			}
			fclose(fin);
		}


		#ifdef _CRS_MPI
				MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
				if (is_aseismic){
					MPI_Bcast(flag_multisnap, 1, MPI_INT, 0, MPI_COMM_WORLD);
				}
		#endif  

		if (fileError){
			return(1);
		}

		#ifdef _CRS_MPI
			MPI_Bcast((*allslipmodels).tmain, Nm0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast((*allslipmodels).cut_surf, Nm0, MPI_INT, 0, MPI_COMM_WORLD);

			nsm = 0;

			if (!is_aseismic) MPI_Bcast((*allslipmodels).mmain, Nm0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast((*allslipmodels).no_slipmodels, Nm0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&size_slipmodels, 1, MPI_INT, 0, MPI_COMM_WORLD);

			if(procId != 0) {

				if (is_aseismic) (*allslipmodels).tsnap= (double *) malloc(size_slipmodels * sizeof(double));

				// If root did reallocation
				if(size_slipmodels > Nm0) {
					(*allslipmodels).slipmodels = realloc((*allslipmodels).slipmodels, size_slipmodels * sizeof(char*));
				}
				for(int nn = 0; nn < Nm0; ++nn) {
					no_slipmod=(*allslipmodels).no_slipmodels[nn];
					(*allslipmodels).disc[nn] = res;

					for(int n = 1; n <= no_slipmod; ++n) {
						if (is_aseismic){
							(*allslipmodels).Nfaults[nsm] = 1;       //actual value found later.
						}
						(*allslipmodels).slipmodels[nsm] = malloc(120 * sizeof(char));
						nsm++;
					}
				}
			}

			nsm = 0;

			if (is_aseismic) MPI_Bcast((*allslipmodels).tsnap, size_slipmodels, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			for(int nn = 0; nn < Nm0; ++nn) {
				no_slipmod=(*allslipmodels).no_slipmodels[nn]; 
				for(int n = 1; n <= no_slipmod; ++n) {
					MPI_Bcast((*allslipmodels).slipmodels[nsm], 120, MPI_CHAR, 0, MPI_COMM_WORLD);
					nsm++;
				}
			}
		#endif
	}

	nsm=0;
	if (is_aseismic) print_logfile("\nAfterslip input file: %s.\n", input_fname);
	else print_logfile("\nSlip input file: %s.\n", input_fname);
	print_logfile("%d %s events:\n", (*allslipmodels).NSM, is_aseismic? "aseismic" : "seismic");
	for (int m=0; m<(*allslipmodels).NSM; m++){
		if (is_aseismic){
			if (m==0) print_logfile("\t time \t name\n");
			print_logfile("\t%.2lf\t%s\n", (*allslipmodels).tmain[m], (*allslipmodels).slipmodels[m]);
		}
		else{
			if (m==0) print_logfile("\t time \t mag \t name\n");
			for (int n=1; n<=(*allslipmodels).no_slipmodels[m]; n++){
				print_logfile("\t%.2lf\t%.2lf\t%s\n", (*allslipmodels).tmain[m], (*allslipmodels).mmain[m], (*allslipmodels).slipmodels[nsm]);
				nsm++;
			}
		}
	}

	return 0;
}

