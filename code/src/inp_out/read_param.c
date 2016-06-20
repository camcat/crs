
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


#include <stdio.h>
#include <time.h>
#include "../defines.h"
//#include "../general/mem_mgmt.h"
#include "../okada/prestress.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int read_modelparameters(char *modelparametersfile, struct crust *crst, struct tm reftime,
						int *fixr, int *fixAsig, int *fixta, double *r0, double *Asig0,
						double *ta0, int *asig_log_step, int *ta_log_step, double *Asig_min, double *Asig_max, double *ta_min,
						double *ta_max, int *nAsig0, int *nta0, double *tw, double *fore_dt,
						int *Nsur, struct flags *flags,
						double *Mc, double *Mag_main, double *Mc_source,
						double *dCFS, double *DCFS_cap, double *dt, double *dM,
						double *xytoll, double *ztoll, double *border, double *res,
						double *gridresxy, double *gridresz, double *smoothing,
						int *LLinversion, int *forecast) {
	/*
	* Reads a parameter file and calcualtes regional stress field.
	*
	* Input:
	*  modelparametersfile: file name.
	*
	* Below is an example model parameter file.
	*
	*
	//#=============================================#
	//#       Coulomb stress parameters             #
	//#=============================================#
	//#
	//# 1. Msource, extra_dist: Minimum magnitude of events to be used as stress sources; extra distance oudside model domain for which sources should be included (both horizontal/vertical distance, km).
	//# 2. source_mode_focmec (iso/fm), source_mode_nofocmec (no/iso/fix):
	//#    source_mode_focmed controls how sources with a known focal mechanisms are treated: as isotropic sources (iso) or with synthetic slip models from the focal mechanisms (fm).
	//#    source_mode_nofocmec controls out a known focal mechanism are treated: ignored (no), as isotropic sources (iso), with synthetic slip models from the uniform regional mechanism (fix).
	//#	note: if source_mode_focmec=iso, source_mode_nofocmec will also automatically be set to "iso" (since there is no distinction between the two event types).
	//# 3. fm_res: Resolution (patch length, km) of the synthetic slip models created from focal mechanisms. If the resolution is small enough that a multiple-patch model is created, this model will be tapered on each side (instead of simply uniform slip).
	//# 4. DCFS_min, DCFS_cap:
	//#    DCFS_min is the smallest Coulomb stress change which should be computed (used to select a subset of grid points for each source and reduce calculations).
	//#    DCFS_cap: values of Coulomb stress change with absolute exceeding DCFS_cap will be capped to DCFS_cap (to avoid singularities in Okada solutions).
	//# 5. res_xy, res_z: internal horizontal/vertical model resolution.
	//5.95	0.0
	//fm	iso
	//3000
	//80      1e6
	//5000.0	1000.0
	//#=============================================#
	//#          Parameter Inversion	              #
	//#=============================================#
	//# Section controling the inversion of Rate-and-State parameters (Asigma, ta, r0).
	//# The inversion is based on maximizing the LogLikelihood, using a simple a grid search algorithm.
	//# For details, see Hainzl et al. (2009), Appendix A1.
	//#
	//# 1. LLinversion: Flag indicating if inversion for RS parameters should be performed.
	//#    If LLinversion is set to 0, the following three lines should indicate the default values of rate and state parameters (and fixX=0 will be ignored).
	//# 2. fixr0, [r0]: fixr is a flag indicating if r0 should be inverted for; if set to 0, it should be followed by the default value of r0.
	//#    r0 is the daily rate of earthquakes in the domain; if a non-homogeneous background rate is used, the rate in each cell is normalized to give total rate r0.
	//# 3. This line has two possible forms (without commas):
	//#	fixAsig, Asig0:	is fixAsig=1 or LLinversion=0, Asig0 is the default value [this option should be used if LLinversion=0]. Unit=Pa.
	//#	fixAsig, Asig_min, Asig_max, nAsig, Asig_stepmode: if fixAsig=0, Asig_(min/max) the lower/upper bound of Asig, and nAsig the number of values to be tested.
	//#    	Asig_step_mode is one of: lin/log, and determines whether the steps between Asig_(min/max) will be linear or logarithmic.
	//# 4. This line has two possible forms:
	//#       fixta, ta0: is fixta=1 or LLinversion=0, ta0 is the default value [this option should be used if LLinversion=0]. Unit=days.
	//#       fixta, ta_min, ta_max, nta, ta_stepmode: if fixta=0, the ta_(min/max) are the lower/upper bound of ta, and the nta number of values to be tested.
	//#    ta_step_mode is one of: lin/log, and determines whether the steps between Asig_(min/max) will be linear or logarithmic.
	//#	fixAsig, Asig_min, Asig_max, nAsig: if fixAsig=1, Asig_(min/max) the lower/upper bound of Asig, and nAsig the number of values to be tested.
	//# 5. Mc: min. magnitude of events to be used for LogLikelihood calculations.
	//#       if Mc>20, the program will estimate the catalog completeness magnitude.
	//# 6. Controls the time windows to be excluded from the LogLikelihood calculation due to incomplete catalog.
	//#    tw, Magmain: tw=length of time window to skip (days); Magmain= min. magnitude of mainshocks following which a time window should be skipped;
	//#    If no time should be skipped, set Mag_main to a large value (setting Mag_main to a small value and tw=0.0 leads to longer computation times).
	//1
	//1      0.022
	//0      8000   10000   1	lin
	//0      7000   10000   1	lin
	//2.0
	//0.0012  5.95
	//#=============================================#
	//#       Treatment of Uncertainties            #
	//#=============================================#
	//# This section controls the treatment of uncertainties in the Coulomb stress field.
	//# For a description of each type uncertainty is treated, see Cattania et al. (2014).
	//# 1. Receiver fault orientation.
	//#    options: oops (optimally oriented planes), fixed (uses the uniform reginal field plane given below), focmec (uses the planes from the focal mechanisms catalog, performing Monte Carlo iterations over the available mechanisms).
	//# 2. grid_error (0/1): flag indicating if the error due to finite grid size should be calculated.
	//# 3. Nsur: number of Monte Carlo iterations.
	//fixed
	//0
	//9
	//#=============================================#
	//#             Other parameters                #
	//#=============================================#
	//#
	//# 1. output_forecast, dt: flag indicating if forecast should be produced; time step used for output of seismicity temporal evolution.
	//# 2-5. Expected difference between the ZMAP earthquake catalog (InputCatalogFile) and the catalog of focal mechanisms (InputCatalogFocMecFile): these values are used as tolerance when associating focal mechanisms to the ZMAP events.
	//# 2. dt: time tolerance (days)
	//# 3. dMag: magnitude tolerance
	//# 4. dxy: horizontal distance tolerance (km)
	//# 5. dz: vertical distance tolerance (km)
	//# 6. smoothing: smoothing distance used for calculating background rate from a catalog (if parameter InputBackgroundRateCatalog is used).
	//1	0.1
	//0.0005
	//0.9
	//50
	//0
	//5.0
	//#=====================================================#
	//# 	Parameters Describing crustal properties	    #
	//#=====================================================#
	//# Elastic paremeters: lambda, mu (MPa)
	//# Friction coefficient, Skeption coefficient
	//# Uniform regional field description:	strike, dip, rake
	//# Choice between "paxis" and "oop", two ways in which the regional stress field can be described:
	//#	if "paxis", the following 3 lines give amplitude (MPa), strike, dip of the principal axis:
	//#		sigma1	strike1	dip1
	//#               sigma2  strike2 dip2
	//#               sigma3  strike3 dip3
	//#	if "oop", the following lines give the amplitudes of the principal axis and the strike, dip and rake
	//#	of a master mechanism, assumed to be optimally oriented according to the Coulomb failure criterion.
	//#		sigma1 sigma2 sigma
	//#
	//31226	26624
	//0.3	0.0
	//330.000  89.000  180.000
	//paxis
	//-5.0	6.646198	-0.596911
	//5.0	96.654556	-0.802279
	//0.0	60.000014	89.000049
	*/



	// Variables used for MPI.

	int procId = 0;
	int fileError = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	FILE * fin;
	char comment[]="#", comm=comment[0];
	int Nchar_long=500;
	char line[Nchar_long];
	int optrake;
	char sourcemode_fm[10], sourcemode_nofm[10], logstep[10];
	char recfault[10];
	struct tm;
	char regstress_mode[120];
	double s[3];	//regional stress field description;
	double st[3];	//regional stress field description;
	double di[3];	//regional stress field description;

	// If there is a file error, only root will know about it.
	// So it is important that the error is broadcast to all
	// processes, and they all return from the function with
	// the same error code.
	if(procId == 0) {
		fin = fopen(modelparametersfile,"r");
		if(fin == NULL) {
			print_screen("Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);
			print_logfile("Error: parameter file %s could not be opened. Exit. \n", modelparametersfile);

			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		return 1;
	}

	//initialize crust structure:
	init_crst(crst);

	if(procId == 0) {
		sprintf(comment,"#");

		//-------------Coulomb stress parameters------------------//
		comm=comment[0];
		line[0]=comm;
		while (line[0]==comm){
			fgets(line,Nchar_long,fin);
		}
		sscanf(line,"%lf %lf", Mc_source, border);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%s %s", sourcemode_fm, sourcemode_nofm);

		if (!strcmp(sourcemode_fm, "iso")){
			(*flags).sources_all_iso=1;
			(*flags).sources_without_focmec=1;
		}
		else {
			if (!strcmp(sourcemode_fm, "fm")){
				(*flags).sources_all_iso=0;
				if (!strcmp(sourcemode_nofm, "no")){
					(*flags).sources_without_focmec=0;
				}
				else{
					if (!strcmp(sourcemode_nofm, "iso")){
						(*flags).sources_without_focmec=1;
					}
					else{
						if (!strcmp(sourcemode_nofm, "fix")){
							(*flags).sources_without_focmec=2;
						}
						else{
							print_screen("Illegal value for source_mode_nofocmec (should be one of: no/iso/fix).\n", modelparametersfile);
							print_logfile("Illegal value for source_mode_nofocmec (should be one of: no/iso/fix).\n", modelparametersfile);
							fileError=1;
						}
					}
				}
			}
			else {
				print_screen("Illegal value for source_mode_focmec (should be one of: iso/fm).\n", modelparametersfile);
				print_logfile("Illegal value for source_mode_focmec (should be one of: iso/fm).\n", modelparametersfile);
				fileError=1;
			}
		}

		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", res);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", dCFS, DCFS_cap);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", gridresxy, gridresz);

		//-------------Parameter inversion------------------//

		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%d", LLinversion);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf", fixr, r0);
		//the following 2 lines (Asig, ta values) can have 2 alternative forms:
		//fixAsig, Asig0: is fixAsig=1 or LLinversion=0, Asig0 is the default value [this option should be used if LLinversion=0]. Unit=Pa.
		//fixAsig, Asig_min, Asig_max, nAsig, Asig_stepmode: if fixAsig=0, Asig_(min/max) the lower/upper bound of Asig, and nAsig the number of values to be tested.//1 Asig0
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf %lf  %d %s", fixAsig, Asig_min, Asig_max, nAsig0, logstep);
		if (*fixAsig || !*LLinversion){
			*Asig0=*Asig_min;
			*Asig_min=*Asig_max=0.0;
			*nAsig0=0;
			*asig_log_step=0;
		}
		else{
			*Asig0=0.0;
			if (!strcmp(logstep, "log")) {
				*asig_log_step=1;
			}
			else if (!strcmp(logstep, "lin")) {
				*asig_log_step=0;
			}
			else{
				print_screen("Illegal value for Asigma step mode (should be one of: lin/log).\n", modelparametersfile);
				print_logfile("Illegal value for Asigma step mode (should be one of: lin/log.\n", modelparametersfile);
				fileError=1;
			}
		}
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d %lf %lf %d %s", fixta, ta_min, ta_max, nta0, logstep);
		if (*fixta || !*LLinversion){
			*ta0=*ta_min;
			*ta_min=*ta_max=0.0;
			*nta0=0;
			*ta_log_step=0;
		}
		else{
			*ta0=0.0;
			if (!strcmp(logstep, "log")) {
				*ta_log_step=1;
			}
			else if (!strcmp(logstep, "lin")) {
				*ta_log_step=0;
			}
			else{
				print_screen("Illegal value for ta step mode (should be one of: lin/log).\n", modelparametersfile);
				print_logfile("Illegal value for ta step mode (should be one of: lin/log.\n", modelparametersfile);
				fileError=1;
			}
		}
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", Mc);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", tw, Mag_main);

		//-------------Treatment of uncertainties------------------//
		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%s %d", recfault, &optrake);
		(*flags).optrake=0;	//default, in case the user forgets to define it.
		if (!strcmp(recfault, "fixed")){
			(*flags).err_recfault=0;
			(*flags).OOPs=0;
			(*flags).optrake=optrake;
		}
		else{
			if (!strcmp(recfault,"oops")){
				(*flags).err_recfault=0;
				(*flags).OOPs=1;
				(*flags).optrake=1;
			}
			else {
				if (!strcmp(recfault,"focmec")){
					(*flags).err_recfault=1;
					(*flags).OOPs=0;
					(*flags).optrake=optrake;
				}
				else{
					print_screen("Illegal value for choice of receiver fault in file %s (should be one of: oops/fixed/focmec).\n", modelparametersfile);
					print_logfile("Illegal value for choice of receiver fault in file %s (should be one of: oops/fixed/focmec).\n", modelparametersfile);
					fileError=1;
				}
			}

		}
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", &((*flags).err_gridpoints));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%d", Nsur);
		*Nsur=MAX(*Nsur,1);


		//-------------Other parameters------------------//

		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%d %lf", forecast, fore_dt);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", dt);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", dM);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", xytoll);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", ztoll);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf", smoothing);


		//-------------Parameters Describing crustal properties------------------//

		line[0]=comm;
		while (line[0]==comm) fgets(line,Nchar_long,fin);
		if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fscanf!\n");
		sscanf(line,"%lf %lf", &((*crst).lambda), &((*crst).mu));
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf", &((*crst).fric), &((*crst).skempton));
		(*crst).fric=(*crst).fric*(1.0-(*crst).skempton);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%lf %lf %lf", (*crst).str0, (*crst).dip0, (*crst).rake0);
		fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
		sscanf(line,"%s", regstress_mode);

		if (!(strcmp(regstress_mode,"oops"))) {
			fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
			sscanf(line,"%lf %lf %lf", s, s+1, s+2);
		}
		else {
			if (!(strcmp(regstress_mode,"paxis"))) {
				fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%lf %lf %lf", s, st, di);
				fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%lf %lf %lf", s+1, st+1, di+1);
				fgets(line,Nchar_long,fin); if (ferror(fin)) fprintf(stderr, "ERROR reading input data using fgets!\n");
				sscanf(line,"%lf %lf %lf", s+2, st+2, di+2);
			}

			else {
				print_logfile("Invalid value for mode of regional stress field: should be 'oops' or 'paxis'. Exit.\n");
				print_screen("Invalid value for mode of regional stress field: should be 'oops' or 'paxis'. Exit.\n");
				fileError=1;
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

		// Copy scalars to the BCast_Model_Parameters struct
		struct BCast_Model_Parameters modelParams;

		if(procId == 0) {
			//Integers:
			modelParams.fixr         =  *fixr;			//1
			modelParams.fixAsig      =  *fixAsig;		//2
			modelParams.fixta        =  *fixta;			//3
			modelParams.nAsig0       =  *nAsig0;		//4
			modelParams.nta0         =  *nta0;			//5
			modelParams.Nsur         =  *Nsur;			//6
			modelParams.LLinversion  =  *LLinversion;	//7
			modelParams.forecast     =  *forecast;		//8
			modelParams.ta_log_step  =  *ta_log_step;	//9
			modelParams.asig_log_step=  *asig_log_step;	//10

			//doubles:
			modelParams.r0  		 =  *r0;			//11
			modelParams.Asig0        =  *Asig0;			//12
			modelParams.ta0  		 =  *ta0;			//13
			modelParams.Asig_min     =  *Asig_min;		//14
			modelParams.Asig_max     =  *Asig_max;		//15
			modelParams.ta_min       =  *ta_min;		//16
			modelParams.ta_max       =  *ta_max;		//17
			modelParams.tw   		 =  *tw;			//18
			modelParams.fore_dt      =  *fore_dt;	//19
			modelParams.Mc   		 =  *Mc;		//20
			modelParams.Mag_main     =  *Mag_main;	//21
			modelParams.Mc_source    =  *Mc_source;	//22
			modelParams.dCFS         =  *dCFS;		//23
			modelParams.DCFS_cap     =  *DCFS_cap;	//24
			modelParams.dt 			 =  *dt;		//25
			modelParams.dM  		 =  *dM;		//26
			modelParams.xytoll       =  *xytoll;	//27
			modelParams.ztoll        =  *ztoll;		//28
			modelParams.border       =  *border;	//29
			modelParams.res 		 =  *res;		//30
			modelParams.gridresxy    =  *gridresxy;	//31
			modelParams.gridresz     =  *gridresz;	//32
			modelParams.smoothing    =  *smoothing;	//33

		}

		//Map the BCast_Model_Parameters struct to MPI struct type
		MPI_Datatype CRS_MPI_BCast_Model_Parameters;
		int blocks_ModelParameters				[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Datatype types_ModelParameters		[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Aint addresses_ModelParameters		[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Aint displacements_ModelParameters	[SIZE_BCAST_MODEL_PARAMETERS];
		MPI_Aint baseAddress_ModelParameters;

		// Set blocks
		for(int i = 0; i < SIZE_BCAST_MODEL_PARAMETERS; ++i) {
			blocks_ModelParameters[i] = 1;
		}

		// Set types
		for(int i = 0; i < 10; ++i) {
			types_ModelParameters[i] = MPI_INT;
		}
		for(int i = 10; i < SIZE_BCAST_MODEL_PARAMETERS; ++i) {
			types_ModelParameters[i] = MPI_DOUBLE;
		}

		// Set addresses
		MPI_Address(&modelParams, &baseAddress_ModelParameters);

		MPI_Address(&(modelParams.fixr),        &addresses_ModelParameters[0]);
		MPI_Address(&(modelParams.fixAsig),     &addresses_ModelParameters[1]);
		MPI_Address(&(modelParams.fixta),       &addresses_ModelParameters[2]);
		MPI_Address(&(modelParams.nAsig0),      &addresses_ModelParameters[3]);
		MPI_Address(&(modelParams.nta0),        &addresses_ModelParameters[4]);
		MPI_Address(&(modelParams.Nsur),        &addresses_ModelParameters[5]);
		MPI_Address(&(modelParams.LLinversion), &addresses_ModelParameters[6]);
		MPI_Address(&(modelParams.forecast),    &addresses_ModelParameters[7]);
		MPI_Address(&(modelParams.ta_log_step), &addresses_ModelParameters[8]);
		MPI_Address(&(modelParams.asig_log_step),&addresses_ModelParameters[9]);

		MPI_Address(&(modelParams.r0), 			&addresses_ModelParameters[10]);
		MPI_Address(&(modelParams.Asig0),       &addresses_ModelParameters[11]);
		MPI_Address(&(modelParams.ta0), 		&addresses_ModelParameters[12]);
		MPI_Address(&(modelParams.Asig_min),    &addresses_ModelParameters[13]);
		MPI_Address(&(modelParams.Asig_max),    &addresses_ModelParameters[14]);
		MPI_Address(&(modelParams.ta_min),      &addresses_ModelParameters[15]);
		MPI_Address(&(modelParams.ta_max),      &addresses_ModelParameters[16]);
		MPI_Address(&(modelParams.tw),  		&addresses_ModelParameters[17]);
		MPI_Address(&(modelParams.fore_dt),     &addresses_ModelParameters[18]);
		MPI_Address(&(modelParams.Mc),  		&addresses_ModelParameters[19]);
		MPI_Address(&(modelParams.Mag_main),    &addresses_ModelParameters[20]);
		MPI_Address(&(modelParams.Mc_source),   &addresses_ModelParameters[21]);
		MPI_Address(&(modelParams.dCFS),        &addresses_ModelParameters[22]);
		MPI_Address(&(modelParams.DCFS_cap),    &addresses_ModelParameters[23]);
		MPI_Address(&(modelParams.dt),  		&addresses_ModelParameters[24]);
		MPI_Address(&(modelParams.dM),  		&addresses_ModelParameters[25]);
		MPI_Address(&(modelParams.xytoll),      &addresses_ModelParameters[26]);
		MPI_Address(&(modelParams.ztoll),       &addresses_ModelParameters[27]);
		MPI_Address(&(modelParams.border),      &addresses_ModelParameters[28]);
		MPI_Address(&(modelParams.res), 		&addresses_ModelParameters[29]);
		MPI_Address(&(modelParams.gridresxy),   &addresses_ModelParameters[30]);
		MPI_Address(&(modelParams.gridresz),    &addresses_ModelParameters[31]);
		MPI_Address(&(modelParams.smoothing),   &addresses_ModelParameters[32]);

		// Set displacements
		for(int i = 0; i < SIZE_BCAST_MODEL_PARAMETERS; ++i) {
			displacements_ModelParameters[i] = addresses_ModelParameters[i] - baseAddress_ModelParameters;
		}

		MPI_Type_struct(SIZE_BCAST_MODEL_PARAMETERS, blocks_ModelParameters,
						displacements_ModelParameters, types_ModelParameters,
						&CRS_MPI_BCast_Model_Parameters);

		MPI_Type_commit(&CRS_MPI_BCast_Model_Parameters);

		MPI_Bcast(&modelParams, 1, CRS_MPI_BCast_Model_Parameters, 0, MPI_COMM_WORLD);

		// For processes other than root, populate the scalars from the struct
		// received from root.
		if(procId != 0) {
			*fixr   		= modelParams.fixr;
			*fixAsig    	= modelParams.fixAsig;
			*fixta  		= modelParams.fixta;
			*nAsig0     	= modelParams.nAsig0;
			*nta0   		= modelParams.nta0;
			*Nsur   		= modelParams.Nsur;
			*LLinversion    = modelParams.LLinversion;
			*forecast       = modelParams.forecast;
			*r0     		= modelParams.r0;
			*Asig0  		= modelParams.Asig0;
			*ta0    		= modelParams.ta0;
			*ta_log_step	= modelParams.ta_log_step;
			*asig_log_step  = modelParams.asig_log_step;
			*Asig_min       = modelParams.Asig_min;
			*Asig_max       = modelParams.Asig_max;
			*ta_min         = modelParams.ta_min;
			*ta_max         = modelParams.ta_max;
			*tw     		= modelParams.tw;
			*fore_dt        = modelParams.fore_dt;
			*Mc     		= modelParams.Mc;
			*Mag_main       = modelParams.Mag_main;
			*Mc_source      = modelParams.Mc_source;
			*dCFS   		= modelParams.dCFS;
			*DCFS_cap       = modelParams.DCFS_cap;
			*dt     		= modelParams.dt;
			*dM     		= modelParams.dM;
			*xytoll         = modelParams.xytoll;
			*ztoll  		= modelParams.ztoll;
			*border         = modelParams.border;
			*res    		= modelParams.res;
			*gridresxy      = modelParams.gridresxy;
			*gridresz       = modelParams.gridresz;
			*smoothing      = modelParams.smoothing;
		}


		// Map struct flags to MPI struct type
		MPI_Datatype CRS_MPI_BCast_Flags;
		int blocks_Flags[BCAST_FLAGS_SIZE];
		MPI_Datatype types_Flags[BCAST_FLAGS_SIZE];
		MPI_Aint baseAddress_Flags;
		MPI_Aint addresses_Flags[BCAST_FLAGS_SIZE];
		MPI_Aint displacements_Flags[BCAST_FLAGS_SIZE];

		// Set blocks and types
		for(int i = 0; i < BCAST_FLAGS_SIZE; ++i) {
			blocks_Flags[i] = 1;
			types_Flags[i] = MPI_INT;
		}

		// Set addresses
		MPI_Address(flags, &baseAddress_Flags);

		MPI_Address(&(flags->err_recfault),      		&addresses_Flags[0]);
		MPI_Address(&(flags->err_gridpoints),    		&addresses_Flags[1]);
		MPI_Address(&(flags->OOPs),      		 		&addresses_Flags[2]);
		MPI_Address(&(flags->aseismic),         		&addresses_Flags[3]);
		MPI_Address(&(flags->aseismic_linear), 		 	&addresses_Flags[4]);
		MPI_Address(&(flags->aseismic_multisnap),   	&addresses_Flags[5]);
		MPI_Address(&(flags->sources_all_iso),   		&addresses_Flags[6]);
		MPI_Address(&(flags->sources_without_focmec),   &addresses_Flags[7]);
		MPI_Address(&(flags->sample_all),        		&addresses_Flags[8]);
		MPI_Address(&(flags->optrake),        			&addresses_Flags[9]);


		// Set displacements
		for(int i = 0; i < BCAST_FLAGS_SIZE; ++i) {
			displacements_Flags[i] = addresses_Flags[i] - baseAddress_Flags;
		}

		MPI_Type_struct(BCAST_FLAGS_SIZE, blocks_Flags, displacements_Flags,
						types_Flags, &CRS_MPI_BCast_Flags);
		MPI_Type_commit(&CRS_MPI_BCast_Flags);

		MPI_Bcast(flags, 1, CRS_MPI_BCast_Flags, 0, MPI_COMM_WORLD);

		// Broadcast variables from the crst struct.
		MPI_Bcast(s,  3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(st, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(di, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Bcast(&(crst->lambda),	 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->mu), 		 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->fric), 	 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->skempton),  1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->str0[0]),  1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->dip0[0]),  1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(crst->rake0[0]), 1,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(regstress_mode, 	 120, MPI_CHAR,   0, MPI_COMM_WORLD);

	#endif

	//calculate regional stress field:
	if (!(strcmp(regstress_mode,"oops"))) {
		prestress(1e6*s[0], 1e6*s[1], 1e6*s[2], (*crst).str0[0], (*crst).dip0[0], (*crst).rake0[0], 0.0,(*crst).fric, &((*crst).S));
	}
	else{
		if (!(strcmp(regstress_mode,"paxis"))) {
			for (int i=0; i<3; i++) s[i]*=-1e6;
			(*crst).S=prestress_eigen(s, st, di);
		}
	}

	return 0;
}
