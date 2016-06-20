
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


#include "okadaDCFS.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif


//---------------------------------------------------------------------//
//-----					Top level functions						  -----//
//---------------------------------------------------------------------//

int resolve_DCFS0(struct pscmp DCFS, struct crust crst, double *strikeRs, double *dipRs, double *rakeRs, int optrake, int cap){
	/* Resolves stress tensor in DCFS on receiver faults given by strikeRs, dipRs, rake, and stores result in DCFS.cmb.
	 * To be called after filling in DCFS.S. (stress tensor).
	 *
	 * Input:
	 *
	 *  DCFS.S: stress tensor.
	 *  strikeRs, dipRs:	strike and dip. It can be a vector if multiple zones with different receiver faults are used
	 *  crst.nofmzones:	contains no. of receiver fault zones, i.e. size of strikeRs,dipRs [0...crst.nofmzones-1]
	 *  rake: pointer to value of rake (a single value, not a vector). if NULL, optimal rake will be used.
	 *  optrake:	flag indicating if optimal rake should be used. If optrake=1, the value of *rake will be ignored.
	 *  cap: flag determining whether coulomb values should be capped to DCFS_cap.
	 * Output:
	 *
	 *  DCFS.cmb: coulomb stress field.
	 *
	 */

	double *sigma0s;
	double **stress0s;
	double **n, **s;	//normal vectors, slip vectors.
	double MaxDCFS=DCFS_cap;
	int Nsel=DCFS.nsel;
	int fm;
	double strikeR, dipR, rakeR;
	int no_fm_zones=crst.nofmzones;

	n=malloc((no_fm_zones)*sizeof(double *));
	s=malloc((no_fm_zones)*sizeof(double *));
	stress0s=malloc((no_fm_zones)*sizeof(double *));
	sigma0s=darray(0,no_fm_zones-1);

	for (int i=0; i<no_fm_zones; i++) {
		//calculate normal vector:
		strikeR=strikeRs[i];
			dipR=dipRs[i];
		if (!optrake && rakeRs) rakeR=rakeRs[i];
		n[i]=normal_vector(strikeR, dipR);

		if (!optrake){
			if (!rakeRs) {
				print_screen("** Warning: optrake=0, but rake is NULL: will use optimal rake (resolve_DCFS).**\n");
				print_logfile("** Warning: optrake=0, but rake is NULL: will use optimal rake (resolve_DCFS).**\n");
				optrake=1;
			}
			else{
			//calculate slip vector:
				s[i]=slip_vector(strikeR, dipR, rakeR);
			}
		}
		else s[i]=NULL;

		//calculated background stress and normal stress:
		stress0s[i]=mtimesv(crst.S,n[i],NULL,3,3);
		sigma0s[i]=-1.0*vdotv(stress0s[i],n[i],3);		// compression is positive (rock mechanics convention).

	}

	#pragma omp parallel for private(fm)
	for (int i=1; i<=Nsel; i++){
		fm= (no_fm_zones==1) ? 0 : crst.fmzone[i];	//zone index.
		DCFS.cmb[i]=resolve_n(DCFS.S[i], n[fm], crst.fric, stress0s[fm], sigma0s[fm], s[fm]);
		if (cap){
		if (DCFS.cmb[i]>MaxDCFS) DCFS.cmb[i]=MaxDCFS;
		if (DCFS.cmb[i]<-MaxDCFS) DCFS.cmb[i]=-MaxDCFS;
		}
	}

	for (int i=0; i<no_fm_zones; i++){
		if (s[i]) free_darray(s[i],1,3);
		free_darray(n[i],1,3);
		free_darray(stress0s[i],1,3);
	}
	free(s);
	free(n);
	free(stress0s);
	free_darray(sigma0s,0,1);

	return(0);
}

#ifdef _CRS_MPI
int okadaCoeff_mpi(float ****Coeffs_st,
				   float ****Coeffs_dip,
				   float ****Coeffs_open,
				   struct eqkfm *eqkfm1,
				   int NF,
				   struct crust crst) {

	/* Allocates and populates arrays containing the okada coefficients between each fault and grid point.
	 * This function is where most of the memory used by the program in most cases is allocated, and it may still be optimized (see comments).
	 *
	 * Input:
	 *  eqkfm1: structure containing slip models. Range [0...NF-1].
	 *  crst: crust structure, containing position of grid points and elastic parameters.
	 */

	int procId = 0, numProcs = 1;
	int start, partitionSize;

	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	double north, east, eqnorth, eqeast;
	double len, width, depth; //for individual patches.
	double depth0; //to differentiate between blind fault (depth0=0) or fault cutting through surface.
	double strike, dip, rake;
	double Sxx, Syy, Szz, Sxy, Syz, Sxz;
	int Nsel = eqkfm1[0].nsel;
	int NP_tot=0, p1, i;
	int noslip_str, noslip_dip, noopen; //these are flags used to reduce no. of calculations (if there is no slip on a subfault).
	int err=0;
	int flag_open, flag_sslip, flag_dslip; //these are flags used to reduce memory allocation (if there is no slip anywhere).
	struct eqkfm *afslip;	//dummy variable equal to pointers in eqkfm1 which point to afterslip elements.

	for(int j=0; j<NF; j++) {
		NP_tot+=eqkfm1[j].np_di*eqkfm1[j].np_st;
	}

	depth0=eqkfm1[0].cuts_surf ? eqkfm1[0].top : 0.0;

	print_logfile("Depth of surface: %.3lf km.\n", depth0);
	print_screen("Depth of surface: %.3lf km.\n", depth0);

	//---------initialize DCFS----------//

	//check if Coeffs tensors should be allocated.
	flag_open=flag_sslip=flag_dslip=0;

	for (int j=0; j<NF; j++){

		//if the element is associated with afterslip, should check whether this afterslip has strike/slip/opening component.
		afslip=	eqkfm1[j].co_aft_pointer;

		//check if there is no slip in this fault:
		noslip_str= (eqkfm1[j].slip_str==NULL && (afslip==NULL || (*afslip).allslip_str==NULL));
		noslip_dip= (eqkfm1[j].slip_dip==NULL && (afslip==NULL || (*afslip).allslip_dip==NULL));
		noopen= (eqkfm1[j].open==NULL && (afslip==NULL || (*afslip).allslip_open==NULL));

		flag_open=max(flag_open, !noopen);
		flag_sslip=max(flag_sslip, !noslip_str);
		flag_dslip=max(flag_dslip, !noslip_dip);

		if (flag_open && flag_sslip && flag_dslip) break;
	}

	//todo could use array of pointers for smarter memory allocation (and avoid allocating if a subfault, or even a single patch, has no slip):
	//todo it may also be better to process each of the 3 tensors (Coeffs_XX) separately: less likely to run out of memory.

	//allocate memory.
	*Coeffs_st = (flag_sslip) ? f3array(1,NP_tot,1,Nsel,1,6) : NULL;
	*Coeffs_dip= (flag_dslip) ? f3array(1,NP_tot,1,Nsel,1,6) : NULL;
	*Coeffs_open= (flag_open) ? f3array(1,NP_tot,1,Nsel,1,6) : NULL;

	for (int i=1; i<=NP_tot; i++){
		for (int j=1; j<=Nsel; j++){
			for (int k=1; k<=6; k++){
				if (flag_sslip) (*Coeffs_st)[i][j][k]=0;
				if (flag_dslip) (*Coeffs_dip)[i][j][k]=0;
				if (flag_open) (*Coeffs_open)[i][j][k]=0;
			}
		}
	}

	//-----------------------------------------------------------------------------------------//
	//-----------Calculate Coulomb stress vector for each patch assuming slip=1.---------------//
	//-----------------------------------------------------------------------------------------//

	print_screen("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);
	print_logfile("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);

	for (int j=0; j<NF; j++) {

		// MPI 	-- 	Flag to indicate if the current
		//		--  fault should be processed in serial.g through surface.
		double strike, dip, rake;
		int processFaultSerially = 0;

		if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
			print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			return(1);
		}

		len   = eqkfm1[j].L*(1.0/eqkfm1[j].np_st);
		width = eqkfm1[j].W*(1.0/eqkfm1[j].np_di);

		// Since MPI parallelization is based on the No. of patches.
		size_t numPatches = eqkfm1[j].np_di*eqkfm1[j].np_st;

		// Create linearized tensors for all patches within the current fault,
		// for use in MPI communication routines.
		size_t fullTensorSize = ((numPatches) * Nsel * 6);
		float *coeffs_st  = (float*) malloc((size_t)(fullTensorSize * sizeof(float)));
		float *coeffs_dip = (float*) malloc((size_t)(fullTensorSize * sizeof(float)));
		float *coeffs_open = (float*) malloc((size_t)(fullTensorSize * sizeof(float)));

		// If the No. of MPI ranks is greater than the number of patches,
		// serially process all patches in the current fault.
		if(numProcs > numPatches) {
			processFaultSerially = 1;

			partitionSize = numPatches;

			start = 0;

			if(procId == 0) {
				printf("\n Number of processes: %d", numProcs);
				printf("\n Number of patches: %d", numPatches);
			}
			print_screen("*** No. of patches is less than the No. of processes. Processing fault in serial ... ***\n",j);
		}
		else {	// Partition the No. of patches for parallel processing by MPI ranks.
			partitionSize = numPatches / numProcs;

			// If partionSize is not large enough to hold all patches, increase
			// the partition size and reallocate the linearized tensors.
			if(partitionSize * numProcs != numPatches) {
				partitionSize += 1;
				coeffs_st  = (float*) realloc(coeffs_st,  (size_t)((partitionSize*Nsel*6*numProcs) * sizeof(float)));
				coeffs_dip = (float*) realloc(coeffs_dip, (size_t)((partitionSize*Nsel*6*numProcs) * sizeof(float)));
				coeffs_open = (float*) realloc(coeffs_open, (size_t)((partitionSize*Nsel*6*numProcs) * sizeof(float)));
			}

			start = (procId * partitionSize);
		}

		// Create linearized tensors for the partition to be processed
		// by the current rank. Linearization is required for MPI
		// communication routines.
		size_t partitionedTensorSize = partitionSize * Nsel * 6;
		float *coeffs_st_partitioned  = (float*) malloc((size_t)(partitionedTensorSize * sizeof(float)));
		float *coeffs_dip_partitioned = (float*) malloc((size_t)(partitionedTensorSize * sizeof(float)));
		float *coeffs_open_partitioned = (float*) malloc((size_t)(partitionedTensorSize * sizeof(float)));

		for(int i = 0; i < partitionedTensorSize; ++i) {
			coeffs_st_partitioned [i] = 0.0;
			coeffs_dip_partitioned[i] = 0.0;
			coeffs_open_partitioned[i] = 0.0;
		}

		int p2 = start, index=0;
		for(int p = 0; p < partitionSize; p++) {
			patch_pos(eqkfm1[j], p2+1, &eqeast, &eqnorth, &depth);
			++p2;

			#pragma omp parallel for private(Sxx, Syy, Szz, Sxy, Syz, Sxz, north, east, i, afslip, noslip_str, noslip_dip, noopen)
			for(int i0=0; i0<Nsel; i0++) {	

				i=eqkfm1[0].selpoints[i0+1];	// Added '+1' to the index
				north=crst.y[i];
				east=crst.x[i];

				//if the element is associated with afterslip, should check whether this afterslip has strike/slip/opening component.
				afslip=	eqkfm1[j].co_aft_pointer;

				//check if both coseismic and afterslip have empty arrays for each component of deformation:
				noslip_str= (eqkfm1[j].slip_str==NULL && (afslip==NULL || (*afslip).allslip_str==NULL));
				noslip_dip= (eqkfm1[j].slip_dip==NULL && (afslip==NULL || (*afslip).allslip_dip==NULL));
				noopen= (eqkfm1[j].open==NULL && (afslip==NULL || (*afslip).allslip_open==NULL));

				if (!noslip_str) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike,  dip, len, width, 1.0, 0.0, 0.0,
							north, east, crst.depth[i]-depth0, &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
							 crst.lambda, crst.mu);

					index = (p * Nsel * 6) + (i0 * 6);

					coeffs_st_partitioned[index + 0] += 1e6*Sxx;
					coeffs_st_partitioned[index + 1] += 1e6*Syy;
					coeffs_st_partitioned[index + 2] += 1e6*Szz;
					coeffs_st_partitioned[index + 3] += 1e6*Sxy;
					coeffs_st_partitioned[index + 4] += 1e6*Syz;
					coeffs_st_partitioned[index + 5] += 1e6*Sxz;
				}

				if (!noslip_dip) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike, dip, len, width, 0.0, -1.0, 0.0,
							 north, east, crst.depth[i]-depth0, &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
							 crst.lambda, crst.mu);

					index = (p * Nsel * 6) + (i0 * 6);

					coeffs_dip_partitioned[index + 0] += 1e6*Sxx;
					coeffs_dip_partitioned[index + 1] += 1e6*Syy;
					coeffs_dip_partitioned[index + 2] += 1e6*Szz;
					coeffs_dip_partitioned[index + 3] += 1e6*Sxy;
					coeffs_dip_partitioned[index + 4] += 1e6*Syz;
					coeffs_dip_partitioned[index + 5] += 1e6*Sxz;
				}

				if (!noopen) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike, dip, len, width, 0.0, 0.0, 1.0,
							 north, east, crst.depth[i]-depth0, &Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,
							 crst.lambda, crst.mu);

					index = (p * Nsel * 6) + (i0 * 6);

					coeffs_open_partitioned[index + 0] += 1e6*Sxx;
					coeffs_open_partitioned[index + 1] += 1e6*Syy;
					coeffs_open_partitioned[index + 2] += 1e6*Szz;
					coeffs_open_partitioned[index + 3] += 1e6*Sxy;
					coeffs_open_partitioned[index + 4] += 1e6*Syz;
					coeffs_open_partitioned[index + 5] += 1e6*Sxz;
				}
			}
		}

		if(processFaultSerially) {
			// Concatenate the partition array into the full patch
			// linearized tensor array, since the fault has been
			// processed serially.
			for(size_t k = 0; k < partitionedTensorSize; ++k) {
				coeffs_st[k]  = coeffs_st_partitioned[k];
				coeffs_dip[k] = coeffs_dip_partitioned[k];
				coeffs_open[k] = coeffs_open_partitioned[k];
			}
		}
		else {
			MPI_Allgather(coeffs_st_partitioned, partitionedTensorSize,
						  MPI_FLOAT, coeffs_st, partitionedTensorSize,
						  MPI_FLOAT, MPI_COMM_WORLD);

			MPI_Allgather(coeffs_dip_partitioned, partitionedTensorSize,
						  MPI_FLOAT, coeffs_dip, partitionedTensorSize,
						  MPI_FLOAT, MPI_COMM_WORLD);

			MPI_Allgather(coeffs_open_partitioned, partitionedTensorSize,
						  MPI_FLOAT, coeffs_open, partitionedTensorSize,
						  MPI_FLOAT, MPI_COMM_WORLD);
		}

		free(coeffs_st_partitioned);
		free(coeffs_dip_partitioned);
		free(coeffs_open_partitioned);

		// Copy data from the linearized tensors to the f3arrays.

		int linearIndex = 0, tensorIndex = 0;

		// Calculate tensorIndex
		for(size_t fault = 0; fault < j; ++fault) {
			// Index should start just after all the patches that
			// have already been processed for previous faults
			tensorIndex += eqkfm1[fault].np_di*eqkfm1[fault].np_st;
		}

		for(int i = 0; i < numPatches; ++i) {
			for(int j = 0; j < Nsel; ++j) {
				for(int k = 0; k < 6; ++k) {
					linearIndex = (i * Nsel * 6) + (j * 6) + k;

					if (flag_sslip) (*Coeffs_st) [tensorIndex + i + 1][j+1][k+1] = coeffs_st [linearIndex];
					if (flag_dslip) (*Coeffs_dip)[tensorIndex + i + 1][j+1][k+1] = coeffs_dip[linearIndex];
					if (flag_open) (*Coeffs_open)[tensorIndex + i + 1][j+1][k+1] = coeffs_open[linearIndex];
				}
			}
		}

		free(coeffs_st);
		free(coeffs_dip);
		free(coeffs_open);
	}

	return(0);
}
#endif

int okadaCoeff(float ****Coeffs_st, float ****Coeffs_dip, float ****Coeffs_open,
		struct eqkfm *eqkfm1, int NF, struct crust crst) {

	/* Allocates and populates arrays containing the okada coefficients between each fault and grid point.
	 * This function is where most of the memory used by the program in most cases is allocated, and it may still be optimized (see comments).
	 *
	 * Input:
	 *  eqkfm1: structure containing slip models. Range [0...NF-1].
	 *  crst: crust structure, containing position of grid points and elastic parameters.
	 */

	double north, east, eqnorth, eqeast;
	double len, width, depth; //for individual patches.
	double depth0; //to differentiate between blind fault (depth0=0) or fault cutting through surface.
	double strike, dip, rake;
	double Sxx, Syy, Szz, Sxy, Syz, Sxz;
	int Nsel=eqkfm1[0].nsel;
	int NP_tot=0, p1, i, noslip_str, noslip_dip, noopen;
	int err=0;
	int flag_open, flag_sslip, flag_dslip;
	struct eqkfm *afslip;	//dummy variable equal to pointers in eqkfm1 which point to afterslip elements.

	for (int j=0; j<NF; j++) NP_tot+=eqkfm1[j].np_di*eqkfm1[j].np_st;

	depth0=eqkfm1[0].cuts_surf ? eqkfm1[0].top : 0.0;

	print_logfile("Depth of surface: %.3lf km.\n", depth0);
	print_screen("Depth of surface: %.3lf km.\n", depth0);
	

	//---------initialize DCFS----------//

	//check if Coeffs tensors should be allocated.
	flag_open=flag_sslip=flag_dslip=0;

	for (int j=0; j<NF; j++){

		//if the element is associated with afterslip, should check whether this afterslip has strike/slip/opening component.
		afslip=	eqkfm1[j].co_aft_pointer;

		//check if there is no slip in this fault:
		noslip_str= (eqkfm1[j].slip_str==NULL && (afslip==NULL || (*afslip).allslip_str==NULL));
		noslip_dip= (eqkfm1[j].slip_dip==NULL && (afslip==NULL || (*afslip).allslip_dip==NULL));
		noopen= (eqkfm1[j].open==NULL && (afslip==NULL || (*afslip).allslip_open==NULL));

		flag_open=max(flag_open, !noopen);
		flag_sslip=max(flag_sslip, !noslip_str);
		flag_dslip=max(flag_dslip, !noslip_dip);

		if (flag_open && flag_sslip && flag_dslip) break;
	}

	*Coeffs_st= (flag_sslip) ? f3array(1,NP_tot,1,Nsel,1,6) : NULL;
	*Coeffs_dip= (flag_dslip) ? f3array(1,NP_tot,1,Nsel,1,6) : NULL;
	*Coeffs_open= (flag_open) ? f3array(1,NP_tot,1,Nsel,1,6) : NULL;

	for (int i=1; i<=NP_tot; i++){
		for (int j=1; j<=Nsel; j++){
			for (int k=1; k<=6; k++){
				if (flag_sslip) (*Coeffs_st)[i][j][k]=0;
				if (flag_dslip) (*Coeffs_dip)[i][j][k]=0;
				if (flag_open) (*Coeffs_open)[i][j][k]=0;
			}
		}
	}

	//-----------------------------------------------------------------------------------------//
	//-----------Calculate Coulomb stress vector for each patch assuming slip=1.---------------//
	//-----------------------------------------------------------------------------------------//

	print_screen("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);
	print_logfile("Calculating Okada solutions (%d patches, %d grid points)...\n", NP_tot, Nsel);

	p1=0;	//count total number of patches (in all faults);
	for (int j=0; j<NF; j++){

		if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
			print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
			return(1);
		}

		len=eqkfm1[j].L*(1.0/eqkfm1[j].np_st);
		width=eqkfm1[j].W*(1.0/eqkfm1[j].np_di);

		for (int p=1; p<=eqkfm1[j].np_di*eqkfm1[j].np_st; p++){
			p1+=1;
			patch_pos(eqkfm1[j], p, &eqeast, &eqnorth, &depth);

			#pragma omp parallel for private(Sxx, Syy, Szz, Sxy, Syz, Sxz, north, east, i, afslip, noslip_str, noslip_dip, noopen)
			for (int i0=1; i0<=Nsel; i0++){
				i=eqkfm1[0].selpoints[i0];
				north=crst.y[i];
				east=crst.x[i];

				//if the element is associated with afterslip, should check whether this afterslip has strike/slip/opening component.
				afslip=	eqkfm1[j].co_aft_pointer;

				//check if both coseismic and afterslip have empty arrays for each component of deformation:
				noslip_str= (eqkfm1[j].slip_str==NULL && (afslip==NULL || (*afslip).allslip_str==NULL));
				noslip_dip= (eqkfm1[j].slip_dip==NULL && (afslip==NULL || (*afslip).allslip_dip==NULL));
				noopen= (eqkfm1[j].open==NULL && (afslip==NULL || (*afslip).allslip_open==NULL));

				if (!noslip_str) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike,  dip, len, width, 1.0, 0.0, 0.0, north, east, crst.depth[i]-depth0,
							&Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz,  crst.lambda, crst.mu);

					(*Coeffs_st)[p1][i0][1]+=1e6*Sxx;
					(*Coeffs_st)[p1][i0][2]+=1e6*Syy;
					(*Coeffs_st)[p1][i0][3]+=1e6*Szz;
					(*Coeffs_st)[p1][i0][4]+=1e6*Sxy;
					(*Coeffs_st)[p1][i0][5]+=1e6*Syz;
					(*Coeffs_st)[p1][i0][6]+=1e6*Sxz;
				}

				if (!noslip_dip) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike, dip, len, width, 0.0, -1.0, 0.0, north, east, crst.depth[i]-depth0,
							&Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz, crst.lambda, crst.mu);

					(*Coeffs_dip)[p1][i0][1]+=1e6*Sxx;
					(*Coeffs_dip)[p1][i0][2]+=1e6*Syy;
					(*Coeffs_dip)[p1][i0][3]+=1e6*Szz;
					(*Coeffs_dip)[p1][i0][4]+=1e6*Sxy;
					(*Coeffs_dip)[p1][i0][5]+=1e6*Syz;
					(*Coeffs_dip)[p1][i0][6]+=1e6*Sxz;
				}

				if (!noopen) {
					pscokada(eqnorth, eqeast, depth-depth0,  strike, dip, len, width, 0.0, 0.0, 1.0, north, east, crst.depth[i]-depth0,
							&Sxx, &Syy, &Szz, &Sxy, &Syz, &Sxz, crst.lambda, crst.mu);

					(*Coeffs_open)[p1][i0][1]+=1e6*Sxx;
					(*Coeffs_open)[p1][i0][2]+=1e6*Syy;
					(*Coeffs_open)[p1][i0][3]+=1e6*Szz;
					(*Coeffs_open)[p1][i0][4]+=1e6*Sxy;
					(*Coeffs_open)[p1][i0][5]+=1e6*Syz;
					(*Coeffs_open)[p1][i0][6]+=1e6*Sxz;
				}
			}
		}
	}

	return(0);
}

int okadaCoeff2DCFS(float ***Coeffs_st, float ***Coeffs_d, float ***Coeffs_open, struct pscmp DCFS, struct eqkfm *eqkfm1){
	/* Calculates stress tensors (DCFS) given a set of okada coefficients (Coeffs_XX) and a slip/opening model (eqkfm1).
	 *
	 * Input:
	 *  Coeffs_XX: okada coefficients, relating unit displacement on a fault patch to stress tensor components at a grid point.
	 *  eqkfm1: structure containing slip models.
	 *
	 * Output:
	 *  DCFS: will be populated with stress tensors at each grid point.
	 */

	double strike, dip, rake;
	int p1, p, j;
	int Nsel;
	int NF=DCFS.NF;
	int err=0, errp=0;

	if ((DCFS.nsel!=(*eqkfm1).nsel)){
		print_screen("**Warning: DCFS.nsel!=eqkfm.nsel in okadaCoeff2DCFS. Will use those from eqkfm1. **\n");
		print_logfile("**Error: DCFS.nsel!=eqkfm1.nsel (%d, %d) in okadaCoeff2DCFS.**\n", DCFS.nsel, (*eqkfm1).nsel);
		DCFS.nsel=(*eqkfm1).nsel;
		DCFS.which_pts=(*eqkfm1).selpoints;
	}

	Nsel=DCFS.nsel;

	for (int i0=1; i0<=DCFS.nsel; i0++) {
		DCFS.cmb[i0]=0.0;
		for (int a=1; a<=3; a++){
				for (int b=1; b<=3; b++) DCFS.S[i0][a][b]=0;
		}
	}

	#pragma omp parallel for private(p1, p, j, strike, dip, rake) reduction(+:errp)
	for (int i=1; i<=Nsel; i++){

		//-----------------------------------------------------------------------------------------//
		//-----------Calculate Coulomb stress vector for each patch and add them up.---------------//
		//-----------------------------------------------------------------------------------------//

		p1=0;
		for (j=0; j<NF; j++){

			if ((err=choose_focmec(eqkfm1[j], &strike, &dip, &rake))!=0){
				print_screen("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
				print_logfile("*** Illegal value for eqkfm[%d].whichfm (okadaDCFS) ***\n",j);
				errp+=1;
			}

			for (p=1; p<=eqkfm1[j].np_di*eqkfm1[j].np_st; p++){
				p1+=1;
				if (eqkfm1[j].slip_str){
					DCFS.S[i][1][1]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][1];
					DCFS.S[i][1][2]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][4];
					DCFS.S[i][1][3]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][6];
					DCFS.S[i][2][2]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][2];
					DCFS.S[i][2][3]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][5];
					DCFS.S[i][3][3]+=eqkfm1[j].slip_str[p]*Coeffs_st[p1][i][3];
				}
				if (eqkfm1[j].slip_dip){
					DCFS.S[i][1][1]+=eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][1];
					DCFS.S[i][1][2]+=eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][4];
					DCFS.S[i][1][3]+=eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][6];
					DCFS.S[i][2][2]+=eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][2];
					DCFS.S[i][2][3]+=eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][5];
					DCFS.S[i][3][3]+=eqkfm1[j].slip_dip[p]*Coeffs_d[p1][i][3];
				}
				if (eqkfm1[j].open){
					DCFS.S[i][1][1]+=eqkfm1[j].open[p]*Coeffs_open[p1][i][1];
					DCFS.S[i][1][2]+=eqkfm1[j].open[p]*Coeffs_open[p1][i][4];
					DCFS.S[i][1][3]+=eqkfm1[j].open[p]*Coeffs_open[p1][i][6];
					DCFS.S[i][2][2]+=eqkfm1[j].open[p]*Coeffs_open[p1][i][2];
					DCFS.S[i][2][3]+=eqkfm1[j].open[p]*Coeffs_open[p1][i][5];
					DCFS.S[i][3][3]+=eqkfm1[j].open[p]*Coeffs_open[p1][i][3];
				}
			}
		}
		DCFS.S[i][2][1]=DCFS.S[i][1][2];
		DCFS.S[i][3][2]=DCFS.S[i][2][3];
		DCFS.S[i][3][1]=DCFS.S[i][1][3];

	}

	return(errp!=0);
}

int isoDCFS(struct pscmp DCFS, struct eqkfm eqkfm1){

	/* Calculates radially decaying, isotropic stress changes given a seismic source (eqkfm1)
	 *
	 * Input:
	 *  eqkfm1: contains information about the earthquake source.
	 *
	 * Output:
	 *  DCFS.cmb contains Coulomb stress changes.
	 */

	double M0, r;
	int Nsel;
	double DCFSmax=DCFS_cap;

	if (DCFS.nsel!=eqkfm1.nsel){
		print_logfile("**Error: DCFS.nsel!=eqkfm.nsel in isoDCFS.**\n");
		print_screen("**Warning: DCFS.nsel!=eqkfm.nsel in isoDCFS. Will choose the one from DCFS.**\n");
		eqkfm1.nsel=DCFS.nsel;
		eqkfm1.selpoints=DCFS.which_pts;
	}

	Nsel=DCFS.nsel;
	#pragma omp parallel for private(r, M0)
	for (int i=1; i<=Nsel; i++){
		r=DCFS.fdist[i];	//r is in km.
		M0=pow(10,1.5*(eqkfm1.mag+6.0));
		DCFS.cmb[i]=(M0/(6.0*PI))*pow(1000*r,-3.0);
		if (fabs(DCFS.cmb[i])>DCFSmax) DCFS.cmb[i]=(DCFS.cmb[i]>0)? DCFSmax : -DCFSmax;
	}
	return(0);
}


//---------------------------------------------------------------------//
//-----					Auxiliary functions						  -----//
//---------------------------------------------------------------------//


int choose_focmec(struct eqkfm eqkfm1, double *strike, double *dip, double *rake){
//assign strike, dip, rake correct value from eqkfm, and convert to radians.

	switch (eqkfm1.whichfm){
		case 1:
			if (strike) *strike=DEG2RAD*eqkfm1.str1;
			if (dip) *dip=DEG2RAD*eqkfm1.dip1;
			if (rake) *rake=DEG2RAD*eqkfm1.rake1;
			break;
		case 2:
			if (strike) *strike=DEG2RAD*eqkfm1.str2;
			if (dip) *dip=DEG2RAD*eqkfm1.dip2;
			if (rake) *rake=DEG2RAD*eqkfm1.rake2;
			break;
		case 0:
			if (strike) *strike=DEG2RAD*eqkfm1.str1;
			if (dip) *dip=DEG2RAD*eqkfm1.dip1;
			if (rake) *rake=DEG2RAD*eqkfm1.rake1;
			break;
		default:
			return(1);
	}
	return (0);
}

void patch_pos(struct eqkfm eqfm, int p, double *east, double *north, double *depth){
//find position of patch in local cartesians.

	double pos_top, dz, dip, strike;

	choose_focmec(eqfm, &strike, &dip, NULL);

	pos_top=eqfm.pos_d[p]-0.5*eqfm.W*(1.0/eqfm.np_di);	//position of top of the patch (along dip).
	*east=eqfm.x+eqfm.pos_s[p]*sin(strike)+pos_top*cos(dip)*cos(strike);
	*north=eqfm.y+eqfm.pos_s[p]*cos(strike)-pos_top*cos(dip)*sin(strike);
	dz=sin(dip)*pos_top;	//depth of top of the patch.
	*depth=eqfm.depth+dz;

	return;
}

double *normal_vector(double strikeR, double dipR){
//calculates normal vector to a plane from its strike, dip.
	double *n;

	n=darray(1,3);
	n[1]=-sin(strikeR*DEG2RAD)*sin(dipR*DEG2RAD);
	n[2]=cos(strikeR*DEG2RAD)*sin(dipR*DEG2RAD);
	n[3]=-cos(dipR*DEG2RAD);

	return n;
}

double *slip_vector(double strikeR, double dipR, double rakeR){
//calculates slip vector from its strike, dip, rake.

	double *s;
	s=darray(1,3);
	s[1]=cos(rakeR*DEG2RAD)*cos(strikeR*DEG2RAD)+sin(rakeR*DEG2RAD)*sin(strikeR*DEG2RAD)*cos(dipR*DEG2RAD);
	s[2]=cos(rakeR*DEG2RAD)*sin(strikeR*DEG2RAD)-sin(rakeR*DEG2RAD)*cos(strikeR*DEG2RAD)*cos(dipR*DEG2RAD);
	s[3]=-sin(rakeR*DEG2RAD)*sin(dipR*DEG2RAD);

	return s;
}

double *opt_s(double *stress, double sigma, double *n, double *result){
	/* returns direction of maximum shear stress on a plane normal to vector n, subject to stress "stress" and normal stress sigma.
	 * if result ==NULL, ignore it and allocate new memory. Otherwise, memory should have been allocated previously.
	 *
	 * stress[1,2,3]= stress vector;
	 * sigma=normal stress (compression is positive);
	 * n=normal vector to plane;
	 * */

	double s_temp[4], *s;
	double s_tempM;

	s= (result) ? result : darray(1,3);

	for (int g=1; g<4; g++) s_temp[g]=stress[g]+sigma*n[g];
	s_tempM=norm(s_temp,3);
	for (int g=1; g<4; g++) s[g]= (s_tempM==0)? s_temp[g] : s_temp[g]*(1.0/s_tempM);

	return s;
}

double *sum_v(double *v1, double *v2, double *sum, int N){
// sums two vectors [1...N].
//if sum==NULL, it is ignored and new memory is allocated.

	double *v3;

	if (!v2) return v1;
	if (!v1) return v2;

	v3= (sum) ? sum : darray(1,N);
	if (v1 && v2) {
		for (int n=1; n<=N; n++) v3[n]=v1[n]+v2[n];
	}
	else {
		if (!v1) v3=v2;
		if (!v2) v3=v1;
	}
	return v3;
}

double resolve_S(double **S, double strikeR, double dipR, double rakeR, double f, double *stress0, double sigma0, int opt_rake){
//Resolves stress tensor S on a mechanism given by strikeR, dipR, rakeR, and calculates the Coulomb stress.
//total stress should be passed if optimal rake is to be used (opt_rake==1).
//resolves stress tensor S on mechanisms with dipR, strikeR and rake.
//sigma0, stress0 are not calculated here for efficiency reasons (they don't change as often as S).

	double *n, *s;
	double cmb;

	n=normal_vector(strikeR, dipR);
	if (opt_rake) s=NULL;
	else s=slip_vector(strikeR, dipR, rakeR);

	cmb=resolve_n(S, n, f, stress0, sigma0, s);

	free_darray(n,1,3);
	if (s) free_darray(s,1,3);

	return cmb;
}

double resolve_n(double **S, double *n, double fric, double *stress0, double sigma0, double *slip_v){
//resolves stress S on focal mechanism with plane normal to vector n[1..3], and slip direction vector s[1..3].
//if slip_v==NULL, will use optimal rake. This is calculated from stress tensor S as well as background stresses stress0, and background normal pressure sigma0.
//is sigma0==0, stress0==NULL, no background stress.
//sigma0, stress0 are not calculated here for efficiency reasons (they don't change as often as S).

	double *s, *stress, *stress_tot=NULL;
	double sigma, tau;
	int optrake= !slip_v;

	s=darray(1,3);
	stress=darray(1,3);

	mtimesv(S,n,stress,3,3);
	sigma=-1.0*vdotv(stress,n,3);	// compression is positive (rock mechanics convention).

	if (optrake) {
		opt_s(stress_tot=sum_v(stress, stress0, NULL, 3), sigma+sigma0, n, s);
		tau=vdotv(stress,s,3);
	}
	else tau=vdotv(stress,slip_v,3);

	free_darray(s,1,3);
	free_darray(stress,1,3);
	if (optrake) free_darray(stress_tot,1,3);

	return tau-fric*sigma;

}


