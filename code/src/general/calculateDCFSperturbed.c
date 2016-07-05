
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



#include "calculateDCFSperturbed.h"

#include <math.h>
//#include <stddef.h>
#include <stdlib.h>

#include "../defines.h"
//#include "../inp_out/write_csep_forecast.h"
#include "../okada/okadaDCFS.h"
#include "../seis/cmbopt.h"
#include "../util/error.h"
#include "../util/moreutil.h"
#include "../util/util1.h"
#include "mem_mgmt.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

//float ran1(long *);

void calculateDCFSperturbed(double **DCFSrand, struct pscmp *DCFS, struct eqkfm *eqkfmAf,
							struct eqkfm *eqkfm0, struct flags flag,
							double *times, int Nmain, int NA, struct crust crst,
							struct Coeff_LinkList *AllCoeff,
							struct Coeff_LinkList *AllCoeffaft, int NTScont,
							double **focmec, int *fmzoneslim, int NFM,
							double tdata0, double tdata1, int refresh, int which_recfault) {

/* Input:
 * 
 * eqkfmAf contains aseismic slip snapshots;
 * eqkfm0 contains seismic sources;
 * times: times to which elements of tevol correspond: S(t=times[j])=S0*tevol[j], where S=slip, S0 is the slip contained in eqkfmAf.
 * NTScont is the total number of time steps for continuous process;
 * Nmain: length of eqkfm0;
 * NA: length of eqkfmAf;
 * AllCoeff: okada coefficients for mainshocks (events in eqkfm0);
 * focmec contains sample of focal mechanisms, NFM is its length (1->NFM)
 * fmzoneslin: indices of focmec corresponding to limits of distinct foc. mec. areas;
 * tdata0, tdata1: start and end time for which data should be used;
 * refresh: flag to be set to 1 if the slip models have changed from previous function call;
 * flag contains various flags, including:
 * 		aseismic and aftershocks are flags indicating if these processes should be included.
 * 		vary_recfault is a flag indicating which received faults to use. 0: use fix planes; 1: vary foc. mec. (use which_fm, or if which_fm==0, choose random one); 2: use OOPs.
 * which_recfault gives index of foc.mec to select; if set to 0, choose random one (useful if NFM>>nsur). if vary_recfault=0, this is meaningless.
 *
 *
 * Output:
 *
 * DCFSrand[i][j] contains the ith stress change at gridpoint j due to continuous processes (modeled linearly between time steps).
 * DCFS[k].cmp[j] contains the stress change due to kth event at gridpoint j (modeled as a step).
 *
 * NB: if fixed receiver faults are used, it will still choose the optimal rake. This can be changed with a flag below.
 *
 */

	// Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	//flags:
	int	aseismic=flag.aseismic, \
		vary_recfault=flag.err_recfault, \
		gridpoints_err=flag.err_gridpoints, \
		multisnap=flag.aseismic_multisnap, \
		full_field=(flag.sources_without_focmec==2), \
		optimal_rake=flag.optrake;

	static double *strike0, *dip0, *rake0;	//strike0, dip0, rake0 are focal mechanism that will change at each iteration.
	double slip;
	double ***Stemp;
	static double lat0, lon0;
	int rand;
	int NgridT=crst.N_allP;
	int NFsofar=0;
	int last, first, i;
	int a_ev=0;	//counter for DCFS_Af
	int nfaults=0;	//counter for eqkfmAf
	int NTSeff;
	static int fm_offset=0;	//offset to be added to crst.str0 (dip0) if a fixed receiver fault is used (more details below).
	static struct eqkfm *eqkfm2;
	static struct eqkfm *eqkfm2A;
	struct Coeff_LinkList *temp, *temp2;
	float ***Coeffs_st, ***Coeffs_dip, ***Coeffs_open;	//coefficients for tensor;
	static float **Coeff_ResS, **Coeff_ResD;	//coefficients for cmb (resolved).
	static struct pscmp *DCFS_Af;
	int DCFS_Af_size;
	int NF_max=0, NP_max=0;
	static int time_in=0;
	static int **nn;	//nearest neighbours points.
	static double **interp_DCFS;
	static double *mycmb=NULL;
	static double **cmb_cumu;
	int n_withslimodel, n_withoutslimodel;

	int cap_aseismic=0;	// This flag determines whether aseismic slip should be capped if its absolute value exceeds DCFS_cap (global variable).

//	FILE *fout;
	time_in+=1;

	//----------------------------------------------------------------------//
	// 		 		This initial part will be executed only once			//
	//----------------------------------------------------------------------//


	if (time_in==1) {
		print_logfile("\nSetting up variables for calculating Coulomb stress fields...\n");
		print_screen("\nSetting up variables for calculating Coulomb stress fields...\n");

		//-----calculate neighbouring points--------//

		if (gridpoints_err==1) {
			nn=i2array(1,NgridT,1,6);
			nearest_neighbours(NgridT, crst.nLat,crst.nLon,crst.nD, nn);
			interp_DCFS=d2array(1,NgridT,1,2);
			mycmb=darray(1,NgridT);
		}


		if (vary_recfault==1){
			//allocate memory to strike0, dip0, rake0.
			//if vary_recfault==1, receiver fault changes at each MC iteration.
			if (which_recfault==0) {
			//different focal mechanisms selected if MC iterations of focal mechanism should be used. However, if which_recfault!=0, a single mechanism should be used (focmec[X][which_recfault]).
				strike0=darray(0,crst.nofmzones-1);
				dip0=darray(0,crst.nofmzones-1);
				rake0=darray(0,crst.nofmzones-1);
			}

			else {
				//in this case a single foc. mec. is needed.
				strike0=malloc(sizeof(double));
				dip0=malloc(sizeof(double));
				rake0=malloc(sizeof(double));
			}
		}
		else{
			// if a fixed mechanism is used, the receiver fault from crst.str0, crst.dip0 should be used.
			// crst.str0[0] contains the regional mechanism, that should be used if no spatially variable mech. is given ((*crst).variable_fixmec=0).
			// crst.str0[1...NP] contain the foc. mech. for individual grid points, which should be used if (*crst).variable_fixmec=0.
			// This is achieved by passing crst.str0+fm_offset to the function "resolve_DCFS" later on.
			fm_offset= crst.variable_fixmec ? 1 : 0;
		}


		//-----------------------------------------------------------------//
		//					Afterslip (initial setup)					   //
		//-----------------------------------------------------------------//

		if (aseismic!=0){

			NTSeff=(multisnap)? NTScont : 1;
			DCFS_Af_size= NTSeff*NA;
			DCFS_Af= pscmp_array(0,DCFS_Af_size);

			a_ev=0;	//counter for DCFS_Af
			nfaults=0;	//counter for eqkfmAf

			if (cap_aseismic) {
				cmb_cumu=d2array(0,NA-1,1,NgridT);
				if (!cmb_cumu) memory_error_quit;
				for (int a=0; a<NA; a++){
					for (int n=1; n<=NgridT; n++) cmb_cumu[a][n]=0.0;
				}
			}

			for (int a=0; a<NA; a++){

				for (int i=0; i<NTSeff; i++){
					DCFS_Af[a_ev+i].NF=AllCoeffaft->NF;
					DCFS_Af[a_ev+i].cmb= darray(1,NgridT);	//only allocated stuff needed by OkadaCoeff2....
					DCFS_Af[a_ev+i].S=d3array(1,NgridT,1,3,1,3);
					DCFS_Af[a_ev+i].nsel=eqkfmAf[nfaults].nsel;
					DCFS_Af[a_ev+i].which_pts=eqkfmAf[nfaults].selpoints;
				}

				Coeffs_st=AllCoeffaft->Coeffs_st;
				Coeffs_dip=AllCoeffaft->Coeffs_dip;
				Coeffs_open=AllCoeffaft->Coeffs_open;

				for (int i=0; i<NTSeff; i++)	{
					for (int nf=0; nf<DCFS_Af[a_ev].NF; nf++){
						//assign slip values for each snapshot and fault into slip_X (this is needed because okadaCoeff2DCFS uses them):
						eqkfmAf[nfaults+nf].slip_str= (eqkfmAf[nfaults+nf].allslip_str)? eqkfmAf[nfaults+nf].allslip_str[i] : NULL;
						eqkfmAf[nfaults+nf].slip_dip= (eqkfmAf[nfaults+nf].allslip_dip)? eqkfmAf[nfaults+nf].allslip_dip[i] : NULL;
						eqkfmAf[nfaults+nf].open= (eqkfmAf[nfaults+nf].allslip_open)? eqkfmAf[nfaults+nf].allslip_open[i] : NULL;
					}
					okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, Coeffs_open, DCFS_Af[a_ev+i], eqkfmAf+nfaults);
					if (vary_recfault==0) {
							if (optimal_rake) resolve_DCFS_nocap(DCFS_Af[a_ev+i], crst, crst.str0+fm_offset, crst.dip0+fm_offset, NULL, 1);
							else resolve_DCFS_nocap(DCFS_Af[a_ev+i], crst, crst.str0+fm_offset, crst.dip0+fm_offset, crst.rake0+fm_offset, 0);

						free_d3array(DCFS_Af[a_ev+i].S, 1,NgridT,1,3,1,3);

							if (cap_aseismic){
								for (int n=1; n<=NgridT; n++) {
									if (multisnap) cmb_cumu[a_ev/NTSeff][n]+=DCFS_Af[a_ev+i].cmb[n];
									else {
									    for (int l=0; l<NTScont; l++) {
										cmb_cumu[a_ev][n]+=DCFS_Af[a_ev].cmb[n]*eqkfmAf[i].tevol[l];
									    }
									}
								}
						}
					}
				}
				nfaults+=DCFS_Af[a_ev].NF;	//counter for eqkfmAf
				a_ev+=NTSeff;	//counter for DCFS_Af
				AllCoeffaft=AllCoeffaft->next;
			}
		}
		print_screen("done.\n");
	}

	//-----------------------------------------------------------------//
	//					Mainshock(initial setup)					   //
	//(this part executed every time a new set of slip models is used) //
	//-----------------------------------------------------------------//

	if (time_in==1 || refresh){
		if (aseismic!=2){
			NFsofar=0;
			temp=AllCoeff;
			for (int i=0; i<Nmain; i++){

				//Don't do anything if event is outside data time period:
				if (DCFS[temp->which_main].t <tdata0 || DCFS[temp->which_main].t>tdata1){
					NFsofar+=temp->NF;
					temp=temp->next;
					continue;
				}

				if (eqkfm0[NFsofar].is_slipmodel){
					//calculate stress tensor at each grid point:
					Coeffs_st=temp->Coeffs_st;
					Coeffs_dip=temp->Coeffs_dip;
					Coeffs_open=temp->Coeffs_open;
					okadaCoeff2DCFS(Coeffs_st, Coeffs_dip, Coeffs_open, DCFS[temp->which_main], eqkfm0+NFsofar);

					//resolve coefficients if receiver faults don't change between iterations:
					switch (vary_recfault){
						case 0:
							if (optimal_rake) resolve_DCFS(DCFS[temp->which_main], crst, crst.str0+fm_offset, crst.dip0+fm_offset, NULL, 1);
							else resolve_DCFS(DCFS[temp->which_main], crst, crst.str0+fm_offset, crst.dip0+fm_offset, crst.rake0+fm_offset, 0);
							if (gridpoints_err){
								int eq1=temp->which_main;

								//need to copy field into cmb0, so that each iteration will smooth the original field (and not the one smoothed in prev. iteration).
								for (int ii=1; ii<=NgridT; ii++) mycmb[ii]=0.0;
								for (int ii=1; ii<=DCFS[eq1].nsel; ii++) mycmb[DCFS[eq1].which_pts[ii]]=DCFS[eq1].cmb[ii];
								interp_nn(NgridT,crst.nLat, crst.nLon, crst.nD, mycmb,interp_DCFS,0,nn);
								DCFS[eq1].cmb0=darray(1,DCFS[eq1].nsel);
								DCFS[eq1].Dcmb=darray(1,DCFS[eq1].nsel);
								for (int ii=1; ii<=DCFS[eq1].nsel; ii++){
									DCFS[eq1].cmb0[ii]=0.5*(interp_DCFS[DCFS[eq1].which_pts[ii]][1]+interp_DCFS[DCFS[eq1].which_pts[ii]][2]);
									DCFS[eq1].Dcmb[ii]=fabs(interp_DCFS[DCFS[eq1].which_pts[ii]][1]-interp_DCFS[DCFS[eq1].which_pts[ii]][2]);
								}
							}

							break;
						case 2:
							DCFScmbopt(DCFS, temp->which_main, crst);	//NB this does not take into account stress from aseismic slip, assuming that stresses from mainshocks are much larger.
							break;
						default:
							break;
					}

				}

				else{
					//prepare isotropic stress fields:
					if (!eqkfm0[NFsofar].is_slipmodel){
						int eq1=temp->which_main;
						isoDCFS(DCFS[eq1], eqkfm0[NFsofar]);
						if (gridpoints_err){
							for (int ii=1; ii<=NgridT; ii++) mycmb[ii]=0.0;
							for (int ii=1; ii<=DCFS[eq1].nsel; ii++) mycmb[DCFS[eq1].which_pts[ii]]=DCFS[eq1].cmb[ii];
							interp_nn(NgridT,crst.nLat, crst.nLon, crst.nD, mycmb,interp_DCFS,0,nn);
							DCFS[eq1].cmb0=darray(1,DCFS[eq1].nsel);
							DCFS[eq1].Dcmb=darray(1,DCFS[eq1].nsel);
							for (int ii=1; ii<=DCFS[eq1].nsel; ii++){
								DCFS[eq1].cmb0[ii]=0.5*(interp_DCFS[DCFS[eq1].which_pts[ii]][1]+interp_DCFS[DCFS[eq1].which_pts[ii]][2]);
								DCFS[eq1].Dcmb[ii]=fabs(interp_DCFS[DCFS[eq1].which_pts[ii]][1]-interp_DCFS[DCFS[eq1].which_pts[ii]][2]);
							}
						}
					}
				}
				NFsofar+=temp->NF;
				temp=temp->next;
			}
		}
	}


	//----------------------------------------------------------------------//
	//	 		From here on will be executed at every iteration			//
	//----------------------------------------------------------------------//


	//At the moment aseismic slip and OOPs at the same time are not implemented.
	if (aseismic==1 && vary_recfault==2) {
		print_screen("*Error: function calculateDCFSperturbed doesn't know how to calculate OOPS when aseismic slip is included!!*\n");
		print_logfile("*Error: function calculateDCFSperturbed doesn't know how to calculate OOPS when aseismic slip is included!!*\n");
		return;
	}

	//pick receiver faults from catalog of focal mechanisms:
	if (vary_recfault==1){
		//if vary_recfault==1, receiver fault changes at each MC iteration.
		if (which_recfault==0) {
			//randomly pick a mechanism for each zone:
			for (int fmzone=0; fmzone<crst.nofmzones; fmzone++){
				first=fmzoneslim[fmzone];
				last=fmzoneslim[fmzone+1]-1;
				rand= (int) ((last-first)*ran1()+first);
				strike0[fmzone]=focmec[1][rand];
				dip0[fmzone]=focmec[2][rand];
				rake0[fmzone]=focmec[3][rand];	//only used for splines==1 (see below).
			}
		}
		else {
			// use the foc. mec. for this iteration (this is done when all focal mechanisms should be sampled; only activated in main.c if nofmzones=1).
			*strike0=focmec[1][which_recfault];
			*dip0=focmec[2][which_recfault];
			*rake0=focmec[3][which_recfault];	//only used for splines==1 (see below).
		}
	}

	//----------------------------------------------//
	//	calculated stress field from aseismic slip: //
	//----------------------------------------------//

	if (aseismic==0){
		for (int l=0; l<NTScont; l++){
			if (times[l] <tdata0 || times[l]>tdata1) continue;
			for (int n=1; n<=NgridT; n++) if (DCFSrand) DCFSrand[l][n]=0.0;
		}
	}
	else {
		i=0;

		for (int l=0; l<NTScont; l++) {
			if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
			for (int n=1; n<=NgridT; n++) if (DCFSrand) DCFSrand[l][n]=0.0;
		}

		for (int a=0; a<NA; a++){

			if (multisnap==0){

				if (vary_recfault==1) {
					if (optimal_rake) resolve_DCFS_nocap(DCFS_Af[a], crst, strike0, dip0, NULL, 1);
					else resolve_DCFS_nocap(DCFS_Af[a], crst, strike0, dip0, rake0, 0);
				}

				for (int l=0; l<NTScont; l++) {
					if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;

					//loop over events with aseismic slip:
					for (int n=1; n<=NgridT; n++) {
						if (cap_aseismic) DCFSrand[l][n]+= (fabs(cmb_cumu[a][n])>DCFS_cap) ? 
										(DCFS_cap/fabs(cmb_cumu[a][n]))*eqkfmAf[i].tevol[l]*DCFS_Af[a].cmb[n] : eqkfmAf[i].tevol[l]*DCFS_Af[a].cmb[n];
						else DCFSrand[l][n]+=eqkfmAf[i].tevol[l]*DCFS_Af[a].cmb[n];
					}
				}
				i+=DCFS_Af[a].NF;
			}

			else{
				if (vary_recfault==1){
					if (cap_aseismic) for (int n=1; n<=NgridT; n++) cmb_cumu[a][n]=0.0;
					for (int l=0; l<NTScont; l++) {
						if (optimal_rake) resolve_DCFS_nocap(DCFS_Af[a*NTScont+l], crst, strike0, dip0, NULL, 1);
						else resolve_DCFS_nocap(DCFS_Af[a*NTScont+l], crst, strike0, dip0, rake0, 0);
						if (cap_aseismic && l<NTScont-1) for (int n=1; n<=NgridT; n++) cmb_cumu[a][n]+=DCFS_Af[a*NTScont+l].cmb[n];
					}
				}

				for (int l=0; l<NTScont; l++) {
					if ((l>0 && times[l-1]) <tdata0 || (l<NTScont-1 && times[l+1]>tdata1)) continue;
					for (int n=1; n<=NgridT; n++) {
						if (cap_aseismic) DCFSrand[l][n]+= (fabs(cmb_cumu[a][n])>DCFS_cap) ? (DCFS_cap/fabs(cmb_cumu[a][n]))*DCFS_Af[a*NTScont+l].cmb[n] : DCFS_Af[a*NTScont+l].cmb[n];
						else DCFSrand[l][n]+= DCFS_Af[a*NTScont+l].cmb[n];
					}
				}
			}
		}
	}


	//------------------------------------------------------------------------------------------------------//
	//	calculated stress field from mainshocks (i.e. events for which a non trivial slip model is used):   //
	//------------------------------------------------------------------------------------------------------//

	temp=AllCoeff;
	NFsofar=0;
	for (int i=0; i<Nmain; i++){

		if (DCFS[temp->which_main].t <tdata0 || DCFS[temp->which_main].t>tdata1){
			NFsofar+=temp->NF;
			temp=temp->next;
			continue;
		}

		else{
			if (aseismic==2){
				for (int n=1; n<=NgridT; n++) DCFS[temp->which_main].cmb[n]=0.0;	//don't need to fill in tensor.
			temp=temp->next;
			}
			else{
				if (eqkfm0[NFsofar].is_slipmodel){
					if (vary_recfault==1){	// In this case stress tensor has to be resolved on a different plane at each iteration:
						if (optimal_rake) resolve_DCFS(DCFS[temp->which_main], crst, strike0, dip0, NULL, 1);
						else resolve_DCFS(DCFS[temp->which_main], crst, strike0, dip0, rake0, 0);
					}
					if (gridpoints_err==1) {
						//if vary_recfault==0, the field hasn't changed and cmb0 should be used.
						smoothen_DCFS(DCFS[temp->which_main], crst.nLat, crst.nLon, crst.nD, vary_recfault==0, nn);
					}
				}
				else {
					//if the earthquake does not have a slip model, the isotropic field does not change between iterations and has already been calculated.
					//only need to calculated error from gridpoints.
					if (gridpoints_err==1) smoothen_DCFS(DCFS[temp->which_main], crst.nLat, crst.nLon, crst.nD, 1, nn);
				}

				NFsofar+=temp->NF;
				temp=temp->next;
			}
		}
	}
}

void smoothen_DCFS(struct pscmp DCFS, int nlat, int nlon, int nd, int use_cmb0, int **nn){
	/* Calculates perturbed stress field by estimating the gradients within each cell from the difference with neighbor ones.
	 * Only works with regular grids.
	 *
	 * Input
	 *  DCFS: already contains stress field (cmb)
	 *  nlat, nlon, nd: number of elements along each dimension.
	 *  use_cmb0: flag. if set to 0, the values in DCFS.cmb are perturbed and overwritten.
	 *  				if set to 1, it will take the values from DCFS.cmb0 and perturb them within a precalculated range given by DCFS.Dcmb.
	 *  						With this option, DCFS.cmb0= mean value of field; DCFS.cmb=range of values. This method is faster, and
	 *  						DCFS.cmb0 is not overwritten: useful if grid point smoothing is only source of uncertainty.
	 *
	 * Output:
	 *  DCFS.cmb is populated.
	 */


	double randcmb;
	double *mycmb, **interp_DCFS;
	int NgridT=nlat*nlon*nd;

	if (!use_cmb0){
		mycmb=darray(1,NgridT);
		interp_DCFS=d2array(1,NgridT,1,2);
		for (int i=1; i<=NgridT; i++) mycmb[i]=0.0;
		for (int i=1; i<=DCFS.nsel; i++) mycmb[DCFS.which_pts[i]]=DCFS.cmb[i];

		interp_nn(NgridT, nlat, nlon, nd, mycmb, interp_DCFS, 0, nn);

		for (int i=1; i<=DCFS.nsel; i++) {
			int p=DCFS.which_pts[i];
			randcmb=interp_DCFS[p][1]+ran1()*(interp_DCFS[p][2]-interp_DCFS[p][1]);
			DCFS.cmb[i]=randcmb;
		}

		free_darray(mycmb,1,NgridT);
		free_d2array(interp_DCFS,1,NgridT,1,2);
	}
	else{
		for (int i=1; i<=DCFS.nsel; i++) {
			randcmb=-0.5*DCFS.Dcmb[i]+ran1()*DCFS.Dcmb[i];
			DCFS.cmb[i]=DCFS.cmb0[i]+randcmb;
		}
	}
}
