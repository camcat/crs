
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

#include "setup.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../defines.h"
#include "../geom/check_same_geometry.h"
#include "../geom/coord_trafos.h"
#include "../geom/top_of_slipmodel.h"
#include "../inp_out/read_eqkfm.h"
#include "../inp_out/read_focmec.h"
#include "../inp_out/read_zmap.h"
#include "../okada/okadaDCFS.h"
#include "../util/error.h"
#include "../util/merge.h"
#include "../util/moreutil.h"

#include "../util/util1.h"
#include "../util/splines_eqkfm.h"
#include "eqkfm_copy.h"
#include "find_timesteps.h"
#include "mem_mgmt.h"
#include "struct_conversions.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int setup_catalogetc(char *catname, char **focmeccat, int nofmcat,
					 struct tm reftime, double dDCFS, double Mag_source, double Mag_main, struct crust crst,
					 struct catalog *cat, struct eqkfm **eqkfm1, double ***focmec,
					 int **firstelements, struct flags flag, int *NFM, int *Ntot,
					 double dt, double dM, double xytoll, double ztoll, double dR, double tw,
					 double tstart, double tend) {

/*	Reads catalog and focal mechanisms catalog, and fills cat, eqkfm structures accordingly.
 *
 * Input:
 *
 * catname:	ZMAP catalog file
 * focmeccat:	list of focal mechanisms catalog files [0...nofmcat-1]
 * nofmcat:	number of focal mechanisms catalogs
 * reftime: IssueTime (times will be calculated with reference to this time)
 * dDCFS:	min. value for which grid points should be selected for calculating stress changes from source events
 * Mag_source: min magnitude of earthquakes to be used as sources of stress.
 * Mag_main:	magnitude of mainshocks (i.e. events which are included as sources also if flags.aftershocks==0)
 * crst:	structure containing domain information
 * flag:	flags structure
 * dt, dM, xytoll, ztoll: max. expected difference (tolerance) between events from catalog and focmec catalog;
 * dR: extra distance to be considered for sources.
 * tw:	time window to be ignored for event selection after each mainshock (NB only ignored in cat, still included as sources in eqkfm).
 * tstart: start time for including sources and catalog events
 * tend: end time for forecast (sourced only included up to IssueTime, i.e. t=0)
 *
 * Output:
 *
 * cat:	catalog structure used for LL calculations
 * eqkfm1:	eqkfm structure containing stress sources. Can be NULL, and will be ignored.	[0...Ntot-1]
 * focmec:	array containing focal mechanisms [1...NFM]
 * first_elements:	indices of focmec elements which correspond to the first element of a new focal mechanism area (i.e. a new focal mechanisms catalog)
 * NFM:	length of focmec
 * Ntot: length of eqkfm1
 * Nmain:	number of mainshocks in eqkfm1
 *
 */

	struct eqkfm *eqkfm1fm;
	double tendS, tendCat;
	double minmag;
	int err=0, NgridT=crst.N_allP;
	int Nfm;

	print_screen("Setting up catalog...\n");

	tendS=0;		//this is the "IssueDate", up to which data can be used for sources.
	tendCat=tend;	//this is the largest between "ForecastEndDate" and "tendLL", up to which data is available. Events after t=0 used to calculate LL of forecast.

	//select events within some tolerance level, since they will have to be matched with focal mechanisms.
	err += readZMAP(cat, eqkfm1, Ntot, catname, crst, reftime, tstart, tendS, tstart, tendCat,
			Mag_main, tw, fmax(xytoll, dR), fmax(ztoll, dR), dDCFS, 1);

	if (err) return (err);

	// read catalog of focal mechanisms and merge it with eqkfm structure from catalog:
	if (focmeccat && (!flag.sources_all_iso || focmec)){
		err+=readmultiplefocmec(focmeccat, nofmcat, crst,fmax(xytoll, dR), fmax(ztoll, dR), dDCFS,
				reftime, tstart, tendS, tendS, (*cat).Mc, focmec, firstelements, NFM, &Nfm,  &eqkfm1fm, 1, 1);
		combine_eqkfm(*eqkfm1, eqkfm1fm, *Ntot, Nfm, dt, dM, xytoll, 1);
	}

	// filter eqkfm according to magnitude, depth.
	minmag= Mag_source;
	eqk_filter(eqkfm1, Ntot, minmag , crst.depmax+fmax(dR,ztoll));

	// calculate distances between source events and grid points.
	eqkfm2dist((*eqkfm1), crst.lat, crst.lon, crst.depth, NgridT, *Ntot, 1);


	print_screen("%d events used for LL calculation, %d events used as sources.\n", (int) (*cat).Z, *Ntot);
	print_logfile("%d events used for LL calculation, %d events used as sources.\n", (int) (*cat).Z, *Ntot);

	return (err!=0);
}

int setup_aseismic_eqkfm(struct slipmodels_list list_slipmodels, struct crust crst, struct eqkfm **eqkfm0res){
/*
 * Reads aseismic files into eqkfm structure.
 *
 * Input:
 * 	list_slipmodel: list of slip model files
 * 	crst: structure containing crust setup information
 *
 * Output:
 *  eqkfm0res: contains all slip models for aseismic processes.
 */

	int Nm;
	double *tmain=list_slipmodels.tmain;	//contains NSM elements (starting time of each aseismic event). NSM=no. of aseismic events.
	double *tsnap=list_slipmodels.tsnap;	//contains nmtot elements, where nmtot=sum(no_slipmodels): times of all slip model files
	char **slipmodels=list_slipmodels.slipmodels;	//contains nmtot elements (all slip model files)
	double *disc=list_slipmodels.disc;		//contains NSM elements (discretization for each aseismic event)
	int *Nfaults=list_slipmodels.Nfaults;	//contains NSM elements (no. of subfaults for each aseismic event)
    int err=0;
    int counter=0;
    int nsm=0, totfaults=0;
    char *cmb_format=list_slipmodels.cmb_format;

    //Find tot. no. of faults by summing over aseismic events:
    for (int N=0; N<list_slipmodels.NSM; N++){

    	//number of snapshots for current aseismic event.
		Nm=list_slipmodels.no_slipmodels[N];

		//read the no. of faults from the first snapshot in each aseismic event ('counter'):
		if (!(strcmp(cmb_format,"farfalle"))) err+=read_farfalle_eqkfm(slipmodels[counter], NULL, Nfaults+N);
		else {
			if (!(strcmp(cmb_format,"pscmp"))) err+=read_pscmp_eqkfm(slipmodels[counter], NULL, Nfaults+N);
			else {
				if (!(strcmp(cmb_format,"fsp"))) err+=read_fsp_eqkfm(slipmodels[counter], NULL, Nfaults+N);
				else {
					print_logfile("Unknown slip model format %s (setup_aseismic_eqkfm).\n", cmb_format);
					return 1;
				}
			}
		}
		if (err) return(err);
		counter+=Nm;
		totfaults+=Nfaults[N];
    }

    //eqkfm0res contains Nfaults[N] elements for each aseismic event N (since the geometry is fixed for all snapshots in the same event). Each element contains the snapshots inside.
    *eqkfm0res=eqkfm_array(0,totfaults-1);
    counter=0;
    totfaults=0;
    for (int N=0; N<list_slipmodels.NSM; N++){
		Nm=list_slipmodels.no_slipmodels[N];	//number of snapshots for current aseismic event.
    	err+=setup_aseismic_element(*eqkfm0res+totfaults, slipmodels+counter, cmb_format, Nm, crst.mu, disc[N], tmain[N], tsnap+counter,
    			crst.N_allP, crst.list_allP, list_slipmodels.cut_surf[counter], crst.lat0, crst.lon0);

		counter+=Nm;
		totfaults+=Nfaults[N];
    }

    return(err);
}


int setup_aseismic_element(struct eqkfm *eqkfm0res, char **slipmodels, char *cmb_format, int no_snap,
						double mu, double disc, double tmain, double *tsnap, int nsel,
						int *sel_pts, int cuts_surf,
						double lat0, double lon0) {

	/* Set up the structure describing a single aseismic event. An event can have multiple subfaults and multiple snapshots.
	 *
	 * Input:
	 *  slipmodels: list of file names corresponding to each snapshot. size [0...no_snap-1].
	 *  cmbformat:	format of input slip model.
	 *  no_snap: no. of snapshots in time.
	 *  mu: shear modulus (used to calculate seismic moments and magnitude)
	 *  disc: final slip model discretization
	 *  tmain: starting time
	 *  tsnap: time of snapshots
	 *  nsel: no. of grid points affected by this event
	 *  sel_points: indices of grid points affected by this event (relative to point lists in crust structure)
	 *  cuts_suft: flag indicating whether the slip model should be assumed to cut the surface
	 *  lat0, lon0: reference lat/lon (used to calculate x,y of slip model in the local coordinate system).
	 *
	 * Output:
	 *  eqkfm0res is populated. size [0...Nfaults] (Nfaults is read from the slip model below).
	 */

	//todo to save memory, may also want to allow individual snapshots of allslip_xxx to be NULL (e.g. eqkfm0res[f].allslip_xxx[t]=NULL if there is no slip at time t).
	int err=0, NF;
	double 	toll=1e-10;
	struct eqkfm *eqkfm0;	//used to read individual slip models (one per snapshot), later copied into eqkfm0res structure.
	double ***allslip_str_temp,***allslip_dip_temp, ***allslip_open_temp;	//store slip values from individual snapshots.
	int is_str, is_dip, is_open;	//flags used to determine which components of deformations are needed (to save memory).

	//read the last snapshot to find NF and initialize temporary variable eqkfm0.
	err=read_eqkfm(slipmodels[no_snap-1], cmb_format, &eqkfm0, &NF, NULL, mu);	//find NF and eqkfm0[x].np_x
	if (err) return (err);

	// allocate temporary storage (since eqkfm0 gets overwritten)
	allslip_str_temp=malloc(NF*sizeof(double **));
	allslip_dip_temp=malloc(NF*sizeof(double **));
	allslip_open_temp=malloc(NF*sizeof(double **));
	for (int nf=0; nf<NF; nf++) {
		allslip_str_temp[nf]=d2array(0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
		allslip_dip_temp[nf]=d2array(0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
		allslip_open_temp[nf]=d2array(0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
		if (!allslip_str_temp | !allslip_dip_temp |! allslip_open_temp) memory_error_quit;

	}

	//read in values for slip:
	for (int m=0; m<no_snap; m++){
		err=read_eqkfm(slipmodels[m], cmb_format, &eqkfm0, &NF, NULL, mu);
		for (int nf=0; nf<NF; nf++){
			copy_vector(eqkfm0[nf].slip_str, &(allslip_str_temp[nf][m]), eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			copy_vector(eqkfm0[nf].slip_dip, &(allslip_dip_temp[nf][m]), eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			copy_vector(eqkfm0[nf].open, &(allslip_open_temp[nf][m]), eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			free_darray(eqkfm0[nf].slip_str,1,0);
			free_darray(eqkfm0[nf].slip_dip,1,0);
			free_darray(eqkfm0[nf].open,1,0);
		}

		if (err) return (err);
	}

	//Needs to do it this way since eqkfm0[nf].tot_slip is overwritten at each read_eqkfm call.
	for (int nf=0; nf<NF; nf++){
		free_darray(eqkfm0[nf].tot_slip,0,0);
		eqkfm0[nf].tot_slip=darray(0,no_snap-1);
		for (int m=0; m<no_snap; m++){
			eqkfm0[nf].tot_slip[m]=eqkfm0[nf].tot_slip[0];
		}
	}

	if (cuts_surf) {
			top_of_slipmodel(eqkfm0, NF);
			for (int i=0; i<NF; i++) eqkfm0[i].cuts_surf=1;
	}

	for (int nf=0; nf<NF; nf++) {
		eqkfm0[nf].tot_slip=darray(0,no_snap-1);
		eqkfm0[nf].ts=darray(1,no_snap);	//shifted by one element because of indexing in copy_vector function.
		copy_vector(tsnap-1, &(eqkfm0[nf].ts), no_snap);	//copy aseismic time steps into eqkfm0 structure.
		eqkfm0[nf].ts+=1;	//since should start from 0th element (not 1st).
		eqkfm0[nf].nosnap=no_snap;
		eqkfm0[nf].t=tmain;
		eqkfm0[nf].nsel=nsel;
		eqkfm0[nf].selpoints=sel_pts;
		eqkfm0[nf].is_slipmodel=1;
		latlon2localcartesian(eqkfm0[nf].lat, eqkfm0[nf].lon, lat0, lon0, &(eqkfm0[nf].north), &(eqkfm0[nf].east));

		//check if all elements are 0, and is so set flag.
		is_str=0; is_dip=0; is_open=0;
		for (int t=0; t<=no_snap-1; t++){
			for (int p=1; p<=eqkfm0[nf].np_st*eqkfm0[nf].np_di; p++){
				if (fabs(allslip_str_temp[nf][t][p])>toll) is_str=1;
				if (fabs(allslip_dip_temp[nf][t][p])>toll) is_dip=1;
				if (fabs(allslip_open_temp[nf][t][p])>toll) is_open=1;

				if (is_str && is_dip && is_open) break;
			}
			if (is_str && is_dip && is_open) break;
		}

		//free memory if elements are all 0.
		if (is_str==0){
			free_d2array(allslip_str_temp[nf],0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			allslip_str_temp[nf]=NULL;
		}
		if (is_dip==0){
			free_d2array(allslip_dip_temp[nf],0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			allslip_dip_temp[nf]=NULL;
		}
		if (is_open==0){
			free_d2array(allslip_open_temp[nf],0,no_snap-1,1, eqkfm0[nf].np_st*eqkfm0[nf].np_di);
			allslip_open_temp[nf]=NULL;
		}

		eqkfm0[nf].allslip_str=allslip_str_temp[nf];
		eqkfm0[nf].allslip_dip=allslip_dip_temp[nf];
		eqkfm0[nf].allslip_open=allslip_open_temp[nf];

		copy_eqkfm_all(eqkfm0[nf], eqkfm0res+nf);
		eqkfm0res[nf].parent_set_of_models=NULL;
	}
	
	return err;
}

int setup_eqkfm_element(struct eqkfm *eqkfm0res, char **slipmodels, char *cmb_format, int no_slipmodels,
						double mu, double tmain, int nsel,
						int *sel_pts, double *mmain, int cuts_surf,
						int *NF0, double lat0, double lon0, int same_geometry) {

	/* Sets up elements of eqkfm0res corresponding to one earthquake.
	 * It may contain several slip models, and several subfaults (with each subfault equal to one struct eqkfm element).
	 *
	 * Input:
	 *  slipmodels: list of files, each containing one slip model. Range: [0...no_slipmodels]
	 *  cmbformat: format of slip model files (pscmp/fsp/farfalle)
	 *  mu: shear modulus (used to calculate magnitude from the slip)
	 *  tmain: time of the eartquake
	 *  nsel: no. of grid points associated with this earthquake.
	 *  sel_pts: list of grid points associated with this earthquake.
	 *  cuts_surf: flag indicating if the model should be forced to cut through the surface.
	 *  lat0, lon0: earthquake coordinates.
	 *  same_geometry: if multiple slip models are given, flag indicating whether they have the same geometry.
	 *
	 * Output:
	 *  eqkfm0res: slip models structure.
	 *  mmain: earthquake magnitude.
	 *  NF0: no. of subfaults (if there are multiple slip models, the largest between them).
	 *
	 */

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct set_of_models setmodels;
	struct eqkfm *eqkfm0;
	int err=0, NF, nftot=0, nfmax=0;
	double 	toll=1e-10;

	setmodels.NF_models=iarray(1,no_slipmodels);
	setmodels.Nmod=no_slipmodels;
	setmodels.current_model=1;
	setmodels.same_geometry=same_geometry;

	for (int m=0; m<no_slipmodels; m++){
		err=read_eqkfm(slipmodels[m], cmb_format, NULL, &NF, NULL, mu);
		if (err){
			print_screen(" ** Error: Input slip model %s could not be read (setup_eqkfm_element). **\n", slipmodels[m]);
			print_logfile(" ** Error: Input slip model %s could not be read (setup_eqkfm_element). **\n", slipmodels[m]);
			return (1);
		}
		setmodels.NF_models[m+1]=NF;
		nfmax=fmax(nfmax,NF);
		nftot+=NF;
	}
	setmodels.NFmax=nfmax;
	setmodels.set_of_eqkfm=eqkfm_array(0,nftot-1);

	nftot=0;
	for (int m=0; m<no_slipmodels; m++){
		err=read_eqkfm(slipmodels[m], cmb_format, &eqkfm0, &NF, mmain, mu);
		if (err) return (err);

		if (cuts_surf) {
			top_of_slipmodel(eqkfm0, NF);
			for (int i=0; i<NF; i++) eqkfm0[i].cuts_surf=1;
		}

		reduce_eqkfm_memory(eqkfm0, NF);	//delete empty slip_str, slip_dip, open arrays.
		for (int nf=0; nf<NF; nf++) {
			eqkfm0[nf].t=tmain;
			eqkfm0[nf].nsel=nsel;
			eqkfm0[nf].selpoints=sel_pts;
			eqkfm0[nf].is_slipmodel=1;
			latlon2localcartesian(eqkfm0[nf].lat, eqkfm0[nf].lon, lat0, lon0, &(eqkfm0[nf].north), &(eqkfm0[nf].east));
			setmodels.set_of_eqkfm[nftot+nf]=eqkfm0[nf];
		}
		nftot+=NF;
	}

	//allocate memory and copy values from setmodels;
	(*eqkfm0res).parent_set_of_models=(struct set_of_models *) malloc((size_t) (sizeof(struct set_of_models)));
	memcpy((*eqkfm0res).parent_set_of_models, &setmodels, (size_t) sizeof(struct set_of_models));

	set_current_slip_model(eqkfm0res,1);
	if (NF0) *NF0=nfmax;
	return err;
}

void set_current_slip_model(struct eqkfm *eqkfm0, int slipmodel_index){
/* Sets variables in eqkfm0 to required slip model.
 *
 * Input:
 *  eqkfm0: contains the slip models. In particular, eqkfm0[0].parent_set_of_models contains pointers to all the alternative slip models (eqkfm0[0].parent_set_of_models.set_of_eqkfm).
 *  slipmodel_index: the index of the slip model that should be set.
 *
 * Output:
 *  The required slip model is activated by setting structures inside eqkfm0 to those stored in "parent_set_of_models".
 */

	// [Fahad] Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int nftot=0;
	struct set_of_models allmod= *(eqkfm0[0].parent_set_of_models);
	struct eqkfm* eqkfmall=allmod.set_of_eqkfm;

	// allmod.Nmod==0 means that no slip model is used for this event: do nothing.
	if (!allmod.Nmod) return;

	// find where slip model no.slipmodel_index starts;
	// copy slip models from eqkfmall to eqkfm0;
	// delete all info from remaining elements of eqkfm0;
	for (int i=1; i<slipmodel_index; i++) nftot+=allmod.NF_models[i];
	for (int nf=0; nf<allmod.NF_models[slipmodel_index]; nf++) copy_eqkfm_all(eqkfmall[nf+nftot], eqkfm0+nf);
	for (int nf=allmod.NF_models[slipmodel_index]; nf<allmod.NFmax; nf++) empty_eqkfm(eqkfm0+nf);

	print_logfile("Slip model set to no. %d.\n", slipmodel_index);

	return;
}

int setup_CoeffsDCFS(struct Coeff_LinkList **Coefficients, struct Coeff_LinkList **Coefficients_aseismic, struct pscmp **DCFS_out,
		struct crust crst, struct eqkfm *eqkfm0, int Nm, int *Nfaults, struct eqkfm *eqkfm_aft, int no_aseismic, int *Nfaults_aft) {

/* Setup structures which will contain okada coefficients.
 * Also matches aseismic slip with seismic slip, and uses pointers to avoid repeating okada calculations if the models have the same geometry (e.g. for afterslip).
 *
 *  crst: crust structure needed to initialize DCFS.
 *  eqkfm0: seismic slip models. Range [0...NFtot-1], where NFtot=sum(Nfaults).
 *  Nfaults: no. of faults for each of the seismic sources. Range [0...Nm-1].
 *
 *  eqkfm_aft: aseismic slip models. Range [0...NFtot-1], where NFtot=sum(Nfaults_aft).
 *  Nfaults_aft: no. of faults for each of the seismic sources. Range [0...No_aseismic-1].
 *
 *
 * Output:
 *  Coefficients, Coefficients_aseismic structures are initialized. The actual coefficients will be calculated (and memory allocated) in okadaDCFS.c.
 *  DCFS (which will contain stresses from seismic sources) is initialized.
 */

	// [Fahad] Variables used for MPI
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	struct pscmp *DCFS;
    struct Coeff_LinkList *AllCoeff, *temp, *temp2, *AllCoeff_aseismic;
    struct set_of_models som_tmp;
    int Nsel, NFtot, NFtotaft;
    int mainshock_withafterslip, same_geometry;
    int aseismic= (eqkfm_aft==NULL) ? 0 : 1;
    double M0;

    //----------set up Coefficients----------------//

    if (Coefficients!=NULL) {
    	//Create elements of structure (allocate memory):
    	AllCoeff= malloc( sizeof(struct Coeff_LinkList));
		temp= AllCoeff;
		for(int i=0; i<Nm; i++) {
			temp->borrow_coeff_from=NULL;	// coefficients must be be calculated for all seismic sources.
			temp->which_main=i;
			if (i<Nm-1) {
				temp->next= malloc(sizeof(struct Coeff_LinkList));
				temp= temp->next;
			}
			else {
				temp->next=(struct Coeff_LinkList *) 0;
			}
		}

		*Coefficients=AllCoeff;
		print_logfile("Okada Coefficients structure set up.\n");
    }

    //--------------set up DCFS-------------------//

    if (DCFS_out!=NULL){

		DCFS=pscmp_arrayinit(crst,0,Nm-1);

		NFtot=0;
		for (int eq=0; eq<Nm; eq++){
			Nsel = eqkfm0[NFtot].nsel;
			DCFS[eq].fdist=eqkfm0[NFtot].distance;
			DCFS[eq].index_cat=eqkfm0[NFtot].index_cat;
			DCFS[eq].which_pts=eqkfm0[NFtot].selpoints;
			DCFS[eq].t=eqkfm0[NFtot].t;
			M0=0.0;
			for (int f=0; f<Nfaults[eq]; f++) M0+=pow(10,1.5*(eqkfm0[NFtot+f].mag+6));
			DCFS[eq].m=(2.0/3.0)*log10(M0)-6;
			DCFS[eq].NF=Nfaults[eq];
			if (DCFS[eq].nsel!=Nsel){
				if (DCFS[eq].nsel>0){
					free_d3array(DCFS[eq].S,1,DCFS[eq].nsel,1,3,1,3);
					free_darray(DCFS[eq].cmb,1,DCFS[eq].nsel);
				}
				DCFS[eq].nsel=Nsel;
				DCFS[eq].S = d3array(1,Nsel,1,3,1,3);
				DCFS[eq].cmb=darray(1,Nsel);
				for (int i=1; i<=Nsel; i++) DCFS[eq].cmb[i]=0.0;
			}
			NFtot+=Nfaults[eq];
		}

		*DCFS_out=DCFS;
		print_logfile("DCFS structure set up.\n");
    }

    //--------------associates afterslip with mainshocks-------------------//
	// uses a ~1 sec tolerance

    if (aseismic){

    	AllCoeff_aseismic= malloc( sizeof(struct Coeff_LinkList));
    	temp=AllCoeff_aseismic;

		NFtot=NFtotaft=0;
    	for (int a=0; a<no_aseismic; a++){
    		//check whether the aseismic slip is associated with coseismic slip.
    		//this is done so that the same okada coefficients are used: less memory/computations are needed.
    		mainshock_withafterslip=closest_element(timesfrompscmp(DCFS, Nm), Nm, eqkfm_aft[NFtotaft].t, 0.000011575);

    		//if a mainshock is found, calculate respective position in eqkfm0 structure (NFtot) and in AllCoeff.
    		if (mainshock_withafterslip!=-1){
    			NFtot=0;
				for (int i=0; i<mainshock_withafterslip; i++) NFtot+=Nfaults[i];
				//if multiple slip models have different geometries, should use its own set of coefficients. Otherwise, check if the geometry is the same as afterslip.
				som_tmp= *(eqkfm0[NFtot].parent_set_of_models);

				if (som_tmp.same_geometry==0 && som_tmp.Nmod>1) {
					same_geometry= 0;
				}
				else{
					same_geometry= check_same_geometry(eqkfm0+NFtot, Nfaults[mainshock_withafterslip], eqkfm_aft+NFtotaft, Nfaults_aft[a]);
				}
    		}

			if (mainshock_withafterslip==-1 || !same_geometry){

				for (int nf=0; nf<Nfaults_aft[a]; nf++){
					eqkfm_aft[NFtotaft+nf].co_aft_pointer=NULL;
				}

				temp->borrow_coeff_from=NULL;
				temp->Coeffs_dip=NULL;
				temp->Coeffs_st=NULL;
				temp->Coeffs_open=NULL;

			}
			else{
				//shift to correct AllCoeff element:
    			temp2=AllCoeff;
    			while (temp2 && temp2->which_main!=mainshock_withafterslip) temp2=temp2->next;

				//set pointers to corresponding co/postseismic elements to each others:
				for (int nf=0; nf<Nfaults[mainshock_withafterslip]; nf++){
					eqkfm0[NFtot+nf].co_aft_pointer=eqkfm_aft+NFtotaft+nf;
					eqkfm_aft[NFtotaft+nf].co_aft_pointer=eqkfm0+NFtot+nf;
				}

				//cuts_surf must be the same for afterslip and coseismic model:
				for (int f=0; f<Nfaults_aft[a]; f++) eqkfm_aft[NFtotaft+f].cuts_surf=eqkfm0[NFtot+f].cuts_surf;

				//set pointers of corresponding elements to each others. Later this will be used to avoid calculating Okada coeff. twice.
				temp->borrow_coeff_from=temp2;
				temp2=temp2->next;
			}

			temp->next= malloc(sizeof(struct Coeff_LinkList));
			temp=temp->next;
			NFtotaft+=Nfaults_aft[a];
		}

		*Coefficients_aseismic=AllCoeff_aseismic;
		print_logfile("Okada Coefficients structure set up for %d aseismic sources.\n", no_aseismic);

    }

    return(0);
}

int update_CoeffsDCFS(struct Coeff_LinkList **Coefficients,
		struct crust crst, struct eqkfm *eqkfm0, int Nm, int *Nfaults) {

	/* Calculates Okada coefficients.
	 * It should be called once at the beginning for each Coeff_LinkList structure, and also after switching slip models in eqkfm0.
	 *
	 * Input:
	 *  eqkfm0: uses slip model geometry to calculate the Okada coefficients. Range [0...NFtot-1], where NFtot=sum(Nfaults).
	 *  		Nfaults: no. of faults for each of the seismic sources. Range [0...Nm-1].
	 */


	// [Fahad] Variables used for MPI
	int procId = 0, numProcs = 1;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	#endif

    struct Coeff_LinkList *temp, *temp2;
    struct set_of_models tmp_setofmodels;
    int NFsofar=0;
    static int first_timein=1, switch_slipmodel;

	//Fill in elements of structure:
	temp= *Coefficients;

	for(int i=0; i<Nm; i++) {
		if (temp->borrow_coeff_from==NULL){	//coefficients should be calculated.
			if (eqkfm0[NFsofar].is_slipmodel) {
				// coefficients should only calculated if more than one slip model is provided (switch_slipmodel):
				// for aseismic slip, multiple slip models are not allowed: parent_set_of_models=NULL
				if (eqkfm0[NFsofar].parent_set_of_models){
					tmp_setofmodels=*(eqkfm0[NFsofar].parent_set_of_models);
					switch_slipmodel= (tmp_setofmodels.Nmod > 1);
				}
				else{
					//with aseismic slip, this function will only be entered once (but first_timein=0 since it was entered with seismic slip), and coefficients should be populated.
					switch_slipmodel=1;
				}

				if (first_timein || switch_slipmodel){
					if (!first_timein){
						if (temp->Coeffs_st) free_f3array(temp->Coeffs_st, 1,0,1,0,1,0);
						if (temp->Coeffs_dip) free_f3array(temp->Coeffs_dip, 1,0,1,0,1,0);
						if (temp->Coeffs_open) free_f3array(temp->Coeffs_open, 1,0,1,0,1,0);
					}

					#ifdef _CRS_MPI
						okadaCoeff_mpi(&(temp->Coeffs_st), &(temp->Coeffs_dip), &(temp->Coeffs_open), eqkfm0+NFsofar, Nfaults[i], crst);
					#else
						okadaCoeff(&(temp->Coeffs_st), &(temp->Coeffs_dip), &(temp->Coeffs_open), eqkfm0+NFsofar, Nfaults[i], crst);
					#endif

					temp->NgridT=eqkfm0[0].nsel;
					temp->NF=Nfaults[i];
					temp->NP=0;
					for(int f=0; f<Nfaults[i]; f++) {
						temp->NP+= eqkfm0[NFsofar+f].np_di*eqkfm0[NFsofar+f].np_st;
					}
					temp->which_main=i;
				}
			}
			else {
				temp->Coeffs_st=temp->Coeffs_dip=temp->Coeffs_open=NULL;
				if (first_timein){
					temp->NgridT=eqkfm0[0].nsel;
					temp->NF=Nfaults[i];
					temp->NP=0;
					for(int f=0; f<Nfaults[i]; f++) {
						temp->NP+= eqkfm0[NFsofar+f].np_di*eqkfm0[NFsofar+f].np_st;
					}
					temp->which_main=i;
				}
			}
		}

		// This is used to avoid calculating the Okada coefficients multiple times if two model (e.g. seismic/aseismic) have the same geometry.
		else{
			temp2=temp->borrow_coeff_from;
			temp->NF=temp2->NF;				// tot. no of faults;
			temp->which_main=temp2->which_main;		// index of pscmp DCFS to which earthquake refer;
			temp->NP=temp2->NP;				// tot. no. of patches (sum of no. of patches of individual faults);
			temp->NgridT=temp2->NgridT;			// no. of grid cells.
			temp->Coeffs_st=temp2->Coeffs_st;
			temp->Coeffs_dip=temp2->Coeffs_dip;
			temp->Coeffs_open=temp2->Coeffs_open;	// Coefficient for strike slip, dip slip displacements.
		}

		NFsofar+=Nfaults[i];


		if (i<Nm-1) {
			temp= temp->next;
		}
	}
	print_logfile("Okada Coefficients structure updated.\n");

	first_timein=0;

    return(0);
}
