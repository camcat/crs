
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


#include "read_eqkfm.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

#include "../defines.h"
#include "../general/eqkfm_copy.h"
#include "../general/mem_mgmt.h"
#include "../general/setup.h"
#include "../general/struct_conversions.h"
#include "../seis/soumod1.h"
#include "../seis/WellsCoppersmith.h"

#include "../util/util1.h"
#include "read_eqkfm_fsp.h"

int eqkfm_addslipmodels(struct eqkfm *eqfm1, struct slipmodels_list all_slipmodels, struct eqkfm **eqfm_comb,
						int N1, int *Ncomb, int **nfout,
						double dt, double dmag, double res, struct crust crst, struct flags flags) {

/* Matches the earthquakes in eqfm1 (presumably from an earthquake catalog) with available slip models. Also takes care of events which do not have a slip model,
 * by creating a synthetic slip from the focal mechanism if possible.
 *
 * Input:
 *  eqfm1: structure containing all sources. Each source is one element of the array. Range [0...N1-1]
 *  all_slipmodels: structure containing the slip models.
 *  dt, dmag: tolerances for matching times and magnitude between catalog and slip model list.
 *  res: final slip model resolution of synthetic slip models derived from a focal mechanism.
 *  crst: crust structure containing grid point information.
 *  flags: contains flags about what to do for events without a slip model and without a focal mechanism available.
 *
 *  eqfm_comb:	structure containing all sources, with final slip models.
 *  			If some of the slip model has multiple subfaults, this array is longer than the total no. of sources. Range [0...NFtot-1], where NFtot=sum(nfout).
 *  nfout: contains the number of subfaults of each source (i.e. no. of eqfm_comb elements associated with it). Range [0...Ncomb-1]
 *  Ncomb: no. of earthquake sources. May be smaller than N1 if some events are excluded because they don't have slip model or foc mec (depends on flags).
 */

	int *which_slipmod;
	int N2, N3=0, c_evfound=0;
	int *nf2, nsm, nfaults;
	int c2=0, c3=0;	//counters.
	int err=0, err0, j;
	int *all_pts, p;
	int no_synthetic_slipmodels=0, no_slipmodels=0;
	static struct set_of_models dummy_parentsetofmodels;
	char *cmb_format=all_slipmodels.cmb_format;
	dummy_parentsetofmodels.Nmod=0;	//this value indicates than no slip model is available (will use synthetic slip model from foc. mec. or isotropic field).

	print_logfile("\nBuilding mainshock slip models (eqkfm_addslipmodels).\n");

	N2=all_slipmodels.NSM;
	nf2=all_slipmodels.Nfaults;
	all_pts=iarray(1,crst.N_allP);
	for (int i=1; i<=crst.N_allP; i++) all_pts[i]=i;

	*Ncomb=0;

	//find no. of faults in each slip model, and save the largest value for each event.
	nsm=0;
	for (int i=0; i<all_slipmodels.NSM; i++){
		nfaults=0;
		all_slipmodels.is_used[i]=0;

		for (int n=1; n<=all_slipmodels.no_slipmodels[i]; n++){
			if (!(strcmp(cmb_format,"farfalle"))) err+=read_farfalle_eqkfm(all_slipmodels.slipmodels[nsm], NULL, all_slipmodels.Nfaults+i);
			else {
				if (!(strcmp(cmb_format,"pscmp"))) err+=read_pscmp_eqkfm(all_slipmodels.slipmodels[nsm], NULL, all_slipmodels.Nfaults+i);
				else {
					if (!(strcmp(cmb_format,"fsp"))) err+=read_fsp_eqkfm(all_slipmodels.slipmodels[nsm], NULL, all_slipmodels.Nfaults+i);
					else {
						print_logfile("Unknown slip model format %s (eqkfm_addslipmodels).\n", cmb_format);
						err=1;
					}
				}
			}
			nsm+=1;
			nfaults=fmax(nfaults,all_slipmodels.Nfaults[i]);
		}
		all_slipmodels.Nfaults[i]=nfaults;
	}

	if (err){
		print_screen("Error in reading input files. Exiting.\n");
		print_logfile("Error in reading input files. Exiting.\n");
		return 1;
	}

	//[check] make sure mmain is already assigned here!! (NB: all this may not be needed if clear association b/w cat and foc mec is given).
	which_slipmod = combine_cats(all_slipmodels.tmain, timesfromeqkfm(eqfm1, N1, (int *) 0),
								 all_slipmodels.mmain, magssfromeqkfm(eqfm1, N1, (int *) 0 ),
								 all_slipmodels.NSM, N1, dt, dmag);

	//calculate tot. no. of faults needed:
	for (int i=0; i<N1; i++){
		N3= (which_slipmod[i]==-1)? N3+1 : N3+ nf2[which_slipmod[i]];
	}

	*eqfm_comb=eqkfm_array(0,N3-1);
	if (nfout) *nfout=iarray(0,N1-1);

	//todo speed up this loop (try to parallelize?)
	for(int i=0; i<N1; i++) {
		(*eqfm_comb)[c3].nsel=crst.N_allP;

		if(which_slipmod[i]==-1) {	//not associated with a complex slip model
			if (flags.sources_without_focmec==0 && !eqfm1[i].is_slipmodel){
				// flags.sources_without_focmec==0 indicates that only events with focal mechanisms should be used.
				continue;
			}
			else {
				copy_eqkfm_all(eqfm1[i], (*eqfm_comb)+c3);
				eqkfm2dist((*eqfm_comb)+c3, crst.lat, crst.lon, crst.depth, crst.N_allP, 1, 1);
				(*eqfm_comb)[c3].parent_set_of_models=&dummy_parentsetofmodels;

				//event has focal mechanis; flags.full_field=0 indicates that an isotropic slip model should be used for all events (also those with foc mec):
				if (eqfm1[i].is_slipmodel && !flags.sources_all_iso) {
					err0 = focmec2slipmodel(crst, (*eqfm_comb)+c3, res, 1, 1);
					if (err0){
						print_screen("Error in creating slip model (function: eqkfm_addslipmodels)\n");
						print_logfile("Error in creating slip model (function: eqkfm_addslipmodels)\n");
						err+=1;
					}
					else {
						no_synthetic_slipmodels+=1;
					}
				}


				else{
					//event does not have foc. mech. but a fixed one should be used (flags.sources_without_focmec==2)
					if (!eqfm1[i].is_slipmodel && flags.sources_without_focmec==2){
						if ( crst.variable_fixmec==1){
							p=find_closest_point(crst.y, crst.x, crst.depth, crst.N_allP, eqfm1[i].y, eqfm1[i].x, eqfm1[i].depth);
							(*eqfm_comb)[c3].str1=crst.str0[p];	
							(*eqfm_comb)[c3].dip1=crst.dip0[p];
							(*eqfm_comb)[c3].rake1=crst.rake0[p];

	
						}
						else{
							(*eqfm_comb)[c3].str1=crst.str0[0];	
							(*eqfm_comb)[c3].dip1=crst.dip0[0];	
							(*eqfm_comb)[c3].rake1=crst.rake0[0];
						}
						(*eqfm_comb)[c3].whichfm=1;
						err0 = focmec2slipmodel(crst, (*eqfm_comb)+c3, res, 1, 1);
						if (err0){
							print_screen("Error in creating slip model (function: eqkfm_addslipmodels)\n");
							print_logfile("Error in creating slip model (function: eqkfm_addslipmodels)\n");
							err+=1;
						}
					}

					else{
						// If none of the conditions above holds, assume isotropic field.
						(*eqfm_comb)[c3].is_slipmodel=0;
						no_slipmodels+=1;
					}
				}

				(*nfout)[*Ncomb]=1;
				c3+=1;
				*Ncomb+=1;
			}
		}

		else {
			j=which_slipmod[i];
			nsm=0;
			for (int n=0; n<j; n++) nsm+=all_slipmodels.no_slipmodels[n];
			c2=0;
			print_screen("Using slip model %s for large event at t=%.5e, mag=%.2lf\n", all_slipmodels.slipmodels[nsm], eqfm1[i].t, eqfm1[i].mag);
			print_logfile("Using slip model %s for large event at t=%.5e, mag=%.2lf\n", all_slipmodels.slipmodels[nsm], eqfm1[i].t, eqfm1[i].mag);

			err += setup_eqkfm_element((*eqfm_comb)+c3, all_slipmodels.slipmodels+nsm, all_slipmodels.cmb_format,
									   all_slipmodels.no_slipmodels[j], crst.mu,
									   all_slipmodels.tmain[j],crst.N_allP, crst.list_allP,
									   all_slipmodels.mmain+j, all_slipmodels.cut_surf[j], NULL, crst.lat0, crst.lon0, all_slipmodels.constant_geometry);

			if (nfout)(*nfout)[*Ncomb]=nf2[which_slipmod[i]];
			for (int cc3=c3; cc3<c3+nf2[which_slipmod[i]]; cc3++) {
				(*eqfm_comb)[c3].distance=eqfm1[i].distance;
				(*eqfm_comb)[c3].index_cat=eqfm1[i].index_cat;
			}
			c3+=nf2[which_slipmod[i]];
			c_evfound+=1;
			*Ncomb+=1;
			all_slipmodels.is_used[j]=1;
		}
	}

	if (no_synthetic_slipmodels){
		print_screen("Using synthetic slip model from focal mechanism for %i earthquakes\n", no_synthetic_slipmodels);
		print_logfile("Using synthetic slip model from focal mechanism for %i earthquakes\n", no_synthetic_slipmodels);
	}

	if (no_slipmodels){
		print_screen("Using isotropic stress field for %i earthquakes\n", no_slipmodels);
		print_logfile("Using isotropic stress field for %i earthquakes\n", no_slipmodels);
	}

	if (c_evfound<N2){
		print_screen("** Warning: some slip models are not used, since events were not found in catalog (or are outside time boundaries)! (function: eqkfm_addslipmodels)\n");
		print_logfile("Warning: some slip models are not used, since events were not found in catalog (or are outside time boundaries)! (function: eqkfm_addslipmodels)\n");
	}

	return (err);
}

int focmec2slipmodel(struct crust crst, struct eqkfm *eqfm1, double res, int refine, int taper){
	/* Creates a synthetic slip model based on the focal mechanism and the empirical relations of Wells and Coppersmith.
	 *
	 * Input:
	 * 	eqfm1: description of the focal mechanism. Should already contain strX, dipX, rakeX, mag, whichfm.
	 * 	res: final resolution of the slip model (patch size, in km). Only used if refine==1.
	 * 	refine: flag indicating is the slip model should be refined.
	 * 	taper: flag indicating if the slip model should be tapered at the edges. Will taper on all 4 edges.
	 *
	 * Output:
	 *  eqfm1 is populated with the slip model.
	 *
	 */

	struct eqkfm eqkfm0;
	struct eqkfm *eqkfmP;
	int err=0;
	double slip;

	if (refine) {
		copy_eqkfm_all((*eqfm1), &eqkfm0);
		eqkfmP=&eqkfm0;
	}
	else eqkfmP=eqfm1;

	(*eqkfmP).is_slipmodel=1;
	(*eqkfmP).np_st=1;
	(*eqkfmP).np_di=1;
	(*eqkfmP).slip_str=darray(1,1);
	(*eqkfmP).slip_dip=darray(1,1);
	(*eqkfmP).open=NULL;
	(*eqkfmP).pos_s=darray(1,1);	//location of patches within fault; [0], [0] for single patch events.
	(*eqkfmP).pos_d=darray(1,1);
	(*eqkfmP).pos_s[1]=0;	//location of patches within fault; [0], [0] for single patch events.
	(*eqkfmP).pos_d[1]=0;

	WellsCoppersmith((*eqkfmP).mag, (*eqkfmP).rake1, &((*eqkfmP).L), &((*eqkfmP).W), &slip);
	slip=(*eqkfmP).tot_slip[0]=pow(10,(1.5*((*eqkfmP).mag+6)))*(1.0/(crst.mu*pow(10,12)*(*eqkfmP).W*(*eqkfmP).L));

	switch ((*eqkfmP).whichfm){
		case 1:
			if ((*eqkfmP).depth<0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1)) {
				if (extra_verbose) {
					print_screen("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
					print_logfile("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
				}
				(*eqkfmP).depth=0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1);
			}
			(*eqkfmP).slip_str[1]=slip*cos(DEG2RAD*(*eqkfmP).rake1);
			(*eqkfmP).slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfmP).rake1);
			break;
		case 2:
			if ((*eqkfmP).depth<0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip2)) {
				if (extra_verbose) {
					print_screen("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
					print_logfile("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
				}
				(*eqkfmP).depth=0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip2);
			}
			(*eqkfmP).slip_str[1]=slip*cos(DEG2RAD*(*eqkfmP).rake2);
			(*eqkfmP).slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfmP).rake2);
			break;
		case 0:
			if ((*eqkfmP).depth<0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1)) {
				if (extra_verbose) {
					print_screen("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
					print_logfile("Shifting center of mainshock fault from depth=%.2lf to depth=%.2lf to avoid fault being above surface (focmec2slipmodel). \n", (*eqkfmP).depth, 0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1));
				}
				(*eqkfmP).depth=0.5*(*eqkfmP).W*sin(DEG2RAD*(*eqkfmP).dip1);
			}
			if (extra_verbose) {
				print_screen("Warning: ambiguous nodal plane (focmec2slipmodel).\n");
				print_logfile("Warning: ambiguous nodal plane (focmec2slipmodel).\n");
			}
			(*eqkfmP).slip_str[1]=slip*cos(DEG2RAD*(*eqkfmP).rake1);
			(*eqkfmP).slip_dip[1]=-slip*sin(DEG2RAD*(*eqkfmP).rake1);
			break;
		default:
			print_screen("Error: illegal value for whichfm (focmec2slipmodel). Exiting. \n");
			print_logfile("Error: illegal value for whichfm (focmec2slipmodel). Exiting. \n");
			return(1);
	}

	if (refine) {
		err+=suomod1_resample(eqkfm0, eqfm1, res);	//create a slip model with right resolution.
		if (taper) err+=suomod1_taper((*eqfm1), eqfm1, 1, 1, 1, 1);
	}

	print_screen("Created slip model with %d patches for event [t, mag, lat, lon, dep]=[%.5e, %.2lf, %.3lf, %.3lf, %.3lf]\n", (*eqfm1).np_di*(*eqfm1).np_st, (*eqfm1).t, (*eqfm1).mag, (*eqfm1).lat, (*eqfm1).lon, (*eqfm1).depth);
	print_logfile("Created slip model with %d patches for event [t, mag, lat, lon, dep]=[%.5e, %.2lf, %.3lf, %.3lf, %.3lf]\n", (*eqfm1).np_di*(*eqfm1).np_st, (*eqfm1).t, (*eqfm1).mag, (*eqfm1).lat, (*eqfm1).lon, (*eqfm1).depth);

	return err;
}


//--------------------------read slip model files----------------------//

int read_eqkfm(char *fname, char *cmb_format, struct eqkfm **eqfm1, int *NF_out, double *Mw, double mu) {
/* Reads file containing a slip model.
 *
 * Input
 *  fname:	name of slip model file
 *  format: string containing the format of the slip model.
 *
 * Output
 *  *eqfm1: structure representing the slip model. Each element of the array is a subfault (i.e. rectangular area discretized into patches).
 *  			Range [0...*NF_out-1].
 * 	*Mw: earthquake magnitude calculated from the slip model
 * 	mu: chear modulus, used to calculate magnitude.
 */

	int NF, NP, err=0;
	double slip, dip_slip=0.0, strike_slip=0.0, M0=0.0;

	if (!(strcmp(cmb_format,"farfalle"))) err=read_farfalle_eqkfm(fname, eqfm1, &NF);
	else {
		if (!(strcmp(cmb_format,"pscmp"))) err=read_pscmp_eqkfm(fname, eqfm1, &NF);
		else{
			if (!(strcmp(cmb_format,"fsp"))) err=read_fsp_eqkfm(fname, eqfm1, &NF);
			else {
				print_screen(" ** Unknown slip model format %s (read_eqkfm).\n ", cmb_format);
				print_logfile("Unknown slip model format %s (read_eqkfm).\n", cmb_format);
				err=1;
			}
		}
	}

	if (err){
		print_screen(" ** Error: could not read slip model file %s (read_slipmodel). ** \n ", fname);
		print_logfile(" ** Error: could not read slip model file %s (read_slipmodel). ** \n ", fname);
		return (1);
	}

	if (eqfm1){
		for (int f=0; f<NF; f++){
			(*eqfm1)[f].str2=(*eqfm1)[f].str1;
			(*eqfm1)[f].dip2=(*eqfm1)[f].dip1;	//in this case there is no ambiguity (correct plane is known).
			(*eqfm1)[f].whichfm=1;
			(*eqfm1)[f].is_slipmodel=1;
			(*eqfm1)[f].tot_slip[0]=0.0;
			NP=(*eqfm1)[f].np_di*(*eqfm1)[f].np_st;
			for (int p=1; p<=NP; p++) {
				strike_slip+=(*eqfm1)[f].slip_str[p];
				dip_slip+=(*eqfm1)[f].slip_dip[p];
				slip=pow((*eqfm1)[f].slip_str[p]*(*eqfm1)[f].slip_str[p]+(*eqfm1)[f].slip_dip[p]*(*eqfm1)[f].slip_dip[p],0.5);
				(*eqfm1)[f].tot_slip[0]+=slip/NP;
			}
			(*eqfm1)[f].rake1=(*eqfm1)[f].rake2=(-RAD2DEG)*atan2(dip_slip, strike_slip);
			if (Mw) M0+=mu*(*eqfm1)[f].tot_slip[0]*(*eqfm1)[f].L*(*eqfm1)[f].W*1e12;
			(*eqfm1)[f].mag=(2.0/3.0)*log10(mu*(*eqfm1)[f].tot_slip[0]*(*eqfm1)[f].L*(*eqfm1)[f].W*1e12)-6.0;
		}
	}
	if (Mw) {
		if (!eqfm1){
			print_screen("*Warning: can not calculate Mw since eqfm1==NULL (read_eqkfm).*\n");
			print_logfile("Warning: can not calculate Mw since eqfm1==NULL (read_eqkfm).\n");
		}
		else *Mw=(2.0/3.0)*log10(M0)-6.0;
	}
	if (NF_out) *NF_out=NF;
	return (0);
}

int read_farfalle_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF_out) {
	/* Reads a slip model file in farfalle format.
	 * At the moment it does not support opening: pscmp format should be used instead.
	 *
	 * Input
	 *  fname:	name of slip model file
	 *
	 * Output
	 *  *eqfm_out: structure representing the slip model. Each element of the array is a subfault (i.e. rectangular area discretized into patches).
	 *  			Range [0...*NF_out-1].
	 *  			It can be passed as NULL, and will be ignored (will just count).
	 */

	// Variables used for MPI
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	FILE *fin;
	struct eqkfm *eqfm;
	int NF, NP, nchar=200;
	char line[nchar];
	double dlen, dwid;
	double *slips, *rakes;
	double plength, pwidth;
	long file_pos;

	if(procId == 0) {
		fin = fopen(fname,"r");
		if(fin == NULL) {
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		print_screen("Invalid input file (%s).\n", fname);
		print_logfile("Invalid input file passed to read_farfalle_eqkfm (%s).\n", fname);
		return (1);
	}
	else {
		if(procId == 0) {
			fgets(line,nchar,fin);
			fgets(line,nchar,fin);
			sscanf(line,"%d",&NF);
		}

		#ifdef _CRS_MPI
			MPI_Bcast(&NF, 1, MPI_INT, 0, MPI_COMM_WORLD);
		#endif

		if (eqfm_out) {
			eqfm=eqkfm_array(0,NF-1);
			for (int f=0; f<NF; f++){
				if(procId == 0) {
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%lf %lf %lf",&(eqfm[f].lat), &(eqfm[f].lon), &(eqfm[f].depth));
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%lf %lf",&(eqfm[f].str1), &(eqfm[f].dip1));
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%lf %lf",&(eqfm[f].L), &plength);
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%lf %lf",&(eqfm[f].W), &pwidth);
					fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					sscanf(line,"%d %d",&(eqfm[f].np_st), &(eqfm[f].np_di));
				}
				eqfm[f].whichfm=1;
			#ifdef _CRS_MPI
				MPI_Bcast(&(eqfm[f].lat), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].lon), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].depth), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].str1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].dip1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].L), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&plength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].W), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&pwidth, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].np_st), 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm[f].np_di), 1, MPI_INT, 0, MPI_COMM_WORLD);
			#endif
				NP=eqfm[f].np_st*eqfm[f].np_di;
				eqfm[f].pos_s=darray(1,NP);
				eqfm[f].pos_d=darray(1,NP);
				eqfm[f].slip_str=darray(1,NP);
				eqfm[f].slip_dip=darray(1,NP);
				eqfm[f].open=NULL;
				if ((f==0) | (NP>eqfm[f-1].np_st*eqfm[f-1].np_di)){
					if (f!=0){
						free_darray(slips, 1,1);
						free_darray(rakes, 1,1);
					}
					slips=darray(1,NP);
					rakes=darray(1,NP);
				}
				if(procId == 0) {
					fgets(line,nchar,fin);
					file_pos=ftell(fin);
				}
				dlen=eqfm[f].L/eqfm[f].np_st;
				dwid=eqfm[f].W/eqfm[f].np_di;
				for (int w=1; w<=eqfm[f].np_di; w++){
					for (int l=1; l<=eqfm[f].np_st; l++){
						if(procId == 0) {
							fscanf(fin, "%lf", slips+(w-1)*eqfm[f].np_st+l);
						}
					/* Abi's email "The hypocenter is the reference point. Then, you draw a line in the strike direction,
					 * and the distance between the hypocenter and the end of the fault is the partial length.
					 * The partial width is the distance to the top of the fault.
					 * It is the end that goes in the opposite direction of the strike, starting from the reference point, if I remember well,
					 * so that you always start drawing the fault in the direction of the strike.*/
					eqfm[f].pos_s[(w-1)*eqfm[f].np_st+l]=(l-0.5)*dlen-plength;
					eqfm[f].pos_d[(w-1)*eqfm[f].np_st+l]=(w-0.5)*dwid-pwidth;
					}
				}

				#ifdef _CRS_MPI
					MPI_Bcast(slips, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				#endif

				if(procId == 0) {
					fseek(fin,file_pos, SEEK_SET);		//since fscanf sometimes reads until next line, sometimes not (depending on small differences in the file, e.g. spaces).
					for (int w=1; w<=eqfm[f].np_di; w++) fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					for (int w=1; w<=eqfm[f].np_di; w++) fgets(line,nchar,fin);
					fgets(line,nchar,fin);
					file_pos=ftell(fin);
					for (int w=1; w<=eqfm[f].np_di; w++){
						for (int l=1; l<=eqfm[f].np_st; l++){
							fscanf(fin, "%lf", rakes+(w-1)*eqfm[f].np_st+l);
						}
					}
					fseek(fin,file_pos, SEEK_SET);		//since fscanf sometimes reads until next line, sometimes not (depending on small differences in the file, e.g. spaces).
					for (int w=1; w<=eqfm[f].np_di; w++) fgets(line,nchar,fin);
				}

				#ifdef _CRS_MPI
					MPI_Bcast(rakes, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				#endif

				for (int p=1; p<=NP; p++) {
					eqfm[f].slip_str[p]=slips[p]*cos(DEG2RAD*rakes[p]);
					eqfm[f].slip_dip[p]=-slips[p]*sin(DEG2RAD*rakes[p]);
				}
			}

			*eqfm_out=eqfm;
		}
		*NF_out=NF;
	}

	return(0);
}

int read_pscmp_eqkfm(char *fname, struct eqkfm **eqfm_out, int *NF2){
/* Reads a slip model file in pscmp format.
 *
 * Input
 *  fname:	name of slip model file
 *
 * Output
 *  *eqfm_out: structure representing the slip model. Each element of the array is a subfault (i.e. rectangular area discretized into patches).
 *  			Range [0...*NF2-1].
 *  			It can be passed as NULL, and will be ignored (will just count).
 */

	// [Fahad] Variables used for MPI.
	int fileError = 0;
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

	int err=0;
	FILE *fin;
	int NP, NF, djunk;
	int nchar=500;
	struct eqkfm *eqfm1;
	char comm[]="#";
	char line[nchar];

	if(procId == 0) {
		fin = fopen(fname,"r");
		if(fin == NULL) {
			fileError = 1;
		}
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&fileError, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if(fileError) {
		print_screen("Invalid input file (%s).\n", fname);
		print_logfile("Invalid input file passed to read_pscmp_eqkfm (%s).\n", fname);

		return (1);
	}

	line[0]=comm[0];
	if(procId == 0) {
		while (line[0]==comm[0])fgets(line, nchar, fin);
		for (int i=1; i<=3; i++){
			while (line[0]!=comm[0])fgets(line, nchar, fin);
			while (line[0]==comm[0]){
				fgets(line, nchar, fin);
			}
		}
		err+=(sscanf(line, "%d", &NF)!=1);
	}

	#ifdef _CRS_MPI
		MPI_Bcast(&NF, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	if (eqfm_out) {
		eqfm1=eqkfm_array(0,NF-1);
		line[0]=comm[0];
		if(procId == 0) {
			while (line[0]==comm[0])fgets(line, nchar, fin);
		}
		for (int f=0; f<NF; f++){
			if(procId == 0) {
				err+=(sscanf(line, "%d   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %d   %d   %lf",
						&djunk, &(eqfm1[f].lat), &(eqfm1[f].lon), &(eqfm1[f].depth), &(eqfm1[f].L),
						&(eqfm1[f].W), &(eqfm1[f].str1), &(eqfm1[f].dip1), &(eqfm1[f].np_st),
						&(eqfm1[f].np_di), &(eqfm1[f].t))!=11);
				fgets(line, nchar, fin);
			}

			#ifdef _CRS_MPI
				MPI_Bcast(&(eqfm1[f].lat), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].lon), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].depth), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].L), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].W), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].str1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].dip1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].np_st), 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].np_di), 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&(eqfm1[f].t), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			#endif

			eqfm1[f].str2=eqfm1[f].str1; eqfm1[f].dip2=eqfm1[f].dip1;	//in this case there is no ambiguity (correct plane is known).
			eqfm1[f].whichfm=1;
			NP=eqfm1[f].np_st*eqfm1[f].np_di;
			eqfm1[f].pos_s=darray(1,NP);
			eqfm1[f].pos_d=darray(1,NP);
			eqfm1[f].slip_str=darray(1,NP);
			eqfm1[f].slip_dip=darray(1,NP);
			eqfm1[f].open=darray(1,NP);
			if(procId == 0) {
				for (int p=1; p<=NP; p++) {
					err+=(sscanf(line, "%lf   %lf    %lf   %lf   %lf",
							&(eqfm1[f].pos_s[p]), &(eqfm1[f].pos_d[p]),
							&(eqfm1[f].slip_str[p]), &(eqfm1[f].slip_dip[p]),
							&(eqfm1[f].open[p]))!=5);
					if (!feof(fin)) fgets(line, nchar, fin);
				}
			}

			#ifdef _CRS_MPI
				MPI_Bcast(eqfm1[f].pos_s, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(eqfm1[f].pos_d, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(eqfm1[f].slip_str, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(eqfm1[f].slip_dip, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(eqfm1[f].open, NP+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			#endif
		}
	}

	if(procId == 0) {
		fclose(fin);
	}

	if (NF2!= NULL) *NF2=NF;
	if (eqfm_out) *eqfm_out=eqfm1;
	if (err) {
                print_screen("**Couldn't read slip model file: %s. \n", fname);
                print_logfile("**Couldn't read slip model file: %s. \n", fname);
	}

	return err;
}
