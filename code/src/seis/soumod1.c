
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


/* Functions written by Katrin Kieling (2010) and Camilla Cattania (2013)
 */

#include "soumod1.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int scale_to_mag(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double * slips, double *rakes){
	/* Rescales the slip from eqkfm2 so that they have the same average slip (i.e. magnitude)
	 * It assumes that the two models have the same area, but different discretizations.
	 *
	 * Input:
	 *  eqkfm1: contains input slip model.
	 *  eqkfm2: contains no. of patches.
	 *  slips, rakes: contains slip values (scalar) and their rakes for eqkfm2
	 *
	 * Output:
	 *  (*eqkfm2).slip_str, (*eqkfm2).slip_dip are populated.
	 */

	double slip, M0old, M0;	//slip, seismic moments.
	int N1= eqkfm1.np_di*eqkfm1.np_st;
	int N2= (*eqkfm2).np_di*(*eqkfm2).np_st;

	M0old=0.0;
	for (int p=1; p<=N1; p++) {
		slip=pow(eqkfm1.slip_str[p]*eqkfm1.slip_str[p]+eqkfm1.slip_dip[p]*eqkfm1.slip_dip[p],0.5);
		M0old+=slip/N1;
	}

	M0=0.0;
	for (int p=1; p<=N2; p++) M0+=slips[p]/N2;

	for (int p=1; p<=N2; p++) {
		slip=slips[p]*M0old/M0;
		(*eqkfm2).slip_str[p]=slip*cos((rakes[p])*pi/180);
		(*eqkfm2).slip_dip[p]=-slip*sin((rakes[p])*pi/180);
		if (isnan((*eqkfm2).slip_dip[p]) || isnan((*eqkfm2).slip_str[p])) return 1;
	}

	return 0;
}

int suomod1_resample(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, double disc){

/* Resamples a slip model to increase resolution.
 *
 * Input:
 *  eqkfm1: original slip model
 *  disc: required resolution
 *
 * Output:
 *  eqkfm2: refined slip model.
 */

  char     outname2[]="source_res.out";
  double   reflat, reflon, refdepth, strike, rake, dip1;
  double   odiscx, odiscy, maxox=-1e10, maxoy=-1e10, minox=1e10, minoy=1e10, start_x, start_y;
  double   *ox, *oy, *REslipo, *rake_v, *strike_v, *dip;
  double   *nx, *ny, *REold, *IMold, *dipold, *rakeold, *strikeold;
  double   oxdim=0.0 /*, oxdim_tot*/;
  double   oydim=0.0;
  double 	*slip;
  double 	ndiscx, ndiscy;
  int	    print_Fourier, printout;
  int      i, j, k, l;
  int      ns, nns, nsx, nsy;

  //can be set manually to print out extra files.
  printout= 0;
  print_Fourier= 0;

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

    ns=eqkfm1.np_st*eqkfm1.np_di;
    REslipo= darray(1,ns);
    rake_v= darray(1,ns);
    strike_v= darray(1,ns);
    dip= darray(1,ns);
	reflat=eqkfm1.lat;
	reflon=eqkfm1.lon;
	refdepth=eqkfm1.depth;
	switch (eqkfm1.whichfm){
		case 1:
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip1=eqkfm1.dip1;
			break;
		case 2:
			rake=eqkfm1.rake2;
			strike=eqkfm1.str2;
			dip1=eqkfm1.dip2;
			break;
		case 0:
			print_screen("Warning: ambiguous focal plane for eqkfm1 (suomod1_resample) -> using first plane!\n");
			print_logfile("Warning: ambiguous focal plane for eqkfm1 (suomod1_resample) -> using first plane!\n");
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip1=eqkfm1.dip1;
			break;
		default:
			print_screen("Error: eqkfm1.whichfm has illegal value in suomod1_resample(%d)!\n",eqkfm1.whichfm);
			print_logfile("Error: eqkfm1.whichfm has illegal value in suomod1_resample(%d)!\n",eqkfm1.whichfm);
			return(1);
	}
	oxdim=eqkfm1.L;
	oydim=eqkfm1.W;
	nsx=eqkfm1.np_st;
	nsy=eqkfm1.np_di;
	ox=eqkfm1.pos_s;
	oy=eqkfm1.pos_d;
	odiscx=oxdim/nsx;
	odiscy=oydim/nsy;

	#pragma omp parallel for
	for (int np=1; np<=ns; np++){
		REslipo[np]=pow(eqkfm1.slip_str[np]*eqkfm1.slip_str[np]+eqkfm1.slip_dip[np]*eqkfm1.slip_dip[np],0.5);
		rake_v[np]=(-180/pi)*atan2(eqkfm1.slip_dip[np],eqkfm1.slip_str[np]);	//check this! (sign)
		if (isnan(rake_v[np])==1) rake_v[np]=rake;	//for example if slip=0 in both directions.
		strike_v[np]=strike;
		dip[np]=dip1;
	}

	for (i=1;i<=ns;i++){
	   if (ox[i]>maxox) maxox=ox[i];
	   if (oy[i]>maxoy) maxoy=oy[i];
	   if (ox[i]<minox) minox=ox[i];
	   if (oy[i]<minoy) minoy=oy[i];
	}

	//oxdim_tot=(maxox-minox)+odiscx;	//could be larger than oxdim if fault is skewed. Is this needed?

	 maxoy=maxoy+odiscy;
	 //resampling of slipmap

	  nns=(int)(ceil(oxdim/disc)*ceil(oydim/disc));
	  nsx=(int)ceil(oxdim/disc);
	  nsy=(int)ceil(oydim/disc);
	  ndiscx=oxdim/nsx;
	  ndiscy=oydim/nsy;
	  start_x=minox-0.5*odiscx+0.5*ndiscx;
	  start_y=minoy-0.5*odiscy+0.5*ndiscy;
	  REold    = darray(1,nns);
	  IMold    = darray(1,nns);
	  dipold   = darray(1,nns);
	  rakeold   = darray(1,nns);
	  strikeold   = darray(1,nns);
	  nx       = darray(1,nns);
	  ny       = darray(1,nns);
      slip=darray(1,nns);

      	  i=1;
	  for (k=1; k<=nsy; k++)
	        { for (j=1;j<=nsx;j++)
	                { nx[i]=start_x+(j-1)*ndiscx; //relative position along strike (note: new fault not skewed).
	                  ny[i]=start_y+(k-1)*ndiscy; //relative position along dip
	                  i++;
	                }
	        }

	  if (ndiscx<odiscx && ndiscy<odiscy) {// refine slip distribution
          for (int i=1; i<=nns; i++){
        	  l=1;
        	  while (l<ns && (fabs(ox[l]-nx[i])>odiscx || fabs(oy[l]-ny[i])>odiscy)) l++;
        	  if (l==ns && (fabs(ox[l]-nx[i])>odiscx || fabs(oy[l]-ny[i])>odiscy)) {	//new point is outside old fault (possible for skewed fault).
            	  REold[i]=0.0;
    			  IMold[i]=0;
    			  dipold[i]=dip1;
    			  rakeold[i]=rake;
    			  strikeold[i]=strike;
        	  }
        	  else {
				  REold[i]=REslipo[l];	//new point is inside new fault
				  IMold[i]=0;
				  dipold[i]=dip[l];
				  rakeold[i]=rake_v[l];
				  strikeold[i]=strike_v[l];
        	  }
          }
      }

	  else if (ndiscx>odiscx && ndiscy>odiscy){
		  if (extra_verbose) {
			  print_screen("New resolution is larger than old one: will not resample (suomod1_resample).\n");
			  print_logfile("New resolution is larger than old one: will not resample (suomod1_resample).\n");
		  }
		  copy_eqkfm_all(eqkfm1, eqkfm2);
		  return (0);
	  }

	  else  {
		  if (extra_verbose) {
			  print_screen("Model has right discretization - will not be resampled. \n");
		  }
		  copy_eqkfm_all(eqkfm1, eqkfm2);
		  return(0);
	  }

  //--------------Fill in eqkfm2.-------------------//

    copy_eqkfm_noslipmodel(eqkfm1, eqkfm2);

 	(*eqkfm2).L=oxdim;
 	(*eqkfm2).W=oydim;
	(*eqkfm2).np_st=nsx;
	(*eqkfm2).np_di=nsy;
	(*eqkfm2).pos_s=darray(1,nns);
	(*eqkfm2).pos_d=darray(1,nns);
	(*eqkfm2).slip_str=darray(1,nns);
	(*eqkfm2).slip_dip=darray(1,nns);

	for (int p=1; p<=nns; p++) {
		(*eqkfm2).pos_s[p]=nx[p];
		(*eqkfm2).pos_d[p]=ny[p];
	}

	scale_to_mag(eqkfm1, eqkfm2, REold, rakeold);		//final slip.
//	(*eqkfm2).tot_slip[0]=tot_slip(*eqkfm2);


  //--------------Print things out-------------------//

  if (printout==1) {
	  if(procId == 0) {
		  print_slipmodel(outname2,eqkfm2,1);
	  }
  }

  //----------------Free memory----------------------//

  free_darray(REold,1,nns);
  free_darray(REslipo,1,ns);
  free_darray(rake_v,1,ns);
  free_darray(strike_v,1,ns);
  free_darray(dip,1,ns);
  free_darray(IMold,1,nns);
  free_darray(strikeold,1,nns);
  free_darray(dipold,1,nns);
  free_darray(rakeold,1,nns);
  free_darray(slip,1,nns);
  free_darray(nx,1,nns);
  free_darray(ny,1,nns);
//  if (ndiscx>odiscx && ndiscy>odiscy){
//	  free_darray(lat,1,nns);
//	  free_darray(lon,1,nns);
//	  free_darray(z,1,nns);
//  }

  return(0);
}

int suomod1_taper(struct eqkfm eqkfm1, struct eqkfm *eqkfm2, int top, int bottom, int right, int left){
/* Tapers a slip model (keeping same geometry/discretization).
 *
 * Input
 *  eqkfm1: original slip model.
 *  top, bottom, right, left: flags indicating which sides should be tapered.
 *
 * Output:
 *  eqkfm2: tapared slip model.
 */

  char     outname2[200];
  double   strike, rake, dip;
  double   odiscx, odiscy, maxox=-1e10, maxoy=-1e10, minox=1e10, minoy=1e10, *start_x, start_y;
  double   *ox, *oy, *REslipo, *slip, *rakes;
  double   *nx, *ny;
  double 	*taper;
  double   alphax, alphay;
  double   oxdim=0.0;
  double   oydim=0.0;
  int	   printout;
  int      ns, nns, nsx, nsy;

    // can be set to 1 to print extra files.
	printout=0;

	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

    ns=nns=eqkfm1.np_st*eqkfm1.np_di;
    if (eqkfm1.np_st==1) left=right=0;
    if (eqkfm1.np_di==1) top=bottom=0;
    if (ns==1) {
    	copy_eqkfm_all(eqkfm1, eqkfm2);	//bug fix: using copy_eqkfm_slipmodel did not copy nsel, so no points were selected.
    	if (extra_verbose) {
    		print_screen("** Warning: model has a single patch, will not be tapered.**\n");
    	   	print_logfile("** Warning: model has a single patch, will not be tapered.**\n");
    	}
    	return 0;
    }

	REslipo=darray(1,ns);
	rakes=darray(1,ns);
	taper=darray(1,ns);
	start_x = darray(1,ns);
	nx = darray(1,nns);
	ny = darray(1,nns);
	slip = darray(1,nns);

    if (printout) sprintf(outname2, "source_tap%d%d%d%d.out",top, bottom, right, left);

	switch (eqkfm1.whichfm){
		case 1:
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip=eqkfm1.dip1;
			break;
		case 2:
			rake=eqkfm1.rake2;
			strike=eqkfm1.str2;
			dip=eqkfm1.dip2;
			break;
		case 0:
			print_screen("Warning: ambiguous focal plane for eqkfm1 (suomod1_taper) -> using first plane!\n");
			print_logfile("Warning: ambiguous focal plane for eqkfm1 (suomod1_taper) -> using first plane!\n");
			rake=eqkfm1.rake1;
			strike=eqkfm1.str1;
			dip=eqkfm1.dip1;
			break;
		default:
			print_screen("Error: eqkfm1.whichfm has illegal value in suomod1_taper(%d)!\n",eqkfm1.whichfm);
			print_logfile("Error: eqkfm1.whichfm has illegal value in suomod1_taper(%d)!\n",eqkfm1.whichfm);
			return(1);
			break;
	}
	oxdim=eqkfm1.L;
	oydim=eqkfm1.W;
	nsx=eqkfm1.np_st;
	nsy=eqkfm1.np_di;
	ox=eqkfm1.pos_s;
	oy=eqkfm1.pos_d;
	odiscx=oxdim/nsx;
	odiscy=oydim/nsy;

	#pragma omp parallel for
	for (int np=1; np<=ns; np++) {
		REslipo[np]=pow(eqkfm1.slip_str[np]*eqkfm1.slip_str[np]+eqkfm1.slip_dip[np]*eqkfm1.slip_dip[np],0.5);
		rakes[np]=(-180/pi)*atan2(eqkfm1.slip_dip[np],eqkfm1.slip_str[np]);
	}

	for (int i=1;i<=ns;i++){
	   if (ox[i]>maxox) maxox=ox[i];
	   if (oy[i]>maxoy) {
		   maxoy=oy[i];
		   start_x[i]=ox[i];	//new row of patches -> new starting point (used for tapering).
	   }
	   else start_x[i]=start_x[i-1];
	   if (ox[i]<minox) minox=ox[i];
	   if (oy[i]<minoy) minoy=oy[i];
	}
	start_y=minoy;	//assume this is the same for all columns (unlike for rows).
	for (int i=1; i<=ns; i++){
	  nx[i]=ox[i]-start_x[i];
	  ny[i]=oy[i]-start_y;
	}

	//-------Tapering-----------//
	/* create taper (similar to tukey window but with onlay 1/4 of the cos period, alpha = fraction of length used for decay, alpha=1 = rectangular window, alpha=0 = similar Hann window)
	 */

	alphax=0.51;
	alphay=0.51;

	#pragma omp parallel for
	for (int i=1;i<=nns;i++){
		taper[i] = ((right==1 && (nx[i]-(nsx-1)/2*odiscx)>alphax*nsx/2*odiscx) || (left==1 && -(nx[i]-(nsx-1)/2*odiscx)>alphax*nsx/2*odiscx))?
			  cos(pi/2*(fabs(fabs(nx[i]-(nsx-1)/2*odiscx)-alphax*nsx/2*odiscx))/((1-alphax)*nsx/2*odiscx)) : 1;
		if ((bottom==1 && (ny[i]-(nsy-1)/2*odiscy)>alphay*nsy/2*odiscy) || (top==1 && -(ny[i]-(nsy-1)/2*odiscy)>alphay*nsy/2*odiscy)){
			  taper[i] *=(cos(pi/2*(fabs(fabs(ny[i]-(nsy-1)/2*odiscy)-alphay*nsy/2*odiscy))/((1-alphay)*nsy/2*odiscy)));
		}
	}

  //Fill in eqkfm2.
  if (eqkfm2!=&eqkfm1){

	copy_eqkfm_noslipmodel(eqkfm1, eqkfm2);

	(*eqkfm2).L=oxdim;
	(*eqkfm2).W=oydim;
	(*eqkfm2).np_st=nsx;
	(*eqkfm2).np_di=nsy;
	(*eqkfm2).pos_s=eqkfm1.pos_s;
	(*eqkfm2).pos_d=eqkfm1.pos_d;
	if ((*eqkfm2).slip_str==NULL) (*eqkfm2).slip_str=darray(1,nns);
	if ((*eqkfm2).slip_dip==NULL) (*eqkfm2).slip_dip=darray(1,nns);

  }

	#pragma omp parallel for
	for (int p=1; p<=nns; p++) {
		slip[p]=REslipo[p]*taper[p];		//final slip.
	}

	scale_to_mag(eqkfm1, eqkfm2, slip, rakes);		//final slip.
//	(*eqkfm2).tot_slip[0]=tot_slip(*eqkfm2);

  //--------------Print things out-------------------//


	if (printout==1) {
		if(procId == 0) {
			print_slipmodel(outname2,eqkfm2,1);
		}
	}


  //--------------Free memory-------------------//

  free_darray(slip,1,ns);
  free_darray(REslipo,1,ns);
  free_darray(start_x,1,ns);
  free_darray(nx,1,ns);
  free_darray(ny,1,ns);
  free_darray(taper,1,ns);
  free_darray(rakes,1,ns);
  return 0;

}
