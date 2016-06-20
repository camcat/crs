memory_error_quit;
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


#include "forecast_stepG.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int rate_state_evolution(struct catalog cat, double *times, double **cmpdata, struct pscmp *DCFS, double tt0, double tt1, double dt_step, double Asig, double ta,
			int points[], double *out_NeX, double *NeT, double *ReT, double **all_NeT, int N, int NTS, int Neqks, double *gamma_init, double *back_rate, double *R, int last){

	/* Calculates seismicity in space and time based on rate-and-state (Dieterich 1994) constitutive law.
	 * Aseismic stresses in cmpdata are assumed to grow linearly between time steps.
	 *
	 * Input:
	 *  cat: catalog containing earthquakes to be included in LogLikelihood calculation.
	 *  times: timesteps corresponding to aseismic loading. [0...NTS]
	 *  cmpdata: stresses due to aseismic loading. [0...NTS-1; 1...N], where cmpdata[x][m] is the stress change at grid point [m] between times[x+1] and times[x].
	 *  		if NULL, it will be ignored.
	 *  DCFS: structure containing seismic sources. [0...Neqks-1]
	 *  tt0, tt1: start, end time of calculation
	 *  Asig, ta: rate-state parameters
	 *  points: list of indices of points to be included (referred to grid vectors in crst structure in main). if NULL, use sequence [1,2,3,...,N]
	 *  N, NTS, Neqks: no. of grid points, time steps for aseismic sources, seismic sources.
	 *  gamma_init: values of gammas at time tt0. NB: indices refer to *entire grid* (not just those elements of "points"): [1...NgridT]
	 *  back_rate: spatially nonuniform background rate. array containing NgridT times the ratio between avg rate and rate at that point. (so that the sum of back_rate is NgridT).
	 *   		if NULL, will assume it's 1 for all points (uniform rate). NB: indices refer to *entire grid* (not just those elements of "points"): [1...NgridT]
	 *  last: flag indicating if the values of gamma_init should be overwritten at the end, with value at tt1.
	 *
	 * Output:
	 *  out_NeX: number of events in each grid cell. [1...N] (i.e. it refers to the subset of points in "points").
	 *  NeT: total no of events in each time step (i.e., sum of out_NeX). [0...nts-1]
	 *  ReT: total seismicity rate at the end of each time step.  [0...nts-1]
	 *  all_NeT: no of events binned in space and time  [0...nts-1][1...N];
	 *  R: seismicity rate at the time of the earthquakes given in cat (used for first term of LogLikelihood).
	 *
	 *  NB: assumes that times, cat, DCFS, NTS, times are always the same (each time function is called), since many static internal structures are initialized in the first function call.
	 */


	// Variables used for MPI.
	int procId = 0;

	#ifdef _CRS_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	#endif

  double  tau, dtau_dt, dtau_dt00=Asig/ta, ta1;
  int     TS0, TS1, n;
  double  gamma, back_rate_n;
  double  tt0step, tt1step;

  //events will contain times, magnitudes of all earthquakes (both from catalog and from DCFS):
  static double **events;
  static int **indices;	//contains indices of events referred to DCFS and cat.
  int Neq;	//size of events[0...Neq];
  int is_incat, is_inDCFS;	//flag for events in events
  int cat_i, DCFS_i; //indices of elements in events referred to cat, DCFS.
  static int *TS_eqk;	//time step before each earthquake in events.

  // variables used to find earthquakes within each grid point (see below):
  static int **which_eqk;
  static int *num_eqk;
  static float ** cat_weights;
  static int ** DCFS_whichpt;

  //internal variables for seismicity calculation:
  double *NeX;
  double **Rprivate;	//private variable for each OMP thread, corresponds to global R.
  double **NeTprivate, **ReTprivate;
  static double *dt, *NeXdum, *ReX;


  //dummy variables:
  int err=0, errtot=0;
  int j0, next_eqk, counter_eqk, next_TS;
  int counter1, counter2, k;
  double t_pre, t_now, t_endstep;
  double a,b;
  int reach_end;
  int step, nts;
  int nthreads, nthreadstot=omp_get_max_threads();

  static int firsttimein=1;
  int warning_printed=0;


  if (firsttimein==1){
		firsttimein=0;

		//allocate memory:
		NeXdum=darray(1,N);
		ReX=darray(1,N);

		//----------------------------------------------------------------------------------------------//
		//					Set up list of earthquakes for each grid point								//
		//----------------------------------------------------------------------------------------------//

		//create a merged set of events from cat and DCFS; this is needed since the both define calculation time steps.
		events=union_cats2(cat, DCFS, Neqks, &indices, &Neq);

		which_eqk=i2array(1,N,0,Neq-1);	//list of earthquakes in each grid point
		num_eqk=iarray(1,N);	//number of earthquakes in each grid point
		cat_weights=f2array(1,N,0,Neq-1);	//catalog weights
		DCFS_whichpt=i2array(1,N,0,Neq-1);	//grid point index in DCFS elements (since they contain stress field for a selection of pts only)

		for (int n=1; n<=N; n++) num_eqk[n]=0;

		for(int eq=0;eq<Neq;eq++){
			counter1=counter2=1;
			cat_i=indices[1][eq];	//	NB:cat_i[i]==-1 means no events selected.
			DCFS_i=indices[2][eq];  //	NB:DCFS_i[i]==-1 means no events selected.
			for (int n=1; n<=N; n++){
				if (counter2 > DCFS[DCFS_i].nsel && counter1 > cat.ngrid[cat_i]) break;	//all points have already been included
				else{
					//find index in cat and DCFS which may refer to this grid point:
					if (cat_i!=-1) while (counter1 <=cat.ngrid[cat_i] && cat.ngridpoints[cat_i][counter1]<n) counter1+=1;
					if (DCFS_i!=-1) while (counter2 <=DCFS[DCFS_i].nsel && DCFS[DCFS_i].which_pts[counter2]<n) counter2+=1;
					//check if index refers to this grid point:
					is_incat= (cat_i!=-1 && counter1 <= cat.ngrid[cat_i] && cat.ngridpoints[cat_i][counter1]==n) ? 1 : 0;
					is_inDCFS= (DCFS_i!=-1 && counter2 <= DCFS[DCFS_i].nsel && DCFS[DCFS_i].which_pts[counter2]==n) ? 1 : 0;
					//add earthquake to list for this grid point, and fill in cat_weights, DCFS_whichpt accordingly
					if (is_incat | is_inDCFS){
						which_eqk[n][num_eqk[n]]=eq;
						cat_weights[n][num_eqk[n]]= (is_incat)? cat.weights[cat_i][counter1] : 0.0;	//only in DCFS
						DCFS_whichpt[n][num_eqk[n]]= (is_inDCFS)? counter2 : 0;	//only in cat: will not change gamma later.
						num_eqk[n]+=1;
					}
				}
			}
		}

		//----------------------------------------------------------------------------------------------//
		//									Time steps initialization									//
		//----------------------------------------------------------------------------------------------//

		//calculate time step duration:
		if (NTS!=0){
			dt=darray(0,NTS-1);
			for (int g=0; g<NTS; g++) dt[g]=times[g+1]-times[g];
		}

		//find last time step before each earthquake:
		TS_eqk=iarray(0,Neq-1);
		for (int i=0; i<Neq; i++){
			if (NTS==0) TS_eqk[i]=1;	//this will make t_pre be correct later (j0<NTS condition).
			else{
				k=0;
				while(k<NTS && times[k]<events[1][i]) k++;
				if(times[k]>=events[1][i]) k--;
				TS_eqk[i]=k;
				if (k<0 && events[1][i]>=tt0) {
					print_screen("**Warning: no time steps available before earthquake no. %d ** (forecast_stepG2_new.c)\n",i);
					print_logfile("**Warning: no time steps available before earthquake no. %d ** (forecast_stepG2_new.c)\n",i);
					return 1;
				}
			}
		}
  }

  if (Asig==0 && ta==0.0) return(0);	//in this case, function has only been called to setup variables above.

  NeX= (out_NeX)? out_NeX : NeXdum;

  // Initialize vector to 0;
  if (ReT) for(int m=1;m<=N;m++) ReX[m]=0.0;
  if (NeX) for(int m=1;m<=N;m++) NeX[m]=0.0;

  if (tt1<tt0) {
	  print_screen("\n*** Warning: tt1<tt0 in forecast_stepG2_new.c  ***\n");
	  print_logfile("\n*** Warning: tt1<tt0 in forecast_stepG2_new.c  ***\n");
	  return 1;
  }
  if (NTS>0 && times[0]>tt0 && times[NTS]<tt1) {
	  print_screen("\n** Warning: time steps in forecast_stepG don't cover entire forecast range!**\n");
	  print_logfile("\n** Warning: time steps in forecast_stepG don't cover entire forecast range!**\n");
	  return 1;
  }

  //TS0= First time step after tt0.
  TS0=0;
  if (NTS==0) TS0=1;	//this will make t_pre be correct later (j0<NTS condition).
  else{
	  while(TS0<NTS && times[TS0]<=tt0) TS0++;
	  if(procId == 0) {
		  if(times[TS0]<=tt0 & !warning_printed) {
			  warning_printed=1;
			  print_screen("\n*** Warning: times[TS0]<=tt0 in forecast_stepG2_new.c  ***\n");
			  print_logfile("\n*** Warning: times[TS0]<=tt0 in forecast_stepG2_new.c  ***\n");
			  return 1;
		  }
	  }
  }

  //TS1= Last time step before tt1.
  TS1=1;	//this will make t_pre be correct later (j0<NTS condition).
  if (NTS>0){
	  while(TS1<NTS-1 && times[TS1]<tt1) TS1++;
	  if(times[TS1]>=tt1) TS1--;
  }

  //Rprivate is used for parallelization (a barrier creates a bottleneck and poor performance):
  Rprivate=d2array(0,nthreadstot-1, 0, cat.Z);
  if (!Rprivate) memory_error_quit;
  for (int t=0; t<omp_get_max_threads(); t++){
	  for (int eq=0; eq<=cat.Z; eq++) Rprivate[t][eq]=0.0;
  }

  //Similar to above, but for individual time steps:
  nts=(dt_step<tol0) ? 0 : ceil((tt1-tt0)/dt_step);	//number of output time steps;
  for (int i=0; i<nts; i++){
    if (ReT) ReT[i]=0.0;
    if (NeT) NeT[i]=0.0;
  }

  NeTprivate=d2array(0,nthreadstot-1, 0, nts-1);
  ReTprivate=d2array(0,nthreadstot-1, 0, nts-1);
  if (!NeTprivate | !ReTprivate) memory_error_quit;
  for (int t=0; t<omp_get_max_threads(); t++){
	  if (NeT) for (int i=0; i<nts; i++) NeTprivate[t][i]=0.0;
	  if (ReT) for (int i=0; i<nts; i++) ReTprivate[t][i]=0.0;
  }

  err=0;

  //loop over grid points:
  #pragma omp parallel for firstprivate(err) private(n, gamma,tau,j0, next_eqk, next_TS, t_now, t_pre, reach_end, t_endstep, cat_i, DCFS_i, a, b, dtau_dt, counter_eqk, back_rate_n, ta1, tt0step, tt1step, step, TS0, TS1) reduction(+:errtot)
  for(int m=1;m<=N;m++){

	nthreads=omp_get_num_threads();

    tt0step=tt0;
    tt1step= (nts==1)? tt1 : tt0step+dt_step;	//to make sure it enters at least once (floating point error may give problem otherwise).
   
	t_now=tt0;

	if (err!=0) continue;
	n=(points==0)? m : points[m];

	//set background rate and starting gamma given as arguments:
	back_rate_n= (back_rate) ? back_rate[n] : 1.0;
	gamma=gamma_init[n];
	if (NeX) NeX[m]=0.0;

	step=0;

	while (tt1step<=tt1){	//overall time span

		if (dt_step<tol0) break;	//to avoid infinite loop if tt0=tt1 and dt_step=0;

		//TS0= First time step after tt0.
		TS0=0;
		if (NTS==0) TS0=1;	//this will make t_pre be correct later (j0<NTS condition).
		else{
		  while(TS0<NTS && times[TS0]<=tt0step) TS0++;
		  if(procId == 0) {
			  if(times[TS0]<=tt0 & !warning_printed) {
				  warning_printed=1;
				  print_screen("\n*** Warning: times[TS0]<=tt0 in forecast_stepG2_new.c  ***\n");
				  print_logfile("\n*** Warning: times[TS0]<=tt0 in forecast_stepG2_new.c  ***\n");
			  }
		  }
		}

		//TS1= Last time step before tt1.
		TS1=1;	//this will make t_pre be correct later (j0<NTS condition).
		if (NTS>0){
		  while(TS1<NTS-1 && times[TS1]<tt1step) TS1++;
		  if(times[TS1]>=tt1step) TS1--;
		}

		j0=TS0;     //next time step (first time step after tt0)
	    counter_eqk=0;
		reach_end=0;

		//find first earthquake after or at tt0:
		while (counter_eqk<num_eqk[n] && events[1][which_eqk[n][counter_eqk]]<tt0step) counter_eqk+=1;

		//loop over time steps until reaching tt1step:
		while (tol0<(tt1step-t_now)){

			if (err!=0) break;	//error in a previous loop;
			if (counter_eqk>=num_eqk[n]) reach_end=1;	//all earthquakes for this grid point have been included (counter_eqk is out of range).
			else {
				next_eqk=which_eqk[n][counter_eqk];
				if (events[1][next_eqk]>=tt1step) reach_end=1;
			}

			t_endstep= reach_end ? tt1step : events[1][next_eqk];	//end time of this while loop iteration;
			next_TS= reach_end ? TS1 : TS_eqk[next_eqk];		//last time step before t_endstep
			cat_i= reach_end ? 0 : indices[1][next_eqk];		//index of next event (occurring at time events[1][next_eqk])
			DCFS_i= reach_end ? 0 : indices[2][next_eqk];		//index of next event (occurring at time events[1][next_eqk])

			// evolve seismicity up to next barrier (which is the smallest between next earthquake (events[1][next_eqk]), next time step (times[j0]) or tt1step)
			t_pre= (j0<=NTS) ? fmin(t_endstep, times[j0]) - t_now : t_endstep - t_now;	//find time left to next barrier

			dtau_dt=(cmpdata && j0-1>=0 && j0<=NTS)? (Asig/ta)+cmpdata[j0-1][n]/dt[j0-1] : (Asig/ta);//stressing rate during current time step (including stress step from cmpdata and background stressing rate):
			if (t_pre>tol0){
				tau=dtau_dt*t_pre;	//stress change during current time step
				ta1=Asig/dtau_dt;	//dummy variable
				if (NeX) {
					//find no. of earthquakes
					a=gamma*dtau_dt-1;	//dummy variable
					b=a*exp(-t_pre/ta1)+1;	//dummy variable
					if (!isinf(fabs(b))) NeX[m]+= fmax(0.0, back_rate_n*(dtau_dt/dtau_dt00)*(t_pre+ta1*log(b/(gamma*dtau_dt))));	//due to numerical error it can give -ve values.
				}
				//update gamma:
				gamma=(fabs(tau/Asig)>1e-10)? (gamma-t_pre/(tau))*exp(-tau/Asig)+t_pre/(tau) : gamma*(1-tau/Asig)+t_pre/Asig;
			}

			//evolve seismicity between time steps:
			t_now+=t_pre;

			for(int j=j0;j<next_TS;j++){
				dtau_dt=(cmpdata)? (Asig/ta)+cmpdata[j][n]/dt[j] : (Asig/ta);
				tau=dtau_dt*dt[j];
				ta1=Asig/dtau_dt;
				if (NeX) {
					a=gamma*dtau_dt-1;
					b=a*exp(-dt[j]/ta1)+1;

					if (!isinf(fabs(b))) NeX[m]+=fmax(0.0, back_rate_n*(dtau_dt/dtau_dt00)*(dt[j]+ta1*log(b/(gamma*dtau_dt))));	//due to numerical error it can give -ve values for gamma->inf
				}
				gamma=(fabs(tau/Asig)>1e-10)? (gamma-dt[j]/(tau))*exp(-tau/Asig)+dt[j]/(tau) : gamma*(1-tau/Asig)+dt[j]/Asig;
				t_now+=dt[j];
			}

			// evolve seismicity till the end:
			t_pre=t_endstep-t_now;

			if (t_pre>tol0) {
				dtau_dt=(cmpdata)? (Asig/ta)+(1.0/dt[next_TS])*cmpdata[next_TS][n] : (Asig/ta);
				tau=dtau_dt*t_pre;
				ta1=Asig/dtau_dt;
				if (NeX) {
					a=gamma*dtau_dt-1;
					b=a*exp(-t_pre/ta1)+1;
					if (!isinf(fabs(b))) NeX[m]+= fmax(0.0, back_rate_n*(dtau_dt/dtau_dt00)*(t_pre+ta1*log(b/(gamma*dtau_dt)))); //due to numerical error it can give -ve values for gamma->inf
				}
				gamma=(fabs(tau/Asig)>1e-10)? (gamma-t_pre/(tau))*exp(-tau/Asig)+t_pre/(tau) : gamma*(1-tau/Asig)+t_pre/Asig;
			}

			t_now+=t_pre;

			if (reach_end==0){
				if (cat_i!=-1) {
					Rprivate[omp_get_thread_num()][cat_i]+= back_rate_n*cat_weights[n][counter_eqk]*(ta/Asig)/gamma;
				}
				gamma= (DCFS_i==-1 | DCFS_whichpt[n][counter_eqk]==0)? gamma : gamma*exp(-DCFS[DCFS_i].cmb[DCFS_whichpt[n][counter_eqk]]/Asig);
				//if (n==1) printf("[%.5e:%.5e][%.5e:%.5e]: gamma=%f\n",tt0,tt1, tt0step, tt1step,gamma);
				if (isinf(gamma)){
					print_screen("*Warning: gamma==Inf, must choose larger Asig!*\n");
					print_logfile("*Warning: gamma==Inf, must choose larger Asig!*\n");
					err=1;
					errtot+=1;
				}
			}
			j0=next_TS+1;
			counter_eqk+=1;
		  }
		if (NeT) NeTprivate[omp_get_thread_num()][step]+=NeX[m];
		if (ReT) ReTprivate[omp_get_thread_num()][step]+=back_rate_n*(ta/Asig)/gamma;

		if (all_NeT) all_NeT[step][m]=NeX[m];

		step+=1;
		tt0step=tt1step;
		tt1step+=dt_step;
		if (fabs(tt1step-tt1)<tol0) tt1step=tt1;	//to avoid floating point error.

	}

	if (last) gamma_init[n]=gamma;	//update gamma_init with final value;
	if (ReT) ReX[m]=back_rate_n*(ta/Asig)/gamma;	//instantaneous seismicity rate at the end of the time step;

  }

  //collect values of R, NeT, ReT from threads:
  //final rate and tot no. of events given by the sum over all grid points:
  for (int i=0; i<nts; i++){
    if (ReT) ReT[i]=0.0;
    if (NeT) NeT[i]=0.0;
  }
  for (int t=0; t<nthreads; t++){
	  if (R) for (int eq=1; eq<=cat.Z; eq++) R[eq]+=Rprivate[t][eq];
	  if (NeT) for (int i=0; i<nts; i++) {
		if (NeT) NeT[i]+= (i==0) ? NeTprivate[t][i] : NeTprivate[t][i]-NeTprivate[t][i-1];
		if (ReT) ReT[i]+= ReTprivate[t][i];
	  }
  }

  free_d2array(Rprivate,0,nthreadstot-1, 0, cat.Z);
  free_d2array(NeTprivate,0,nthreadstot-1, 0, nts-1);
  free_d2array(ReTprivate,0,nthreadstot-1, 0, nts-1);
  return(errtot);

}
