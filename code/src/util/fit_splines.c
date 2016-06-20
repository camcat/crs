
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

#include "fit_splines.h" 
#include "../defines.h"
#include "error.h"

#ifdef _no_numericalrecipes

  void fit_splines(double *t, double *t2, int TS, int TS2, int N, double **slip_before, double *slip_before_err, double ***slip_after,
                int early_inter_mode, long *seed){
	
	print_logfile("Error: spline mode can not be activated if Numerical Recipes are not installed. Exit.\n");
	print_screen("Error: spline mode can not be activated if Numerical Recipes are not installed. Exit.\n");	
	return;

   }

#else

   #include "../nr/localnr.h"
   #include <math.h>   
   #include "interp_quad.h"

    void fit_splines(double *t, double *t2, int TS, int TS2, int N, double **slip_before, double *slip_before_err, double ***slip_after,
    		int early_inter_mode, long *seed){
    // t spans: 1...TS1. Contains initial points.
    // t2 spans: 1...TS2+1. Contains final sampled points.
    // early_inter_mode: 1 for linear interpolation between 0 and first time step, 2 for quadratic.
    // slip_before_err: error associated with each patch. if not given, assumed 1m.
    
    int NIT=500;
    double *e, *et;
    double f=0.01, sumas;	//error is fraction of maximum value: e(t)=f*max(t);
    float yp1,ypn;
    double sign, sum=0;
    int normalize=0, mc=1, no_oscillations=1;	//mc (=monte carlo) controls if errors (e) are included. space_error controls if they are space dependent (only available for Parkfield).
    double grenz0=1.5, grenz, dum;
    float *tf, *t2f, *sp, *sp2, **s2f, **ds2dt2, sp_value, sp_value2; //same as t, but floats(needed by spline function).
    
    if ((*slip_after)==(double **) 0) {
	(*slip_after)=d2array(1,N,1,TS2);
	if (!(*slip_after)) memory_error_quit;
    }
    
    float **a, **b;
    int n=3;
    
    a = matrix(1, n, 1, n);
    b = matrix(1, n, 1, N);
    
    yp1=1e30;
    ypn=1e30;
    
    e = darray(1,N);
    et = darray(1,TS);
    tf=vector(0,TS);
    sp = vector(0,TS);
    sp2 = vector(1,TS2);
    ds2dt2 = matrix(1, N, 1, TS);
    t2f=vector(1,TS2);
    s2f=matrix(1,N,1,TS2);
    
    for (int tt=1; tt<=TS; tt++) tf[tt]=(float) t[tt];
    for (int tt=1; tt<=TS2; tt++) t2f[tt]=(float) t2[tt];
    
    //calculate largest value at each point in time (error is set to be a fraction of this: e(t)=f*max_slip(t)).
    
    if (mc==1){
    	for (int tt=1; tt<=TS; tt++){
    		et[tt]=0.0;
    		for (int p=1; p<=N; p++) et[tt]=fmax(et[tt],f*fabs(slip_before[p][tt]));
    	}
    	if (slip_before_err== (double *) 0) for (int h=1; h<=N; h++) e[h]=1.0;			//no space dependent error.
    	else e = slip_before_err;
    }
    else{
    	for(int h=1; h<=N; h++) e[h]=0.0;	//no errors.
    	NIT=2;
    }
    
    //-----find first derivative for each point at x=0 (slope of quadratic that fits first 3 points).
    
    for (int h=1; h<=N; h++){
    	for (int y1=1;y1<=n;y1++){
    		b[y1][h]=0;
    		for (int y2=n;y2>=1;y2--){
    			a[y1][n-y2+1]=(float) pow(t[y1],y2-1);
    			b[y1][h]=(float) slip_before[h][y1];
    		}
    	}
    	gaussj(a, n, b, 1);
    }
    
    for (int h=1; h<=N; h++) sum+=slip_before[h][TS];
    sign=sum/fabs(sum);
    
    tf[0]=sp[0]=0.0;		//no errors on the requirement y(0)=0; (this is not used in finding splines, for stability).
    
    
    for (int h=1; h<=N; h++){
    
    	sumas=0.0;
    	for (int ts=1; ts<=TS2; ts++) sp2[ts]=0;
    	for (int ts=1; ts<=TS; ts++)	sumas+=slip_before[h][ts];
    
    	if (fabs(sumas)<=0.025) e[h]=0;		// if the patch has almost 0 slip at all times, fitted monte carlo curve is overestimated -> force curve to be close to 0.
    	for (int ts=1;ts<=TS2; ts++) (*slip_after)[h][ts]=0.0;
    
    	for (int iter=1; iter<=NIT; iter++){
    		for (int t0=1; t0<=TS; t0++){
    			grenz=grenz0;
    			dum=100*grenz;
    			while(fabs(dum)>grenz) {
    				dum=2.0*gasdev(seed)-1.0;	//-1 to 1.
    				(*seed)=-1.0*fabs(*seed);
    			}
    			// fill them backwards so that later value is kept if sp[t]>sp[t-1].
    			sp[TS+1-t0]= mc ? (float)(slip_before[h][TS+1-t0]+dum*e[h]*et[TS+1-t0]) : slip_before[h][TS+1-t0];		// fill them backwards so that later value is kept if sp[t]>sp[t-1].
    			//rules to avoid oscillations (different for last point to avoid overestimating the curve).
    			if (no_oscillations){
    				if (t0>2 && sign*sp[TS+1-t0]>sign*sp[TS+2-t0] && mc==1){
    					if (iter % 2 == 0) sp[TS+1-t0]=sp[TS+2-t0];
    					else sp[TS+2-t0]=sp[TS+1-t0];
    				}
    				if (sign*sp[TS+1-t0]<0  && mc==1) sp[TS+1-t0]=0;
    			}
    		 }
    
    //		spline(tf-1,sp-1,TS+1,b[n-1][h],ypn,ds2dt2[h]);					// y(0)=0, dydt(0) fitted to polynomial (degree 2 or 3);
    //		spline(tf-1,sp-1,TS+1,yp1,ypn,ds2dt2[h]);						// y(0)=0, dydt(0) fitted to straight line (natural spline);
    		spline(tf,sp,TS,yp1,ypn,ds2dt2[h]);						//no constrain on y(0), dydt(0) fitted to line (natural spline). Da best.;
    
    		int early_fit=1;
    		double *sp_values;
    		while (early_fit<TS2 && t2[early_fit]<t[1]) early_fit++;
    		early_fit--;
    		if (early_fit==TS2) {
    			print_screen("* Warning: all time steps are earlier than first input step, spline information not used*\n");
    			print_logfile("* Warning: all time steps are earlier than first input step, spline information not used*\n");
    		}
    		sp_values=darray(1,early_fit);
    
    		switch (early_inter_mode){
    			case 1:
    				for (int ts=1;ts<=early_fit; ts++) (*slip_after)[h][ts]+=(sp[1]*(t2[ts]/t[1])*(1.0/NIT));
    				break;
    			case 2:
    				splint(tf,sp,ds2dt2[h],TS,t2f[early_fit+1],&sp_value);
    				splint(tf,sp,ds2dt2[h],TS,t2f[early_fit+1]+0.0001,&sp_value2); //to obtain derivative.
    				interp_quad(0.0, t2f[early_fit+1], sp_value, (sp_value2-sp_value)/0.0001, t2, sp_values, early_fit);
    				for (int ts=1;ts<=early_fit; ts++) (*slip_after)[h][ts]+=(double) (sp_values[ts]*(1.0/NIT));
    				break;
    			case 0:
    				break;
    		}
    
    
    		for (int ts=early_fit+1;ts<=TS2; ts++) {
    			splint(tf,sp,ds2dt2[h],TS,t2f[ts],&sp_value);
    			(*slip_after)[h][ts]+=(double) (sp_value*(1.0/NIT));
    
    		}
    
    		}
    }
    
    	if (normalize){
    		int t0=TS, t1=0;
    		double final_cumslip0=0.0, final_cumslip1=0.0, dslip;
    
    		while (t0>0 && t[t0]>t2[TS2]) t0--;	//find last input time step within output time span.
    		for (int n=1; n<=N; n++) final_cumslip0+= slip_before[n][t0];
    		while(t1<TS2 && t2[t1+1]<t[t0]) t1++; //find last input time step before selected output time.
    		for (int h=1; h<=N; h++) {
    			if (t2[t1+1]==t[t0]) {
    				final_cumslip1+= (*slip_after)[h][t1+1];
    			}
    			else {
    				dslip=((*slip_after)[h][t1+1]-(*slip_after)[h][t1+1])*(t[t0]-t2[t1])/((t2[t1+1]-t2[t1]));	//linear interpolation.
    				final_cumslip1+= (*slip_after)[h][t1]+dslip;
    			}
    		}
    
    		for (int h=1; h<=N; h++) for (int ts=1; ts<=TS2; ts++) (*slip_after)[h][ts]*=final_cumslip0/final_cumslip1;	//normalize to keep final total slip constant.
    	}
    }
#endif
