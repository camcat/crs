
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


#include "GR.h"

double *assign_GRnorm(double *mags, int N, double b, int Minf){
/* Calculates normalized weights of seismicity within magnitude bins mags, based on a Gutenberg-Richter distribution.
 *
 * Input:
 *  mags: vector containing the center of sorted, equally spaced magnitude bins. Indices: [1...N].
 *  b: b-value of GR distrib: N(M>=m)=10^(a-b*m).
 *  Minf: flag indicating if last bin is open.
 */

	double dm;
	double norm;
	double *weights;
	double one;
	double mmax, mmin;
	double a=1;	//doesn't matter since weights are normalized.

	if (N==1) {
		weights=darray(1,1);
		weights[1]=1;
	}

	else {
		weights=darray(1,N);
		dm=mags[2]-mags[1];
		mmin=mags[1]-0.5*dm;
		mmax=mags[N]+0.5*dm;
		norm=Minf? pow(10,a-b*mmin) : pow(10,a-b*mmin)-pow(10,a-b*mmax);
		for (int i=1; i<=N; i++){
			mmin=mags[i]-0.5*dm;
			mmax=mags[i]+0.5*dm;
			weights[i]=(pow(10,a-b*mmin)-pow(10,a-b*mmax))/norm;
		}
		if (Minf) weights[N]=(pow(10,a-b*mmin))/norm;
	}

	return weights;

}

double Mc_maxcurv(double *mags, int N){
/* Returns magnitude of completeness using maximum curvature method;
 * NB: events are binned so that each bin has same no of counts, and normalized by bin size.
 *
 * Input:
 *  mags: list of magnitudes. Range [0...N-1].
 */

	double *bin_c, *bin_h;
	double f=0.0;
	double *mags2;
	int Nbin=floor(sqrt(N)), ind;

	mags2=darray(0,N-1);
	for (int i=0; i<N; i++) mags2[i]=mags[i];	//since next function sorts (and overwrites) array.

	// Bin events so that each magnitude bin has the same no. of events:
	// To get a pdf, bin_h is the bin count normalized by bin width.
	bin_equnumber(mags2,N, Nbin, &bin_c, &bin_h);

	//find value with largest count:
	for (int i=0; i<Nbin; i++){
		if (f<=bin_h[i]){
			f=bin_h[i];
			ind=i;
		}
	}

	free_darray(mags2, 0, N);
	return bin_c[ind];
}

double calculatebvalue(double *mags, int N, double Mc){
/* Calculates the Gutenberg-Richter b-value for a distribution of magnitudes, selecting only events above completeness Mc.
 *
 * Input:
 *  mag: list of magnitudes (Range: [0...N-1])
 *  Mc: completeness magnitude.
 *
 * Returns:
 *  b-value.
 */

  double  *mags2;
  double dM, dum, mean, b, db;
  int Z;

  mags2=darray(0,N-1);

  Z=0;
  mean=0.0;
  //pick magnitudes above completeness, calculate mean magnitude and fill in new vector.
  for (int i=0; i<N; i++){
	  if (mags[i]>=Mc){
		  mags2[Z]=mags[i];
		  mean+=mags[i];
		  Z+=1;
	  }
  }
  if (Z) Z-=1;
  mean*=(1.0/Z);

  qsort (mags2, Z, sizeof(double), &compare);

  dM=1000.0;
  for(int i=1;i<Z;i++) {
	  dum=(mags2[i]-mags2[i-1]);
	  if(dum!=0.0 && dum<dM) dM=dum;
  }

// Maximum Likelihood fit result:
  b  = log10(exp(1.0))/(mean-(Mc-0.5*dM));
  db = 1.96*b/sqrt(1.0*Z);

  return b;
}

int bin_equnumber(double *v, int N, int Nbin, double **bin_c, double **norm_count){
/* Bins a sample into Nbin, such that each bin has the same no. of elements. Outputs normalized count and bin center.
 *
 * Input:
 *  v: sample. Range [0...N-1].
 *  Nbin: no. of bins.
 *
 * Output:
 *  bin_c: center of each bin.
 *  norm_count: the no. of events per bin divided by bin size.
 */

	double bot, top;
	int Nel, Nellast;

	qsort (v, N, sizeof(double), &compare);
	Nel=floor(N/Nbin);

	*bin_c=darray(0,Nbin-1);
	*norm_count=darray(0,Nbin-1);

	for (int i=0; i<Nbin-1; i++){

		bot= (i==0) ? v[0] : 0.5*(v[Nel*i]+v[Nel*i-1]);
		top= 0.5*(v[Nel*(i+1)-1]+v[Nel*(i+1)]);

		(*bin_c)[i]=0.5*(top+bot);
		(*norm_count)[i]=Nel/(top-bot);
	}

	bot=0.5*(v[Nel*(Nbin-1)]+v[Nel*(Nbin-1)-1]);
	top=v[N-1];
	Nellast=N-Nel*(Nbin-1);
	(*bin_c)[Nbin-1]=0.5*(top+bot);
	(*norm_count)[Nbin-1]=Nellast/(top-bot);

	return (0);

}

int compare (const void * a, const void * b){
	/* Compares two numbers and returns [-1 0 1] depending on which number is larger.
	 */
	if (*(double*)a==*(double*)b) return 0;
	else if (*(double*)a>*(double*)b) return 1;
	else return -1;
}
