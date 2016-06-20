
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


#include "find_timesteps.h"

int findtimestepsomori(double te, double t0, double t1, double K, double p, double c, double *t, double *K0, int *L) {
/*
 * Calculates time steps with an increasing spacing, a function of the form dt~(t+c)^p is used: for p=1 and a logarithmic stressing history,
 * this gives equal stresses within each time step. (NB: t refers to the start of each aseismic event).
 *
 * Input:
 * 	te: earthquake time (or in general, start time of the aseismic event).
 * 	t0,t1: start and end time of the period for which time steps are calculated.
 *
 * 	K, p, c: parameters controlling the shape of the function used to estimate time steps: *
 *
 * Output:
 *  t is populated with the time steps (NB must be allocated beforehand).
 *  K0: value such that t_{i}=t_{i-1}+K(t+c-teq)^p. Ignored if NULL.
 *  L: number of time steps.
 */

	double dfdt, dt, tnow;
	int N, j=0, err=0;

	dt=K*pow(t0+c-te,p);
	N=(int) (t1-t0)*(1.0/dt);			//N is always > than the needed number of elements since the first derivative is monotonically decreasing.

	tnow=t0;
	if (t) t[0]=t0;
	while (tnow<=t1){
		dt=K*pow(tnow+c-te,p);
		tnow+=dt;
		if (t) t[j+1]=tnow;
		j+=1;
	}

	if (j>1000) print_screen("\n ** Warning: findtimesteps.m produced %d time steps! **\n",j);

	*L=j;
	if (K0) *K0= K;	//so that t_{i}=t_{i-1}+K(t+c-teq)^p

	return err;
}
