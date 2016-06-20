
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


#include "decluster.h"

int *decluster_catalog(struct catalog cat, double Mmain, double **weights, int d3){
/* Declusters catalog using  Knopoff-Gardner 74 criterium.
 *
 * Input:
 * 	cat:	catalog
 * 	Mmain:	magnitude above which events are selected as mainshocks;
 * 	d3: flag indicating if 3d distance (instead of horizontal distance) should be used.
 *
 * Output:
 *   Returns array containing flag (1/0) for selected/excluded events.
 *
 *   if *weights==NULL, memory will be allocated.
 *
 * NB: events must be sorted chronologically.
 *
 * Since Knopoff-Gardner declustering removes all events in a certain space-time window (including background events), the remaining earthquakes
 * should be rescaled to account for the time windows being removed by declustering. This is explained in more detail below.
 *
 */

	double D, T, d, *tnow, *time_missing, dt;
	int *sel=iarray(1,cat.Z);
	for (int i=1; i<=cat.Z; i++) sel[i]=1;

	if (weights){
		tnow=darray(1,cat.Z);
		time_missing=darray(1,cat.Z);
		if (!(*weights)) *weights=darray(1,cat.Z);
		for (int p=1; p<=cat.Z; p++) {
			tnow[p]=cat.tstart;
		}
	}

	for (int i=1; i<=cat.Z; i++) {
		sel[i]=1;
		if (weights) time_missing[i]=0.0;
	}

	for (int i=1; i<=cat.Z; i++){
		if (cat.mag[i]>=Mmain){
			KG74(cat.mag[i], &D, &T);	//calculate spatial and temporal window.
			for (int j=1; j<=cat.Z; j++){
				if (sel[j]==0) continue;	//event has already been removed (aftershock of a previous mainshock).
				d= pow(cat.x0[j]-cat.x0[i],2)+pow(cat.y0[j]-cat.y0[i],2);	//horizontal distance
				if (d3) d+=pow(cat.depths0[j]-cat.depths0[i],2);			//3D distance
				d=sqrt(d);
				if (d<=D){	//decluster catalog:
					if((cat.t[j]-cat.t[i])>0 && (cat.t[j]-cat.t[i])<=T && j!=i) sel[j]=0;
					else {
						//calculate period missing at the location of this earthquake:
						dt=fmin(cat.tend, cat.t[i]+T)-fmax(cat.t[i], tnow[j]);	//time missing due to declustering; fmax(..) to avoid double counting it Ts overlap.
						if (dt>0) time_missing[j]+= dt;
						tnow[j]=fmax(fmin(cat.tend, cat.t[i]+T),tnow[j]);	//since cat.t[i]+T may be smaller than a previous cat.t[i]+T (tnow).
					}
				}
			}
		}
	}

	if (weights) {
		// weights are given to earthquakes in order to account for the time period missing by declustering. This is necessary since declustering removes
		// all the events in a given time period, i.e. not only the true aftershocks but also the background events.
		// For example if 5% of the total period T was removed by declustering at a given location, the number of earthquakes
		// should be rescaled: Ntrue=N*(1/0.95), since the earthquakes left after declustering occur in a period of 0.95T.
		// This is achieved by assigning a weight of T/(T-Tmissing) to the earthquakes left after declustering.
		// However, if the time period removed by declustering is too large, T/(T-Tmissing) may go to infinity; as Tmissing ->T, the time period T-Tmissing
		// may be too small to extend the seismicity to the entire period. Therefore the rescaling factor is capped to 10.0
		// (i.e. 90% of the time period is removed by declustering).
		for (int i=1; i<=cat.Z; i++) (*weights)[i]= (sel[i]==0)? 0.0 : fmin((cat.tend-cat.tstart)/(cat.tend-cat.tstart-time_missing[i]),10.0);
		free_darray(time_missing, 1, cat.Z);
		free_darray(tnow, 1, cat.Z);
	}
	return sel;
}

void KG74(double M, double *D, double *T){
/* calculation of Knopoff-Gardner 74 criterium for aftershock declustering
 * function modified from one written by Olga Zakharova.
 *
 * Input:
 * M= magnitude.
 *
 * Output:
 * D= spatial window.
 * T= temporal window.
 */

  if (D) *D = pow(10,(0.1238*M+0.983));
  if (T) *T = (M >= 6.5) ? pow(10,(0.032*M+2.7389)) : pow(10,(0.5409*M-0.547));

  return;
}
