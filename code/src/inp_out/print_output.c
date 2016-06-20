
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


#include "print_output.h"

int sum_DCFS(struct pscmp *DCFS, double **cmb, int N, int Ntot){
/* Adds up all DCFS[0...N-1], and returns vector cmb containing cumulative field.
 *
 * Input:
 *  Ntot: total no. of grid points.
 *
 * Output:
 *  cmb: array containing the sum of stresses from DCFS[x].cmb.
 *  if *cmb==NULL memory is allocated, otherwise not (make sure vector has correct no. of elements).
 */

	int i;

	if (*cmb==NULL) *cmb=darray(1,Ntot);
	for (int k=1; k<=Ntot; k++) (*cmb)[k]=0.0;


	for (int n=0; n<N; n++){
		i=1;
		for (int k=1; k<=Ntot; k++){
			if (i<=DCFS[n].nsel && DCFS[n].which_pts[i]==k) {
				(*cmb)[k]+=DCFS[n].cmb[i];
				i++;
			}
		}
	}

	return 0;
}

int sum_DCFSrand(double **DCFSrand, double **cmb, int TS, int N){

	/* Sums elements in DCFSrand[0...TS-1][0...N-1]
	 * Returns them in array cmb[1...N]
	 *
	 * Memory for cmb will be allocated if *cmb==NULL.
	 */

	if (*cmb==NULL) *cmb=darray(1,N);
	for (int k=1; k<=N; k++) (*cmb)[k]=0.0;

	for (int n=0; n<N; n++){
		for (int k=0; k<TS; k++){
			(*cmb)[n]+=DCFSrand[k][n];
		}
	}

	return 0;
}

int print_grid(char *fname, struct pscmp DCFS, struct crust crst, double *rate){
/* if rate is a null pointer, prints out the coulomb stress field (DCFS.cmb).
 * uses refined grid geometry.
 */

	double *r;
	double Lon1, Lon2, Lat1, Lat2, D1, D2;
	FILE *fout;
	int i, Ntot=crst.uniform? (crst.nLat*crst.nLon*crst.nD) : crst.N_allP;

	r=(rate==NULL) ? DCFS.cmb : rate;

	fout=fopen(fname,"w");
	if (fout==NULL){
		print_screen("Error: file %s could not be opened (print_cmb). \n",fname);
		print_logfile("Error: file %s could not be opened (print_cmb). \n",fname);
		return(1);
	}
	i=1;
	for (int k=1; k<=Ntot; k++){
		Lon1=crst.lon[k]-crst.dlon/2.0;
		Lon2=crst.lon[k]+crst.dlon/2.0;
		Lat1=crst.lat[k]-crst.dlat/2.0;
		Lat2=crst.lat[k]+crst.dlat/2.0;
		D1=crst.depth[k]-crst.ddepth/2.0;
		D2=crst.depth[k]+crst.ddepth/2.0;
		if (i<=DCFS.nsel && DCFS.which_pts[i]==k) {
			//rate would contain all elements, but DCFS.cmb only those selected.
			if (rate) fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, 3.0, 8.0, r[k]);
			else fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, 3.0, 8.0, r[i]);
			i++;
		}
		else fprintf(fout,"%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \n", Lon1, Lon2, Lat1, Lat2, D1,D2, 3.0, 8.0, 0.0);
	}
	fclose(fout);

	return(0);
}

int print_slipmodel(char* filename, struct eqkfm *eqfm1, int NF){
// filename with extension; NF=no. of faults (elements of eqkfm).

	FILE *fout;
	int err=0;

	fout=fopen(filename,"w");
	if (fout==NULL) {
		print_screen("Warning: file %s could not be written (by function: print_slipmode).\n",filename);
		return 1;
	}
	else{
		for (int f=0; f<NF; f++){
			switch (eqfm1[f].whichfm){
				case 1:
					fprintf(fout, "%d   %.4lf   %.4lf   %.3lf   %.2lf   %.2lf   %.3lf   %.3lf   %d   %d   %.5lf\n", f+1, eqfm1[f].lat, eqfm1[f].lon, eqfm1[f].depth, eqfm1[f].L, eqfm1[f].W, eqfm1[f].str1, eqfm1[f].dip1, eqfm1[f].np_st, eqfm1[f].np_di, eqfm1[f].t);
					break;
				case 2:
					fprintf(fout, "%d   %.4lf   %.4lf   %.3lf   %.2lf   %.2lf   %.3lf   %.3lf   %d   %d   %.5lf\n", f+1, eqfm1[f].lat, eqfm1[f].lon, eqfm1[f].depth, eqfm1[f].L, eqfm1[f].W, eqfm1[f].str2, eqfm1[f].dip2, eqfm1[f].np_st, eqfm1[f].np_di, eqfm1[f].t);
					break;
				case 0:
					print_screen("Warning: whichfm=0, using first plane (print_slipmodel).\n");
					fprintf(fout, "%d   %.4lf   %.4lf   %.3lf   %.2lf   %.2lf   %.3lf   %.3lf   %d   %d   %.5lf\n", f+1, eqfm1[f].lat, eqfm1[f].lon, eqfm1[f].depth, eqfm1[f].L, eqfm1[f].W, eqfm1[f].str1, eqfm1[f].dip1, eqfm1[f].np_st, eqfm1[f].np_di, eqfm1[f].t);
					break;
				default:
					print_screen("Warning: ambiguos focal plane  -> output not written (print_slipmodel).\n");
					err+=1;
					continue;		//associated with for loop (not switch).
			}
			for (int p=1; p<=eqfm1[f].np_st*eqfm1[f].np_di; p++) {
				fprintf(fout, "%12.5lf\t%12.5lf\t%12.5lf\t%12.5lf\t%12.5lf\n", eqfm1[f].pos_s[p], eqfm1[f].pos_d[p], eqfm1[f].slip_str[p], eqfm1[f].slip_dip[p], 0.0);
			}
		}
		fclose(fout);
	}
	return err;
}
