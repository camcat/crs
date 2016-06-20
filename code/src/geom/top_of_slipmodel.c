
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


#include <math.h>

#include "../defines.h"

void top_of_slipmodel(struct eqkfm* eqkfm, int N){
	/* Finds the top of a slip model (eqkfm) and writes into eqkfm.top
	 *
	 * Input:
	 * 	eqkfm: slip model structure [0...N-1]
	 * 	N: size of eqkfm;
	 *
	 * Output:
	 *  eqkfm[i].top is filled with depth of the top of the model.
	 */

	double z0, dz_dip;	//depth of reference point, distance along dip of closest patch;
	double len_dip;	//length of patch along dip;
	double dz;	//depth of top of shallowest patch w.r.t. reference point.
	int np;

	z0=1e30;
	for (int i=0; i<N; i++){
		np=eqkfm[i].np_di*eqkfm[i].np_st;
		len_dip=eqkfm[i].W/eqkfm[i].np_di;

		//find along dip distance of center of shallowest patch:
		dz_dip=1e30;
		for (int j=1; j<=np; j++) dz_dip=fmin(dz_dip, eqkfm[i].pos_d[j]);
		if (eqkfm[i].whichfm!=2) dz=(dz_dip-0.5*len_dip)*sin(DEG2RAD*eqkfm[i].dip1);
		else dz=(dz_dip-0.5*len_dip)*sin(DEG2RAD*eqkfm[i].dip2);

		z0=fmin(z0, eqkfm[i].depth+dz);
	}

	for (int i=0; i<N; i++) eqkfm[i].top=z0;

	return;
}
