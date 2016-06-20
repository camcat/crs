
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


#include "WellsCoppersmith.h"

#ifdef _CRS_MPI
	#include "mpi.h"
#endif

int WellsCoppersmith(double M, double rake, double *L, double *W, double *S) {

/* Estimate earthquake size and slip based on:
 *
 * Wells, D. L., & Coppermith, K. J. (1994). New Empirical Relationships among Magnitude, Rupture Length, Rupture Width, Rupture Area, and Surface Displacement.
 * Bull. Seism. Soc. Am., 84(4), 974–1002.
 *
 * It uses the value of surface area, and assumes a square model.
 * L,W in km; S in m.
 *
 */

	double a_RLD, sd_a_RLD,	b_RLD, sd_b_RLD, a_RW, sd_a_RW, b_RW, \
		   sd_b_RW, a_AD, sd_a_AD, b_AD, sd_b_AD, a_RA, sd_a_RA, \
		   b_RA, sd_b_RA;

	if(rake >=45 && rake <135) {	//reverse
		a_RLD=-2.42; sd_a_RLD=0.21;
		b_RLD=0.58; sd_b_RLD=0.03;

		a_RW=-1.61;	sd_a_RW=0.2;
		b_RW=0.41;	sd_b_RW=0.03;

		a_AD=-0.74;	sd_a_AD=1.4;	// Not significant at 95% level (see paper).
		b_AD=0.09;	sd_b_AD=0.21;	// Not significant at 95% level (see paper).

		a_RA=-3.99;	sd_a_RA=0.36;
		b_RA=0.98;	sd_b_RA=0.06;
	}
	else {
		if(rake>=225 && rake <315){	//normal
			a_RLD=-1.88;sd_a_RLD=0.37;
			b_RLD=0.50;	sd_b_RLD=0.06;

			a_RW=-1.14;	sd_a_RW=0.28;
			b_RW=0.35;	sd_b_RW=0.05;

			a_AD=-4.45;	sd_a_AD=1.59;
			b_AD=0.63;	sd_b_AD=0.24;

			a_RA=-2.87;	sd_a_RA=0.5;
			b_RA=0.82;	sd_b_RA=0.08;
		}
		else {
			if (rake>=0 && rake <360){
				a_RLD=-2.57; sd_a_RLD=0.12;
				b_RLD=0.62;	sd_b_RLD=0.02;

				a_RW=-0.76; sd_a_RW=0.12;
				b_RW=0.27;	sd_b_RW=0.02;

				a_AD=-6.32;	sd_a_AD=0.61;
				b_AD=0.9;	sd_b_AD=0.09;

				a_RA=-3.42;	sd_a_RA=0.18;
				b_RA=0.9;	sd_b_RA=0.03;
			}
			else {
				print_screen("Illegal value of rake in function WellsComppersmith!\n");
				print_logfile("Illegal value of rake in function WellsComppersmith!\n");
				return(1);
			}
		}
	}

//	use values for width and depth of rupture:
//	(*L)=pow(10.0,a_RLD + b_RLD*M);
//	(*W)=pow(10.0,a_RW + b_RW*M);

// use values for rupture area and assume square patch.
	(*L)=pow(10.0,0.5*(a_RA + b_RA*M));
	(*W)=(*L);
	(*S)=pow(10.0,a_AD + b_AD*M);

	return(0);
}
