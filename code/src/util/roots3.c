
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


/* Translated from fortran functions written by Ronjiang Wang, 2008,
 * part of the software PSCMP.
 */


#include "roots3.h"

#include <math.h>

#include "../defines.h"

void roots3(double b, double c, double d, double *x){

//finding 3 real roots of eq: x**3 + b*x**2 + c*x + d = 0
// x is vector with indices [1,2,3].

 double p,q,n,u,delta;

 if(d==0.0){
	x[1]=0.0;
	p=b*b-4.0*c;
	if(p<0.0 & extra_verbose) {
		print_screen("* Error in roots3: not all roots are real!*");
		print_logfile("* Error in roots3: not all roots are real!*");
	}
	else{
	  x[2]=0.5*(-b+sqrt(p));
	  x[3]=0.5*(-b-sqrt(p));
	}
 }
  else{
	p=c-b*b/3.0;
	q=d-b*(c-2.0*b*b/9.0)/3.0;
	if(4.0*pow(p,3)+27.0*pow(q,2)>0.0  & extra_verbose) {
		print_screen("* Error in roots3: not all roots are real!*");
		print_logfile("* Error in roots3: not all roots are real!*");
	}
	else{
	  n=sqrt(-4.0*p/3.0);
	  u=acos(-0.5*q/sqrt(-pow(p,3)/27.0))/3.0;
	  delta=8.0*atan(1.0)/3.0;
	  x[1]=n*cos(u)-b/3.0;
	  x[2]=n*cos(u+delta)-b/3.0;
	  x[3]=n*cos(u+2.0*delta)-b/3.0;
	}
  }
  return;
 }
