
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

#include "prestress.h"
#include <math.h>
#include <stddef.h>
#include "../defines.h"
#include "../util/moreutil.h"

#include "../util/util1.h"

void prestress(double s1, double s2, double s3, double strike, double dip, double rake, double p,double f, double ***s){
/*
 * Determines regional stress tensor using the known principal stresses and master fault mechanism, assuming that the master
 * fault follows the Coulomb failure criterion.
 *
 *
 * Input:
 *
 * s1, s2, s3: principal stresses.
 * strike, dip, rake: mechanism of optimally oriented planes.
 * p, f: pore pressure, friction coefficient
 *
 * Output:
 * 	s : stress tensor Range [1...3][1...3].
 */

//	local memories:
  
      double  alpha,st,di,ra,cmb1,cmb2,cmb3, cmb;
      double  ns[3],ts[3],rst[3],rdi[3];
      double  sig[3],rot[3][3];

      *s=d2array(1,3,1,3);
      if(s1==0.0 && s2==0.0 && s3==0.0) return;

      cmb1=0.5*fabs(s2-s3)*sqrt(1+f*f)+f*(0.5*(s2+s3)+p);
      cmb2=0.5*fabs(s3-s1)*sqrt(1+f*f)+f*(0.5*(s3+s1)+p);
      cmb3=0.5*fabs(s1-s2)*sqrt(1+f*f)+f*(0.5*(s1+s2)+p);
      cmb=fmax(fmax(cmb1,cmb2),cmb3);

      if(s1==s2 && s2==s3)return;

      if(cmb==cmb1){
        sig[2]=s1;
        sig[0]=fmax(s2,s3);
        sig[1]=fmin(s2,s3);
      }
      else {
		if(cmb==cmb2){
		  sig[0]=fmax(s3,s1);
		  sig[1]=fmin(s3,s1);
		  sig[2]=s2;
		}
		else{
		  sig[0]=fmax(s1,s2);
		  sig[1]=fmin(s1,s2);
		  sig[2]=s3;
		}
      }
      
////determine principal stress orientation
      alpha=0.5*atan2(1.0,f);
      st=strike*DEG2RAD;
      di=dip*DEG2RAD;
      ra=rake*DEG2RAD;
//	 normal vector:
      ns[0]=sin(di)*cos(st+0.5*PI);
      ns[1]=sin(di)*sin(st+0.5*PI);
      ns[2]=-cos(di);
//	 strike vector:
      rst[0]=cos(st);
      rst[1]=sin(st);
      rst[2]=0.0;
//	 dip vector:
      rdi[0]=cos(di)*cos(st+0.5*PI);
      rdi[1]=cos(di)*sin(st+0.5*PI);
      rdi[2]=sin(di);
//	 slip vector:
      for (int i=0; i<3; i++) ts[i]=rst[i]*cos(ra)-rdi[i]*sin(ra);
	      
      for (int i=0; i<3; i++){
        rot[i][0]=ns[i]*cos(alpha)+ts[i]*sin(alpha);
        rot[i][1]=ns[i]*sin(alpha)-ts[i]*cos(alpha);
      }
      
      rot[0][2]=ts[1]*ns[2]-ts[2]*ns[1];
      rot[1][2]=ts[2]*ns[0]-ts[0]*ns[2];
      rot[2][2]=ts[0]*ns[1]-ts[1]*ns[0];
      
      for (int j=0; j<3; j++){
        for (int i=0; i<=j; i++){
          (*s)[i+1][j+1]=0.0;
          for (int k=0; k<3; k++){
        	  (*s)[i+1][j+1]=(*s)[i+1][j+1]+sig[k]*rot[i][k]*rot[j][k];
		  }
		}
      }
      for (int j=0; j<3; j++) for (int i=2; i>j; i--)(*s)[i+1][j+1]=(*s)[j+1][i+1];

     return;
}

double **prestress_eigen(double *s, double *str, double *dip){
	/* Calculates stress tensor knowing magnitude and orientation of principal stresses.
	 *
	 * Input:
	 *  s[1...3]: principal axis
	 *  str, dip: strike and dips of principal stresses axis. Range [1...3]
	 */

	double **Q, **QT, **L, **S, **ST;

	Q=d2array(1,3,1,3);
	QT=d2array(1,3,1,3);
	L=d2array(1,3,1,3);

	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++){
			L[i][j]=0;
		}
	}

	for (int i=1; i<=3; i++){
		L[i][i]=s[i-1];
		Q[1][i]=QT[i][1]=cos(DEG2RAD*str[i-1])*cos(DEG2RAD*dip[i-1]);
		Q[2][i]=QT[i][2]=sin(DEG2RAD*str[i-1])*cos(DEG2RAD*dip[i-1]);
		Q[3][i]=QT[i][3]=sin(DEG2RAD*dip[i-1]);
	}

	S=mtimesm3(ST=mtimesm3(Q, L, NULL), QT, NULL);
	free_d2array(ST, 1,3,1,3);
	return S;

}
