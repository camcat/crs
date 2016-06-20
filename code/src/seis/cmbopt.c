
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

#include "cmbopt.h"

void DCFScmbopt(struct pscmp *DCFS, int ind, struct crust crst){
	/* Calculates optimally oriented planes in a stress field given by the sum of background stress field (crst.S) and stress steps contained in DCFS.
	 * Also calculates Coulomb stress change due to DCFS[ind] on these planes.
	 * ind: gives the index of last element of DCFS that should be included in stress field calculation; (e.g. if 0, only one element).
	 *
	 * Results stored in: DCFS[ind].st1, di1, ra1[2]; largest stress changes from DCFS[ind] stored in DCFS[ind].cmb
	 */

	int j, k, ev;
	double cmb1, cmb2;
	double sxx, syy, szz, sxy, syz, sxz;
	double MaxDCFS=DCFS_cap;

	DCFS[ind].st1=darray(1,DCFS[ind].nsel);
	DCFS[ind].di1=darray(1,DCFS[ind].nsel);
	DCFS[ind].ra1=darray(1,DCFS[ind].nsel);
	DCFS[ind].st2=darray(1,DCFS[ind].nsel);
	DCFS[ind].di2=darray(1,DCFS[ind].nsel);
	DCFS[ind].ra2=darray(1,DCFS[ind].nsel);

	#pragma omp parallel for private(j, sxx, sxy, syy, syz, szz, sxz, cmb1, cmb2, k, ev)
	for (int i=1; i<=DCFS[ind].nsel; i++){
		j=DCFS[ind].which_pts[i];

		sxx=crst.S[1][1]; 		syy=crst.S[2][2];		szz=crst.S[3][3];
		sxy=crst.S[1][2]; 		syz=crst.S[2][3];		sxz=crst.S[1][3];

		//loop over events un to ind:
		for (ev=0; ev<=ind; ev++){
			k=0;
			//find if grid point is affected by event ev, and if so add its stresses:
			while (k<DCFS[ev].nsel && DCFS[ev].which_pts[k]<j) k++;
			if (DCFS[ev].which_pts[k]==j){
				sxx+=DCFS[ev].S[k][1][1];
				syy+=DCFS[ev].S[k][2][2];
				szz+=DCFS[ev].S[k][3][3];
				sxy+=DCFS[ev].S[k][1][2];
				syz+=DCFS[ev].S[k][2][3];
				sxz+=DCFS[ev].S[k][1][3];
			}
		}

		//find OOPs and resolve total stress (DCFS+S0):
		//out of the two, the oop closest to the regional mechanism (crst.str0[0], crst.dip0[0], crst.rake0[0]) is chosen.
		cmbopt(sxx, syy, szz, sxy, syz, sxz, 0, crst.fric, crst.str0[0], crst.dip0[0], crst.rake0[0], DCFS[ind].cmb+i, DCFS[ind].st1+i, DCFS[ind].di1+i, DCFS[ind].ra1+i, DCFS[ind].st2+i, DCFS[ind].di2+i, DCFS[ind].ra2+i);
		//resolve background stress on both planes and calculate largest DCFS:
		cmb1=resolve_S(DCFS[ind].S[i], DCFS[ind].st1[i], DCFS[ind].di1[i], DCFS[ind].ra1[i], crst.fric, NULL, 0.0, 0);
		cmb2=resolve_S(DCFS[ind].S[i], DCFS[ind].st2[i], DCFS[ind].di2[i], DCFS[ind].ra2[i], crst.fric, NULL, 0.0, 0);
		DCFS[ind].cmb[i]=fmax(cmb1,cmb2);
		if (DCFS[ind].cmb[i]>MaxDCFS) DCFS[ind].cmb[i]=MaxDCFS;
		if (DCFS[ind].cmb[i]<-MaxDCFS) DCFS[ind].cmb[i]=-MaxDCFS;
	}

	free_darray(DCFS[ind].st1,1,DCFS[ind].nsel);
	free_darray(DCFS[ind].di1,1,DCFS[ind].nsel);
	free_darray(DCFS[ind].ra1,1,DCFS[ind].nsel);
	free_darray(DCFS[ind].st2,1,DCFS[ind].nsel);
	free_darray(DCFS[ind].di2,1,DCFS[ind].nsel);
	free_darray(DCFS[ind].ra2,1,DCFS[ind].nsel);

}

void cmbopt(double sxx, double syy, double szz, double sxy, double syz, double szx, double p, double f, double st0, double di0, double ra0,
		double *cmb,double *st1, double *di1, double *ra1, double *st2, double *di2, double *ra2){
/*  Calculates coulomb stress with the optimal orientation. Taken from Ronjiang fortran function in pscmp program.
 *
 * 	Input:
 * 	 sxx, sxy, ...: stress tensor components
 * 	 p=pore pressure
 * 	 f=friction coefficient
 * 	 str0, di0, rake0: a reference focal mechanism: output focal mechanisms str1, di1, rake1 refer to mechanism closer to this.
 *
 * 	Output:
 * 	 cmb: max. Coulomb stress at the two optimally oriented fault planes
 * 	 st1, di1, ra1, [2]: strike, dip, rake of the two OOPs.
 */

      int j0,j1,j2,jmin,jmax;
      double b,c,d,s1,s2,s3,snn,alpha,am,swap, sig;
      double cmb1,cmb2,cmb3,det1,det2,det3,detmax,rmax;
      double *s,**r,**ns,**ts;

      s=darray(1,3);
      r=d2array(1,3,1,3);
      ns=d2array(1,3,1,2);
      ts=d2array(1,3,1,3);

      if(sxy==0.0 && syz==0.0 && szx==0.0){
        s[1]=sxx;
        s[2]=syy;
        s[3]=szz;
      }
      else{
        b=-(sxx+syy+szz);
        c=sxx*syy+syy*szz+szz*sxx-pow(sxy,2)-pow(syz,2)-pow(szx,2);
        d=sxx*pow(syz,2)+syy*pow(szx,2)+szz*pow(sxy,2)-2.0*sxy*syz*szx-sxx*syy*szz;
        roots3(b,c,d,s);
      }
      cmb1=0.5*fabs(s[2]-s[3])*sqrt(1+f*f)+f*(0.5*(s[2]+s[3])+p);
      cmb2=0.5*fabs(s[3]-s[1])*sqrt(1+f*f)+f*(0.5*(s[3]+s[1])+p);
      cmb3=0.5*fabs(s[1]-s[2])*sqrt(1+f*f)+f*(0.5*(s[1]+s[2])+p);
      *cmb=fmax(cmb1,fmax(cmb2,cmb3));
      *st1=0.0;
      *di1=0.0;
      *ra1=0.0;
      *st2=0.0;
      *di2=0.0;
      *ra2=0.0;
      if(*cmb==cmb1){
        s3=s[1];
        s1=fmax(s[2],s[3]);
        s2=fmin(s[2],s[3]);
      }
      else if(*cmb==cmb2){
		s1=fmax(s[3],s[1]);
		s2=fmin(s[3],s[1]);
		s3=s[2];
      }
   	  else{
		s1=fmax(s[1],s[2]);
		s2=fmin(s[1],s[2]);
		s3=s[3];
      }
      sig=0.5*((s1-s2)*f/sqrt(1+f*f)+s1+s2);
      s[1]=s1;
      s[2]=s2;
      s[3]=s3;
//determine eigenvectors (the principal stress directions)
      j0=0;
      if(s[1]==s[2]){
        j0=3;
        j1=1;
        j2=2;
      }
      else if(s[2]==s[3]){
		j0=1;
		j1=2;
		j2=3;
	  }
      else if(s[3]==s[1]){
		j0=2;
		j1=1;
		j2=3;
      }

      if(j0==0){
        jmin=1;
        jmax=3;
      }
      else{
        jmin=j0;
        jmax=j0;
        if (extra_verbose){
			print_screen("* Warning: more than two optimal rupture orientations! *");
			print_logfile("* Warning: more than two optimal rupture orientations! *");
        }
      }
      for (int j=jmin; j<=jmax; j++){
        det1=syz*syz-(syy-s[j])*(szz-s[j]);
        det2=szx*szx-(sxx-s[j])*(szz-s[j]);
        det3=sxy*sxy-(sxx-s[j])*(syy-s[j]);
        detmax=fmax(fabs(det1),fmax(fabs(det2),fabs(det3)));
        if(fabs(det1)==detmax){
          r[1][j]=det1;
          r[2][j]=(szz-s[j])*sxy-syz*szx;
          r[3][j]=(syy-s[j])*szx-syz*sxy;
        }
        else if(fabs(det2)==detmax){
          r[1][j]=(szz-s[j])*sxy-szx*syz;
          r[2][j]=det2;
          r[3][j]=(sxx-s[j])*syz-szx*sxy;
        }
        else{
          r[1][j]=(syy-s[j])*szx-sxy*syz;
          r[2][j]=(sxx-s[j])*syz-sxy*szx;
          r[3][j]=det3;
        }
      }
/*if any two eigenvalues are identical, their corresponding
eigenvectors should be redetermined by orthogonalizing
them to the 3. eigenvector as well as to each other*/

      if(j0 > 0){
        rmax=fmax(fabs(r[1][j0]),fmax(fabs(r[2][j0]),fabs(r[3][j0])));
        if(fabs(r[1][j0])==rmax){
          r[1][j1]=-r[2][j0];
          r[2][j1]=r[1][j0];
          r[3][j1]=0.0;
          r[1][j2]=-r[3][j0];
          r[2][j2]=0.0;
          r[3][j2]=r[1][j0];
          am=r[1][j1]*r[1][j2]/(pow(r[1][j1],2)+pow(r[2][j1],2));
          for (int i=1; i<=3; i++) r[i][j2]=r[i][j2]-am*r[i][j1];
        }
        else if(fabs(r[2][j0])==rmax){
          r[1][j1]=r[2][j0];
          r[2][j1]=-r[1][j0];
          r[3][j1]=0.0;
          r[1][j2]=0.0;
          r[2][j2]=-r[3][j0];
          r[3][j2]=r[2][j0];
          am=r[2][j1]*r[2][j2]/(pow(r[1][j1],2)+pow(r[2][j1],2));
          for (int i=1; i<=3; i++) r[i][j2]=r[i][j2]-am*r[i][j1];
        }
        else if(fabs(r[3][j0])==rmax){
          r[1][j1]=r[3][j0];
          r[2][j1]=0.0;
          r[3][j1]=-r[1][j0];
          r[1][j2]=0.0;
          r[2][j2]=r[3][j0];
          r[3][j2]=-r[2][j0];
          am=r[3][j1]*r[3][j2]/(pow(r[1][j1],2)+pow(r[3][j1],2));
          for (int i=1; i<=3; i++) r[i][j2]=r[i][j2]-am*r[i][j1];
        }
      }
      for (int j=1; j<=3; j++){
        am=sqrt(pow(r[1][j],2)+pow(r[2][j],2)+pow(r[3][j],2));
        for (int i=1; i<=3; i++) r[i][j]=r[i][j]/am;
      }

      alpha=0.5*atan(1.0/f);
      snn=s[1]*pow(cos(alpha),2)+s[2]*pow(sin(alpha),2);
//determine the two optimal fault-plane normals
      for (int i=1; i<=3; i++){
        ns[i][1]=r[i][1]*cos(alpha)+r[i][2]*sin(alpha);
        ns[i][2]=r[i][1]*cos(alpha)-r[i][2]*sin(alpha);
      }
//determine the direction of max. shear stress
      for (int j=1; j<=2; j++){
        am=sqrt(pow(ns[1][j],2)+pow(ns[2][j],2)+pow(ns[3][j],2));
        if (ns[3][j] > 0.0) am=-am;
        for (int i=1; i<=3; i++) ns[i][j]=ns[i][j]/am;

        ts[1][j]=(sxx-snn)*ns[1][j]+sxy*ns[2][j]+szx*ns[3][j];
        ts[2][j]=sxy*ns[1][j]+(syy-snn)*ns[2][j]+syz*ns[3][j];
        ts[3][j]=szx*ns[1][j]+syz*ns[2][j]+(szz-snn)*ns[3][j];
        am=sqrt(pow(ts[1][j],2)+pow(ts[2][j],2)+pow(ts[3][j],2));
        for (int i=1; i<=3; i++) ts[i][j]=ts[i][j]/am;
      }

//determine the two optimal focal mechanisms

      *st1=fmod(atan2(ns[2][1],ns[1][1])*180.0/PI+270.0,360.0);
	  *di1=acos(-ns[3][1])*180.0/PI;
      s1=cos((*st1)*PI/180.0);
      s2=sin((*st1)*PI/180.0);
      *ra1=acos(fmin(fmax(s1*ts[1][1]+s2*ts[2][1],-1.0),1.0))*180.0/PI;

      if (ts[3][1] > 0.0) *ra1=-(*ra1);
      *st2=fmod(atan2(ns[2][2],ns[1][2])*180.0/PI+270.0,360.0);
	  *di2=acos(-ns[3][2])*180.0/PI;
      s1=cos((*st2)*PI/180.0);
      s2=sin((*st2)*PI/180.0);
      *ra2=acos(fmin(fmax(s1*ts[1][2]+s2*ts[2][2],-1.0),1.0))*180.0/PI;
      if(ts[3][2] > 0.0) *ra2=-(*ra2);
      double msc1=mscorr(st0,di0,ra0,*st1,*di1,*ra1);
      double msc2=mscorr(st0,di0,ra0,*st2,*di2,*ra2);
      if(msc1<msc2){
        swap=*st1;
        *st1=*st2;
        *st2=swap;
        swap=*di1;
        *di1=*di2;
        *di2=swap;
        swap=*ra1;
        *ra1=*ra2;
        *ra2=swap;
      }

      free_darray(s,1,3);
      free_d2array(r,1,3,1,3);
      free_d2array(ns,1,3,1,2);
      free_d2array(ts,1,3,1,2);
      return;
}
