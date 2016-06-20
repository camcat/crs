
/*   Original author (fortran code):
 *
 *   Yoshimitsu Okada, Natl.Res.Inst. for Earth Sci. & Disas.Prev, Tsukuba, Japan
 *   Reference: Okada, Y., 1992, Internal deformation due to shear and tensile faults in a half-space,
 * 		Bull. Seism. Soc. Am., 82, 1018-1040.
 *   Source code available at: http://www.bosai.go.jp/study/application/dc3d/download/DC3Dfortran.txt
 *
 *   Translated in C by Christoph Bach (2010).
 */


void UA(double XI, double ET, double Q, double DISL1, double DISL2, double DISL3, double *U,
		double Y11, double X11, double ALP2, double ALP1, double TT, double R, double ALE,
		double XI2, double Y32, double Q2, double SD, double R3, double FY, double D,
		double EY, double CD, double FZ, double Y, double EZ, double ALX, double GY, double GZ,
		double HY, double HZ);

void UB(double XI, double ET, double Q, double DISL1, double DISL2, double DISL3, double *U,
		double R, double D, double Y, double CD, double XI2, double Q2,
		double CDCD, double SDCD, double SD, double ALE, double Y11, double X11,
		double ALP3, double TT, double Y32, double R3, double FY, double EY, double FZ, double EZ,
		double GY, double GZ, double HY, double HZ, double SDSD);

void UC(double XI, double ET, double Q, double Z, double DISL1, double DISL2, double DISL3, double *U,
		double D, double R2, double R, double XI2, double X11, double Y11, double ET2, double CD,
		double SD, double R3, double Y32, double R5, double Y, double ALP4, double ALP5, double X32,
		double Q2, double SDSD, double SDCD, double CDCD);

void DCCON0(double ALPHA, double DIP,
		    double *ALP1, double *ALP2, double *ALP3, double *ALP4, double *ALP5,
		    double *SD, double *CD, double *SDSD, double *CDCD, double *SDCD, double *S2D, double *C2D);

void DCCON2(double XI, double ET, double Q, double KXI, double KET, double SD, double CD,
		    double *XI2, double *ET2, double *Q2, double *R2, double *R,
		    double *R3, double *R5, double *Y, double *D, double *TT,
		    double *ALX, double *X11, double *X32, double *ALE, double *Y11, double *Y32,
		    double *EY, double *EZ, double *FY, double *FZ, double *GY, double *GZ, double *HY, double *HZ);

void  DC3D(double ALPHA, double X, double Y, double Z, double DEPTH,
		   double DIP, double AL1, double AL2, double AW1, double AW2,
		   double DISL1, double DISL2, double DISL3,
		   double *UX,double *UY, double *UZ, double *UXX, double *UYX,
		   double *UZX, double *UXY, double *UYY, double *UZY,
		   double *UXZ, double *UYZ, double *UZZ, int *IRET);
