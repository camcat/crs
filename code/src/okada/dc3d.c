
/*   Original author (fortran code):
 *
 *   Yoshimitsu Okada, Natl.Res.Inst. for Earth Sci. & Disas.Prev, Tsukuba, Japan
 *   Reference: Okada, Y., 1992, Internal deformation due to shear and tensile faults in a half-space,
 * 		Bull. Seism. Soc. Am., 82, 1018-1040.
 *   Swq
 *   ource code available at: http://www.bosai.go.jp/study/application/dc3d/download/DC3Dfortran.txt
 *
 *   Translated in C by Christoph Bach (2010).
 */


/* C********************************************************************   */
/* C*****                                                          *****   */
/* C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   */
/* C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   */
/* C*****              CODED BY  Y.OKADA ... SEP.1991              *****   */
/* C*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****   */
/* C*****                          JUL.1993                        *****   */
/* C********************************************************************   */
/* C                                                                       */
/* C***** INPUT                                                            */
/* C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           */
/* C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          */
/* C*****   DEPTH : DEPTH OF REFERENCE POINT                               */
/* C*****   DIP   : DIP-ANGLE (DEGREE)                                     */
/* C*****   AL1,AL2   : FAULT LENGTH RANGE (-STRIKE,+STRIKE)               */
/* C*****   AW1,AW2   : FAULT WIDTH RANGE  ( DOWNDIP, UPDIP)               */
/* C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              */
/* C                                                                       */
/* C***** OUTPUT                                                           */
/* C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               */
/* C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /              */
/* C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) ) */
/* C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                      */
/* C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )   */

#include <math.h>

#include "../defines.h"
#include "../util/error.h"

void UA(double XI, double ET, double Q, double DISL1, double DISL2, double DISL3, double *U, double Y11, double X11, double ALP2, double ALP1, double TT, double R, double ALE, double XI2, double Y32, double Q2, double SD, double R3,
		double FY, double D, double EY, double CD, double FZ, double Y, double EZ, double ALX, double GY, double GZ, double HY, double HZ) {

	/* C********************************************************************  */
	/* C*****    DISPLACEMENT AND STRAIN AT DEPTH [PART-A]             *****  */
	/* C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****  */
	/* C********************************************************************  */
	/* C                                                                      */
	/* C***** INPUT                                                           */
	/* C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                 */
	/* C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS             */
	/* C***** OUTPUT                                                          */
	/* C*****   U[12] : DISPLACEMENT AND THEIR DERIVATIVES                    */

	double DU[13];
	double XY, QX, QY;
	double F2, PI2, dum1, dum2, dum3, dum4, ALP2Q;
	int i;

	F2 = 2.0;
	PI2 = 6.283185307179586;

	for (i = 1; i <= 12; i++)
		U[i] = 0.0;

	XY = XI * Y11;
	QX = Q * X11;
	QY = Q * Y11;
	dum2 = ALP2 * XI;
	ALP2Q = ALP2 * Q;
	/* C====================================== */
	/* C=====  STRIKE-SLIP CONTRIBUTION  ===== */
	/* C====================================== */
	if (DISL1 != 0.0) {
		dum1 = ALP1 * XY;
		DU[1] = TT / F2 + dum2 * QY;
		DU[2] = ALP2Q / R;
		DU[3] = ALP1 * ALE - ALP2Q * QY;
		DU[4] = -ALP1 * QY - ALP2 * XI2 * Q * Y32;
		DU[5] = -dum2 * Q / R3;
		DU[6] = dum1 + dum2 * Q2 * Y32;
		DU[7] = dum1 * SD + dum2 * FY + D / F2 * X11;
		DU[8] = ALP2 * EY;
		DU[9] = ALP1 * (CD / R + QY * SD) - ALP2Q * FY;
		DU[10] = dum1 * CD + dum2 * FZ + Y / F2 * X11;
		DU[11] = ALP2 * EZ;
		DU[12] = -ALP1 * (SD / R - QY * CD) - ALP2Q * FZ;
		dum1 = DISL1 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum1 * DU[i];
	}
	/* C====================================== */
	/* C=====    DIP-SLIP CONTRIBUTION   ===== */
	/* C====================================== */
	if (DISL2 != 0.0) {
		dum3 = ALP1 * D * X11;
		dum4 = ALP1 * Y * X11;
		DU[1] = ALP2Q / R;
		DU[2] = TT / F2 + ALP2 * ET * QX;
		DU[3] = ALP1 * ALX - ALP2Q * QX;
		DU[4] = -dum2 * Q / R3;
		DU[5] = -QY / F2 - ALP2 * ET * Q / R3;
		DU[6] = ALP1 / R + ALP2 * Q2 / R3;
		DU[7] = ALP2 * EY;
		DU[8] = dum3 + XY / F2 * SD + ALP2 * ET * GY;
		DU[9] = dum4 - ALP2Q * GY;
		DU[10] = ALP2 * EZ;
		DU[11] = dum4 + XY / F2 * CD + ALP2 * ET * GZ;
		DU[12] = -dum3 - ALP2Q * GZ;
		dum1 = DISL2 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum1 * DU[i];
	}
	/* C======================================== */
	/* C=====  TENSILE-FAULT CONTRIBUTION  ===== */
	/* C======================================== */
	if (DISL3 != 0.0) {
		DU[1] = -ALP1 * ALE - ALP2Q * QY;
		DU[2] = -ALP1 * ALX - ALP2Q * QX;
		DU[3] = TT / F2 - ALP2 * (ET * QX + XI * QY);
		DU[4] = -ALP1 * XY + dum2 * Q2 * Y32;
		DU[5] = -ALP1 / R + ALP2 * Q2 / R3;
		DU[6] = -ALP1 * QY - ALP2Q * Q2 * Y32;
		DU[7] = -ALP1 * (CD / R + QY * SD) - ALP2Q * FY;
		DU[8] = -ALP1 * Y * X11 - ALP2Q * GY;
		DU[9] = ALP1 * (D * X11 + XY * SD) + ALP2Q * HY;
		DU[10] = ALP1 * (SD / R - QY * CD) - ALP2Q * FZ;
		DU[11] = ALP1 * D * X11 - ALP2Q * GZ;
		DU[12] = ALP1 * (Y * X11 + XY * CD) + ALP2Q * HZ;
		dum1 = DISL3 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum1 * DU[i];
	}
}

void UB(double XI, double ET, double Q, double DISL1, double DISL2, double DISL3, double *U, double R, double D, double Y, double CD, double XI2, double Q2, double CDCD, double SDCD, double SD, double ALE, double Y11, double X11,
		double ALP3, double TT, double Y32, double R3, double FY, double EY, double FZ, double EZ, double GY, double GZ, double HY, double HZ, double SDSD) {

	/* C******************************************************************** */
	/* C*****    DISPLACEMENT AND STRAIN AT DEPTH [PART-B]             ***** */
	/* C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
	/* C******************************************************************** */
	/* C                                                                     */
	/* C***** INPUT                                                          */
	/* C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                */
	/* C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS            */
	/* C***** OUTPUT                                                         */
	/* C*****   U[12] : DISPLACEMENT AND THEIR DERIVATIVES                   */

	double DU[13];
	double QX, QY, X, XY, RD, RD2, D11, AJ1, AJ2, AJ3, AJ4, AJ5, AJ6, AI1, AI2, AI3, AI4, AK1, AK2, AK3, AK4;
	double F1, F2, PI2, dum1, dum2, dum3;
	int i;

	F1 = 1.0;
	F2 = 2.0;
	PI2 = 6.283185307179586;

	RD = R + D;
	D11 = F1 / (R * RD);
	AJ2 = XI * Y / RD * D11;
	AJ5 = -(D + Y * Y / RD) * D11;
	if (CD != 0.0) {
		if (XI == 0.0) AI4 = 0.0;
		else {
			X = sqrt(XI2 + Q2);
			AI4 = F1 / CDCD * (XI / RD * SDCD + F2 * atan((ET * (X + Q * CD) + X * (R + X) * SD) / (XI * (R + X) * CD)));
		}
		AI3 = (Y * CD / RD - ALE + SD * log(RD)) / CDCD;
		AK1 = XI * (D11 - Y11 * SD) / CD;
		AK3 = (Q * Y11 - Y * D11) / CD;
		AJ3 = (AK1 - AJ2 * SD) / CD;
		AJ6 = (AK3 - AJ5 * SD) / CD;
	}
	else {
		RD2 = RD * RD;
		AI3 = (ET / RD + Y * Q / RD2 - ALE) / F2;
		AI4 = XI * Y / RD2 / F2;
		AK1 = XI * Q / RD * D11;
		AK3 = SD / RD * (XI2 * D11 - F1);
		AJ3 = -XI / RD2 * (Q2 * D11 - F1 / F2);
		AJ6 = -Y / RD2 * (XI2 * D11 - F1 / F2);
	}

	XY = XI * Y11;
	AI1 = -XI / RD * CD - AI4 * SD;
	AI2 = log(RD) + AI3 * SD;
	AK2 = F1 / R + AK3 * SD;
	AK4 = XY * CD - AK1 * SD;
	AJ1 = AJ5 * CD - AJ6 * SD;
	AJ4 = -XY - AJ2 * CD + AJ3 * SD;

	for (i = 1; i <= 12; i++)
		U[i] = 0.0;
	QX = Q * X11;
	QY = Q * Y11;
	/* C====================================== */
	/* C=====  STRIKE-SLIP CONTRIBUTION  ===== */
	/* C====================================== */
	if (DISL1 != 0.0) {
		dum1 = ALP3 * SD;
		DU[1] = -XI * QY - TT - dum1 * AI1;
		DU[2] = -Q / R + dum1 * Y / RD;
		DU[3] = Q * QY - dum1 * AI2;
		DU[4] = XI2 * Q * Y32 - dum1 * AJ1;
		DU[5] = XI * Q / R3 - dum1 * AJ2;
		DU[6] = -XI * Q2 * Y32 - dum1 * AJ3;
		DU[7] = -XI * FY - D * X11 + dum1 * (XY + AJ4);
		DU[8] = -EY + dum1 * (F1 / R + AJ5);
		DU[9] = Q * FY - dum1 * (QY - AJ6);
		DU[10] = -XI * FZ - Y * X11 + dum1 * AK1;
		DU[11] = -EZ + dum1 * Y * D11;
		DU[12] = Q * FZ + dum1 * AK2;
		dum1 = DISL1 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum1 * DU[i];
	}
	/* C====================================== */
	/* C=====    DIP-SLIP CONTRIBUTION   ===== */
	/* C====================================== */
	if (DISL2 != 0.0) {
		dum2 = ALP3 * SDCD;
		DU[1] = -Q / R + dum2 * AI3;
		DU[2] = -ET * QX - TT - dum2 * XI / RD;
		DU[3] = Q * QX + dum2 * AI4;
		DU[4] = XI * Q / R3 + dum2 * AJ4;
		DU[5] = ET * Q / R3 + QY + dum2 * AJ5;
		DU[6] = -Q2 / R3 + dum2 * AJ6;
		DU[7] = -EY + dum2 * AJ1;
		DU[8] = -ET * GY - XY * SD + dum2 * AJ2;
		DU[9] = Q * GY + dum2 * AJ3;
		DU[10] = -EZ - dum2 * AK3;
		DU[11] = -ET * GZ - XY * CD - dum2 * XI * D11;
		DU[12] = Q * GZ - dum2 * AK4;
		dum2 = DISL2 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum2 * DU[i];
	}
	/* C======================================== */
	/* C=====  TENSILE-FAULT CONTRIBUTION  ===== */
	/* C======================================== */
	if (DISL3 != 0.0) {
		dum3 = ALP3 * SDSD;
		DU[1] = Q * QY - dum3 * AI3;
		DU[2] = Q * QX + dum3 * XI / RD;
		DU[3] = ET * QX + XI * QY - TT - dum3 * AI4;
		DU[4] = -XI * Q2 * Y32 - dum3 * AJ4;
		DU[5] = -Q2 / R3 - dum3 * AJ5;
		DU[6] = Q * Q2 * Y32 - dum3 * AJ6;
		DU[7] = Q * FY - dum3 * AJ1;
		DU[8] = Q * GY - dum3 * AJ2;
		DU[9] = -Q * HY - dum3 * AJ3;
		DU[10] = Q * FZ + dum3 * AK3;
		DU[11] = Q * GZ + dum3 * XI * D11;
		DU[12] = -Q * HZ + dum3 * AK4;
		dum3 = DISL3 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum3 * DU[i];
	}
}

void UC(double XI, double ET, double Q, double Z, double DISL1, double DISL2, double DISL3, double *U, double D, double R2, double R, double XI2, double X11, double Y11, double ET2, double CD, double SD, double R3, double Y32, double R5,
		double Y, double ALP4, double ALP5, double X32, double Q2, double SDSD, double SDCD, double CDCD) {

	/* C******************************************************************** */
	/* C*****    DISPLACEMENT AND STRAIN AT DEPTH [PART-C]             ***** */
	/* C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   ***** */
	/* C******************************************************************** */
	/* C                                                                     */
	/* C***** INPUT                                                          */
	/* C*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM            */
	/* C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS            */
	/* C***** OUTPUT                                                         */
	/* C*****   U[12] : DISPLACEMENT AND THEIR DERIVATIVES                   */

	double DU[13];
	double C, H, Y0, Z0, XY, YY0, PPY, PPZ, QX, QY, QR, CQX, CDR, QQ, QQZ, QQY, Z32, Z53;
	double F1, F2, F3, PI2, X53, Y53, dum1, dum2, dum3, dum4;
	int i;

	F1 = 1.0;
	F2 = 2.0;
	F3 = 3.0;
	PI2 = 6.283185307179586;

	C = D + Z;
	X53 = (8.0 * R2 + 9.0 * R * XI + F3 * XI2) * X11 * X11 * X11 / R2;
	Y53 = (8.0 * R2 + 9.0 * R * ET + F3 * ET2) * Y11 * Y11 * Y11 / R2;
	H = Q * CD - Z;
	Z32 = SD / R3 - H * Y32;
	Z53 = F3 * SD / R5 - H * Y53;
	Y0 = Y11 - XI2 * Y32;
	Z0 = Z32 - XI2 * Z53;
	PPY = CD / R3 + Q * Y32 * SD;
	PPZ = SD / R3 - Q * Y32 * CD;
	QQ = Z * Y32 + Z32 + Z0;
	QQY = F3 * C * D / R5 - QQ * SD;
	QQZ = F3 * C * Y / R5 - QQ * CD + Q * Y32;
	XY = XI * Y11;
	QX = Q * X11;
	QY = Q * Y11;
	QR = F3 * Q / R5;
	CQX = C * Q * X53;
	CDR = (C + D) / R3;
	YY0 = Y / R3 - Y0 * CD;

	for (i = 1; i <= 12; i++)
		U[i] = 0.0;
	/* C====================================== */
	/* C=====  STRIKE-SLIP CONTRIBUTION  ===== */
	/* C====================================== */
	if (DISL1 != 0.0) {
		dum1 = ALP4 * XI;
		dum2 = ALP5 * XI;
		dum3 = ALP5 * C;
		dum4 = ALP4 * CD;
		DU[1] = dum4 * XY - dum2 * Q * Z32;
		DU[2] = ALP4 * (CD / R + F2 * QY * SD) - dum3 * Q / R3;
		DU[3] = dum4 * QY - ALP5 * (C * ET / R3 - Z * Y11 + XI2 * Z32);
		DU[4] = dum4 * Y0 - ALP5 * Q * Z0;
		DU[5] = -dum1 * (CD / R3 + F2 * Q * Y32 * SD) + dum3 * XI * QR;
		DU[6] = -dum1 * Q * Y32 * CD + dum2 * (F3 * C * ET / R5 - QQ);
		DU[7] = -dum1 * PPY * CD - dum2 * QQY;
		DU[8] = ALP4 * F2 * (D / R3 - Y0 * SD) * SD - Y / R3 * CD - ALP5 * (CDR * SD - ET / R3 - C * Y * QR);
		DU[9] = -ALP4 * Q / R3 + YY0 * SD + ALP5 * (CDR * CD + C * D * QR - (Y0 * CD + Q * Z0) * SD);
		DU[10] = dum1 * PPZ * CD - dum2 * QQZ;
		DU[11] = ALP4 * F2 * (Y / R3 - Y0 * CD) * SD + D / R3 * CD - ALP5 * (CDR * CD + C * D * QR);
		DU[12] = YY0 * CD - ALP5 * (CDR * SD - C * Y * QR - Y0 * SDSD + Q * Z0 * CD);
		dum1 = DISL1 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum1 * DU[i];
	}
	/* C====================================== */
	/* C=====    DIP-SLIP CONTRIBUTION   ===== */
	/* C====================================== */
	if (DISL2 != 0.0) {
		dum1 = ALP5 * C;
		dum2 = F2 * Q;
		DU[1] = ALP4 * CD / R - QY * SD - dum1 * Q / R3;
		DU[2] = ALP4 * Y * X11 - dum1 * ET * Q * X32;
		DU[3] = -D * X11 - XY * SD - dum1 * (X11 - Q2 * X32);
		DU[4] = -ALP4 * XI / R3 * CD + dum1 * XI * QR + XI * Q * Y32 * SD;
		DU[5] = -ALP4 * Y / R3 + dum1 * ET * QR;
		DU[6] = D / R3 - Y0 * SD + dum1 / R3 * (F1 - F3 * Q2 / R2);
		DU[7] = -ALP4 * ET / R3 + Y0 * SDSD - ALP5 * (CDR * SD - C * Y * QR);
		DU[8] = ALP4 * (X11 - Y * Y * X32) - dum1 * ((D + dum2 * CD) * X32 - Y * ET * Q * X53);
		DU[9] = XI * PPY * SD + Y * D * X32 + dum1 * ((Y + dum2 * SD) * X32 - Y * Q2 * X53);
		DU[10] = -Q / R3 + Y0 * SDCD - ALP5 * (CDR * CD + C * D * QR);
		DU[11] = ALP4 * Y * D * X32 - dum1 * ((Y - dum2 * SD) * X32 + D * ET * Q * X53);
		DU[12] = -XI * PPZ * SD + X11 - D * D * X32 - dum1 * ((D - dum2 * CD) * X32 - D * Q2 * X53);
		dum1 = DISL2 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum1 * DU[i];
	}
	/* C======================================== */
	/* C=====  TENSILE-FAULT CONTRIBUTION  ===== */
	/* C======================================== */
	if (DISL3 != 0.0) {
		dum1 = ALP5 * C;
		dum2 = F2 * Q;
		DU[1] = -ALP4 * (SD / R + QY * CD) - ALP5 * (Z * Y11 - Q2 * Z32);
		DU[2] = ALP4 * F2 * XY * SD + D * X11 - dum1 * (X11 - Q2 * X32);
		DU[3] = ALP4 * (Y * X11 + XY * CD) + ALP5 * Q * (C * ET * X32 + XI * Z32);
		DU[4] = ALP4 * XI / R3 * SD + XI * Q * Y32 * CD + ALP5 * XI * (F3 * C * ET / R5 - F2 * Z32 - Z0);
		DU[5] = ALP4 * F2 * Y0 * SD - D / R3 + dum1 / R3 * (F1 - F3 * Q2 / R2);
		DU[6] = -ALP4 * YY0 - ALP5 * (C * ET * QR - Q * Z0);
		DU[7] = ALP4 * (Q / R3 + Y0 * SDCD) + ALP5 * (Z / R3 * CD + C * D * QR - Q * Z0 * SD);
		DU[8] = -ALP4 * F2 * XI * PPY * SD - Y * D * X32 + dum1 * ((Y + dum2 * SD) * X32 - Y * Q2 * X53);
		DU[9] = -ALP4 * (XI * PPY * CD - X11 + Y * Y * X32) + ALP5 * (C * ((D + dum2 * CD) * X32 - Y * ET * Q * X53) + XI * QQY);
		DU[10] = -ET / R3 + Y0 * CDCD - ALP5 * (Z / R3 * SD - C * Y * QR - Y0 * SDSD + Q * Z0 * CD);
		DU[11] = ALP4 * F2 * XI * PPZ * SD - X11 + D * D * X32 - dum1 * ((D - dum2 * CD) * X32 - D * Q2 * X53);
		DU[12] = ALP4 * (XI * PPZ * CD + Y * D * X32) + ALP5 * (C * ((Y - dum2 * SD) * X32 + D * ET * Q * X53) + XI * QQZ);

		dum1 = DISL3 / PI2;
		for (i = 1; i <= 12; i++)
			U[i] = U[i] + dum1 * DU[i];
	}
}

void DCCON0(double ALPHA, double DIP, double *ALP1, double *ALP2, double *ALP3, double *ALP4, double *ALP5, double *SD, double *CD, double *SDSD, double *CDCD, double *SDCD, double *S2D, double *C2D) {

	/* C*******************************************************************  */
	/* C*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****  */
	/* C*******************************************************************  */
	/* C                                                                     */
	/* C***** INPUT                                                          */
	/* C*****   ALPHA : MEDIUM CONSTANT  [LAMBDA+MYU]/[LAMBDA+2*MYU]          */
	/* C*****   DIP   : DIP-ANGLE [DEGREE]                                    */
	/* C### CAUTION ### IF COS[DIP] IS SUFFICIENTLY SMALL, IT IS SET TO ZERO  */

	double F1, F2, PI2, EPS, P18;

	F1 = 1.0;
	F2 = 2.0;
	PI2 = 6.283185307179586;
	EPS = 1e-6;

	*ALP1 = (F1 - ALPHA) / F2;
	*ALP2 = ALPHA / F2;
	*ALP3 = (F1 - ALPHA) / ALPHA;
	*ALP4 = F1 - ALPHA;
	*ALP5 = ALPHA;

	P18 = PI2 / 360.0;
	*SD = sin(DIP * P18);
	*CD = cos(DIP * P18);
	if (fabs(*CD) < EPS) {
		*CD = 0.0;
		if (*SD > 0.0) *SD = F1;
		if (*SD < 0.0) *SD = -F1;
	}
	*SDSD = (*SD) * (*SD);
	*CDCD = (*CD) * (*CD);
	*SDCD = (*SD) * (*CD);
	*S2D = F2 * (*SDCD);
	*C2D = (*CDCD) - (*SDSD);
}

void DCCON2(double XI, double ET, double Q, double KXI, double KET, double SD, double CD, double *XI2, double *ET2, double *Q2, double *R2, double *R, double *R3, double *R5, double *Y, double *D, double *TT, double *ALX, double *X11,
		double *X32, double *ALE, double *Y11, double *Y32, double *EY, double *EZ, double *FY, double *FZ, double *GY, double *GZ, double *HY, double *HZ) {

	/* C**********************************************************************  */
	/* C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   *****  */
	/* C**********************************************************************  */
	/* C                                                                        */
	/* C***** INPUT                                                             */
	/* C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                   */
	/* C*****   SD,CD   : SIN, COS OF DIP-ANGLE                                 */
	/* C*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY   */
	/* C                                                                        */
	/* C### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER0 */

	double RXI, RET;
	double F1, F2, EPS;

	F1 = 1.0;
	F2 = 2.0;
	EPS = 1e-6;

	if (fabs(XI) < EPS) XI = 0.0;
	if (fabs(ET) < EPS) ET = 0.0;
	if (fabs(Q) < EPS) Q = 0.0;
	*XI2 = XI * XI;
	*ET2 = ET * ET;
	*Q2 = Q * Q;
	*R2 = (*XI2) + (*ET2) + (*Q2);
	*R = sqrt(*R2);
	if (*R == 0.0) return;
	*R3 = (*R) * (*R2);
	*R5 = (*R3) * (*R2);
	*Y = ET * CD + Q * SD;
	*D = ET * SD - Q * CD;
	if (Q == 0.0) *TT = 0.0;
	else *TT = atan(XI * ET / (Q * (*R)));

	if (KXI == 1) {
		*ALX = -log((*R) - XI);
		*X11 = 0.0;
		*X32 = 0.0;
	}
	else {
		RXI = (*R) + XI;
		*ALX = log(RXI);
		*X11 = F1 / ((*R) * RXI);
		*X32 = ((*R) + RXI) * (*X11) * (*X11) / (*R);
	}
	if (KET == 1) {
		*ALE = -log((*R) - ET);
		*Y11 = 0.0;
		*Y32 = 0.0;
	}
	else {
		RET = (*R) + ET;
		*ALE = log(RET);
		*Y11 = F1 / ((*R) * RET);
		*Y32 = ((*R) + RET) * (*Y11) * (*Y11) / (*R);
	}
	*EY = SD / (*R) - (*Y) * Q / (*R3);
	*EZ = CD / (*R) + (*D) * Q / (*R3);
	*FY = (*D) / (*R3) + (*XI2) * (*Y32) * SD;
	*FZ = (*Y) / (*R3) + (*XI2) * (*Y32) * CD;
	*GY = F2 * (*X11) * SD - (*Y) * Q * (*X32);
	*GZ = F2 * (*X11) * CD + (*D) * Q * (*X32);
	*HY = (*D) * Q * (*X32) + XI * Q * (*Y32) * SD;
	*HZ = (*Y) * Q * (*X32) + XI * Q * (*Y32) * CD;
}

void DC3D(double ALPHA, double X, double YY, double Z, double DEPTH, double DIP, double AL1, double AL2, double AW1, double AW2, double DISL1, double DISL2, double DISL3, double *UX, double *UY, double *UZ, double *UXX, double *UYX,
		double *UZX, double *UXY, double *UYY, double *UZY, double *UXZ, double *UYZ, double *UZZ, int *IRET) {

//for geometry, see	http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html

	double XI[3], ET[3], KXI[3], KET[3], U[13], DU[13], DUA[13], DUB[13], DUC[13];
	double EPS, AALPHA, DDIP, ZZ, DD1, DD2, DD3, R12, R21, R22, P, Q, Y;
	int i, j, k;
	double ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD, S2D, C2D;
	double XI2, ET2, Q2, R, R2, R3, R5, D, TT, ALX, ALE, X11, Y11, X32, Y32, EY, EZ, FY, FZ, GY, GZ, HY, HZ;
	static int warning_notprintedyet=1;

	EPS = 1e-6;

	if (Z > 0.0) {
		if (warning_notprintedyet) {
			print_logfile("** Warning: POSITIVE Z WAS GIVEN IN SUB-DC3D. **\n");
			warning_notprintedyet=0;
		}
	}
	for (i = 1; i <= 12; i++) {
		U[i] = 0.0;
		DUA[i] = 0.0;
		DUB[i] = 0.0;
		DUC[i] = 0.0;
	}
	AALPHA = ALPHA;
	DDIP = DIP;

	DCCON0(AALPHA, DDIP, &ALP1, &ALP2, &ALP3, &ALP4, &ALP5, &SD, &CD, &SDSD, &CDCD, &SDCD, &S2D, &C2D);

	ZZ = Z;
	DD1 = DISL1;
	DD2 = DISL2;
	DD3 = DISL3;
	XI[1] = X - AL1;
	XI[2] = X - AL2;
	if (fabs(XI[1]) < EPS) XI[1] = 0.0;
	if (fabs(XI[2]) < EPS) XI[2] = 0.0;
	/* C======================================     */
	/* C=====  REAL-SOURCE CONTRIBUTION  =====     */
	/* C======================================     */
	D = DEPTH + Z;
	P = YY * CD + D * SD;
	Q = YY * SD - D * CD;
	ET[1] = P - AW1;
	ET[2] = P - AW2;
	if (fabs(Q) < EPS) Q = 0.0;
	if (fabs(ET[1]) < EPS) ET[1] = 0.0;
	if (fabs(ET[2]) < EPS) ET[2] = 0.0;
	/* C--------------------------------  */
	/* C----- REJECT SINGULAR CASE -----  */
	/* C--------------------------------  */
	/* C----- ON FAULT EDGE               */
	if (Q == 0.0 && ((XI[1] * XI[2] <= 0.0 && ET[1] * ET[2] == 0.0) || (ET[1] * ET[2] <= 0.0 && XI[1] * XI[2] == 0.0))) goto raus;
	/* C----- ON NEGATIVE EXTENSION OF FAULT EDGE  */
	KXI[1] = 0;
	KXI[2] = 0;
	KET[1] = 0;
	KET[2] = 0;
	R12 = sqrt(XI[1] * XI[1] + ET[2] * ET[2] + Q * Q);
	R21 = sqrt(XI[2] * XI[2] + ET[1] * ET[1] + Q * Q);
	R22 = sqrt(XI[2] * XI[2] + ET[2] * ET[2] + Q * Q);
	if (XI[1] < 0.0 && R21 + XI[2] < EPS) KXI[1] = 1;
	if (XI[1] < 0.0 && R22 + XI[2] < EPS) KXI[2] = 1;
	if (ET[1] < 0.0 && R12 + ET[2] < EPS) KET[1] = 1;
	if (ET[1] < 0.0 && R22 + ET[2] < EPS) KET[2] = 1;

	for (k = 1; k <= 2; k++) {
		for (j = 1; j <= 2; j++) {
			DCCON2(XI[j], ET[k], Q, KXI[k], KET[j], SD, CD, &XI2, &ET2, &Q2, &R2, &R, &R3, &R5, &Y, &D, &TT, &ALX, &X11, &X32, &ALE, &Y11, &Y32, &EY, &EZ, &FY, &FZ, &GY, &GZ, &HY, &HZ);
			UA(XI[j], ET[k], Q, DD1, DD2, DD3, DUA, Y11, X11, ALP2, ALP1, TT, R, ALE, XI2, Y32, Q2, SD, R3, FY, D, EY, CD, FZ, Y, EZ, ALX, GY, GZ, HY, HZ);

			for (i = 1; i <= 10; i += 3) {
				DU[i] = -DUA[i];
				DU[i + 1] = -DUA[i + 1] * CD + DUA[i + 2] * SD;
				DU[i + 2] = -DUA[i + 1] * SD - DUA[i + 2] * CD;
				if (i == 10) {
					DU[i] = -DU[i];
					DU[i + 1] = -DU[i + 1];
					DU[i + 2] = -DU[i + 2];
				}
			}
			for (i = 1; i <= 12; i++) {
				if (j + k != 3) U[i] = U[i] + DU[i];
				if (j + k == 3) U[i] = U[i] - DU[i];
			}
		}
	}
	/* C======================================= */
	/* C=====  IMAGE-SOURCE CONTRIBUTION  ===== */
	/* C======================================= */
	D = DEPTH - Z;
	P = YY * CD + D * SD;
	Q = YY * SD - D * CD;
	ET[1] = P - AW1;
	ET[2] = P - AW2;
	if (fabs(Q) < EPS) Q = 0.0;
	if (fabs(ET[1]) < EPS) ET[1] = 0.0;
	if (fabs(ET[2]) < EPS) ET[2] = 0.0;
	/* C-------------------------------- */
	/* C----- REJECT SINGULAR CASE ----- */
	/* C-------------------------------- */
	/* C----- ON FAULT EDGE              */
	if (Q == 0.0 && ((XI[1] * XI[2] <= 0.0 && ET[1] * ET[2] == 0.0) || (ET[1] * ET[2] <= 0.0 && XI[1] * XI[2] == 0.0))) goto raus;
	/* C----- ON NEGATIVE EXTENSION OF FAULT EDGE */
	KXI[1] = 0;
	KXI[2] = 0;
	KET[1] = 0;
	KET[2] = 0;
	R12 = sqrt(XI[1] * XI[1] + ET[2] * ET[2] + Q * Q);
	R21 = sqrt(XI[2] * XI[2] + ET[1] * ET[1] + Q * Q);
	R22 = sqrt(XI[2] * XI[2] + ET[2] * ET[2] + Q * Q);
	if (XI[1] < 0.0 && R21 + XI[2] < EPS) KXI[1] = 1;
	if (XI[1] < 0.0 && R22 + XI[2] < EPS) KXI[2] = 1;
	if (ET[1] < 0.0 && R12 + ET[2] < EPS) KET[1] = 1;
	if (ET[1] < 0.0 && R22 + ET[2] < EPS) KET[2] = 1;

	for (k = 1; k <= 2; k++) {
		for (j = 1; j <= 2; j++) {
			DCCON2(XI[j], ET[k], Q, KXI[k], KET[j], SD, CD, &XI2, &ET2, &Q2, &R2, &R, &R3, &R5, &Y, &D, &TT, &ALX, &X11, &X32, &ALE, &Y11, &Y32, &EY, &EZ, &FY, &FZ, &GY, &GZ, &HY, &HZ);
			UA(XI[j], ET[k], Q, DD1, DD2, DD3, DUA, Y11, X11, ALP2, ALP1, TT, R, ALE, XI2, Y32, Q2, SD, R3, FY, D, EY, CD, FZ, Y, EZ, ALX, GY, GZ, HY, HZ);
			UB(XI[j], ET[k], Q, DD1, DD2, DD3, DUB, R, D, Y, CD, XI2, Q2, CDCD, SDCD, SD, ALE, Y11, X11, ALP3, TT, Y32, R3, FY, EY, FZ, EZ, GY, GZ, HY, HZ, SDSD);
			UC(XI[j], ET[k], Q, ZZ, DD1, DD2, DD3, DUC, D, R2, R, XI2, X11, Y11, ET2, CD, SD, R3, Y32, R5, Y, ALP4, ALP5, X32, Q2, SDSD, SDCD, CDCD);

			for (i = 1; i <= 10; i += 3) {
				DU[i] = DUA[i] + DUB[i] + Z * DUC[i];
				DU[i + 1] = (DUA[i + 1] + DUB[i + 1] + Z * DUC[i + 1]) * CD - (DUA[i + 2] + DUB[i + 2] + Z * DUC[i + 2]) * SD;
				DU[i + 2] = (DUA[i + 1] + DUB[i + 1] - Z * DUC[i + 1]) * SD + (DUA[i + 2] + DUB[i + 2] - Z * DUC[i + 2]) * CD;
				if (i == 10) {
					DU[10] = DU[10] + DUC[1];
					DU[11] = DU[11] + DUC[2] * CD - DUC[3] * SD;
					DU[12] = DU[12] - DUC[2] * SD - DUC[3] * CD;
				}
			}
			for (i = 1; i <= 12; i++) {
				if ((j + k) != 3) U[i] = U[i] + DU[i];
				if ((j + k) == 3) U[i] = U[i] - DU[i];
			}
		}
	}

	*UX = U[1];
	*UY = U[2];
	*UZ = U[3];
	*UXX = U[4];
	*UYX = U[5];
	*UZX = U[6];
	*UXY = U[7];
	*UYY = U[8];
	*UZY = U[9];
	*UXZ = U[10];
	*UYZ = U[11];
	*UZZ = U[12];
	*IRET = 0;

	return;

	/* C=========================================== */
	/* C=====  IN CASE OF SINGULAR [ON EDGE]  ===== */
	/* C=========================================== */

	raus: *UX = 0.0;
	*UY = 0.0;
	*UZ = 0.0;
	*UXX = 0.0;
	*UYX = 0.0;
	*UZX = 0.0;
	*UXY = 0.0;
	*UYY = 0.0;
	*UZY = 0.0;
	*UXZ = 0.0;
	*UYZ = 0.0;
	*UZZ = 0.0;
	*IRET = 1;

	return;

}

