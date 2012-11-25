#include "SEBS_kb_1_yang_bulktransfercoefficients.h"

float Psim(float Zref_m, float z0m,float L)
{
	float psim = 0.0;
	float betam = 5.3;
	float gammam = 19.0;

	if (L<0.0)
	{
		/* unstable */
		float xx2 = pow(( 1.0-gammam*Zref_m/L ),0.25); 
		float xx1 = pow(( 1.0-gammam*z0m/L ),0.25);
		psim = 2.0*log((1.0+xx2)/(1.0+xx1))+log((1.0+xx2*xx2)/(1.0+xx1*xx1))-2.0*atan(xx2)+2.0*atan(xx1);
	}
	else
	{
		/* stable */
		psim = max( -betam , -betam*(Zref_m-z0m)/L ); 
	}
	
	return(psim);
}

float Psih(float z0h,float Zref_h,float L)
{
	float psih = 0.0;
	float betah = 8.0;
	float gammah = 11.6;

	if (L<0.0)
	{
		/* Unstable */
		psih = 2.0*log((1.0+sqrt(1.0-gammah*Zref_h/L))/(1.0+sqrt(1.0-gammah*z0h/L))); 
	}
	else
	{
		/* stable */
		psih = max( -betah , -betah*(Zref_h-z0h)/L );
	}
	
	return(psih)	;
}


void SEBS_kb_1_yang_bulktransfercoefficients(
		float L,
		float z0m,
		float z0h,
		float Zref_m,
		float Zref_h,
		float * C_D,
		float * C_H)
{
/*
	% Note: 
	% PURPOSE:
	% Compute C_u and C_pt for both unstable and stable surface layer based on Monin-Obukhov similarity theory. The univeral profile form 
	% was proposed Hogstrom (1996) and the analytical solution was developed by Yang (2000)1
	% C_D: Bulk Transfer Coefficient of momentum
	% C_H: Bulk Transfer Coefficient of heat
*/	
	*C_D=L;
	*C_H=L;

	float psim = Psim(Zref_m,z0m,L);
	float psih = Psih(z0h,Zref_h,L);
	
	float log_z_z0m = log(Zref_m/z0m);
	float log_z_z0h = log(Zref_h/z0h);
	
	float uprf,ptprf;
	if (L<0.0)
	{
		/* Unstable case */
		uprf = ( log_z_z0m-psim >= 0.50*log_z_z0m ) ? log_z_z0m-psim : 0.50*log_z_z0m;
		ptprf = ( log_z_z0h-psih >= 0.33*log_z_z0h ) ? log_z_z0h-psih : 0.33*log_z_z0h; 
	}
	else
	{
		/* Stable case */
		uprf = ( log_z_z0m-psim <= 2.0*log_z_z0m ) ? log_z_z0m-psim : 2.0*log_z_z0m;
		ptprf = ( log_z_z0h-psih <= 2.0*log_z_z0h ) ? log_z_z0h-psih : 2.0*log_z_z0h;
	}
	*C_D = k*k/( uprf*uprf ); /* eq 160a p, 1651 */
	
	*C_H = k*k/Pr_u/( uprf*ptprf); /* eq 160b p, 1651 */
}

