#include "SEBS_kb_1_yang_MOlength.h"

float PFunc(
		float alpha, 
		float beta, 
		float VarGamma)
{		
	float p;
	p= 0.03728 - 
	   0.093143*log(-alpha) + 0.017131*log(-alpha)*log(-alpha) -
	   0.240690*log(beta) - 0.084598*log(beta)*log(beta) +
	   0.306160*log(VarGamma) - 0.125870*log(VarGamma)*log(VarGamma) +
	   0.037666*log(-alpha)*log(beta) - 0.016498*log(-alpha)*log(VarGamma) + 0.1828*log(beta)*log(VarGamma) ;
	return (max(0.0,p));
}

float SEBS_kb_1_yang_MOlength(
		float Zref_m,
		float Zref_h,
		float z0m,
		float z0h,
		float Uref,
		float Theta_s,
		float Theta_a)
{
/*  
	Monin Obuhkov length
	Computes L for both unstable and stable surface layer based on Monin-Obukhov similarity theory. 
	The univeral profile form was proposed Hogstrom (1996) and the anlytical solution was 
	developed by Yang (2001)
	
	REFERENCE: Similarity analytical solution:Yang, K., Tamai, N. and Koike, T., 2001: Analytical Solution of Surface Layer Similarity Equations, J. Appl. Meteorol. 40, 
	1647-1653. Profile form: Hogstrom (1996,Boundary-Layer Meteorol.)
*/
	
	/*  Allocate Memory */
	float zeta=0.0;
	
	/*  Parameters */
	float gammam = 19.0;  /* Stability function coefficient */
	float gammah = 11.6; /* Stability function coefficient; */
	
	float betam = 5.3; /* Stability function coefficient */
	float betah = 8.0; /* Stability function coefficient */
	
	/*  Bulk Richardson number eq 3, p1649 (T0 is replace by Theta_a) */
	float Ri = g*(Zref_m-z0m)/(Uref*2)*(Theta_a-Theta_s)/Theta_a;  
	
	float log_z_z0m  = log(Zref_m/z0m);
	float log_z_z0h  = log(Zref_h/z0h);
	
	
	/* Calcualtions */
	
	/*  Unstable case: calculated by an aproximate analytical solution proposed by Yang (2001,JAM) */	
	
	if (Ri<0.0)
	{
		float Ri_u = max(Ri,-10.0);
		float numerator = (Ri_u/Pr_u)*(log_z_z0m*log_z_z0m / log_z_z0h ) * (Zref_m / (Zref_m-z0m )); /* Part of eq 13 */
		
		float p = PFunc( Ri_u/Pr_u, log_z_z0m, log_z_z0h);
		float denominator = 1.0- ((Ri_u/Pr_u) * gammam*gammam/(8.0 *gammah)*(Zref_m-z0m)/(Zref_h-z0h))*p;
		zeta= numerator/denominator; /*  eq 13, p1650 */
	}
	else
	{
		/*  Stable case: calculated by the exact analytical solution proposed by Yang (2001,JAM) */ 
		float Ri_s2 =(0.2 <= Pr_s*betah*(1.0-z0h/Zref_h)/betam*betam/(1.0-z0m/Zref_m) - 0.05) 
					 ? 0.2
					 : Pr_s*betah*(1.0-z0h/Zref_h)/betam*betam/(1.0-z0m/Zref_m) - 0.05;
		float Ri_s = (Ri <= Ri_s2) ? Ri : Ri_s2;
		
		float a=(Ri_s/Pr_s)*betam*betam*(1.0-z0m/Zref_m)*(1.0-z0m/Zref_m)-betah*(1.0-z0m/Zref_m)*(1.0-z0h/Zref_h);  /*  eq 9a, p1649 */
		float b=(2.0*(Ri_s/Pr_s)*betam*log_z_z0m-log_z_z0h)*(1.0-z0m/Zref_m); /*  eq 9b, p1649 */
		float c=(Ri_s/Pr_s)*log_z_z0m*log_z_z0m; /* eq 9c, p1649 */
		zeta = (-b - sqrt(b*b-4.0*a*c))/(2.0*a); /* eq 8, p1649 */
	}
	
	
	/* Monin Obukhov length */
	if (zeta==0.0)
	{
		zeta=1e-6;
	}
	float L  = Zref_m/zeta; /* p1648. */
	return(L);
}


