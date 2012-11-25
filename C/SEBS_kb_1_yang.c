#include "SEBS_kb_1_yang.h"

float z0mz0h(float Zref,float Nu,float ustar,float Thetastar)
{
/*	% Developed by River and Environmental Engineering Laboratory, University of Tokyo
	% input
	% Zref ! reference level of air temperature (m)
	% z0m ! aerodynamic roughness length
	% ustar ! friction velocity
	% Thetastar ! friction temperature (=-H/(rhoair*cp*ustar))
	% Nu ! kinematic viscousity 
	% 
	% real z0h	 ! roughness length for heat
	%
	% 2002, Zref_hang et al, Improvement of surface Flux parametrizations with a turbulence-related length
*/
	float beta = 7.2; /* calibration coefficient [m^-1/2 s^1/2 K^1/4] */
	float z0h = 70.0 * Nu / ustar * exp(-beta*sqrt(ustar)* pow(fabs(-Thetastar),0.25) ); /* eq 17, p2081. */
	z0h = min(Zref/10.0,max(z0h,1.0E-10));
	
	return(z0h);
}

float SEBS_kb_1_yang(
		float z0m,
		float Zref,
		float Uref,
		float LST_K,
		float Tref_K,
		float qa_ref,
		float Pref,
		float P0,
		float Ps)
{
	/* use yang's model (2002,QJRMS) */
	float Zref_m = Zref;
	float Zref_h = Zref;

	/* calculate parameters */
	float Theta_a = Tref_K *pow((P0/Pref ),0.286); /* potential Air temperature [K] */
	float Theta_s = LST_K *pow((P0/Ps ),0.286); /* potential surface temperature [K] */
	
	float c_u = k/log(Zref_m/z0m);
	float c_pt = k/log(Zref_h/z0m);
	float Thetastar = c_pt*(Theta_a - Theta_s);
	float ustar = c_u * Uref; 
	
	float Nu = 1.328e-5 * (P0/Pref) * pow(Theta_a/T0,1.754); /* viscosity of air */

	/* iteration over z0h */
	int i;
	float z0h,L,C_D,C_H;
	for(i=0;i<3;i++)
	{
		z0h = z0mz0h(Zref_h,Nu,ustar,Thetastar);
		L = SEBS_kb_1_yang_MOlength(Zref_m,Zref_h,z0m,z0h,Uref,Theta_s,Theta_a);
		SEBS_kb_1_yang_bulktransfercoefficients(L,z0m,z0h,Zref_m,Zref_h,&C_D,&C_H);
		ustar = sqrt(C_D*Uref); /* eq 16a, p1651 update ustar */
		Thetastar = C_H*(Theta_a-Theta_s)/sqrt(C_D); /* eq 16b, p1651 update Thetastar */
	}
	
	return(z0h);
}







