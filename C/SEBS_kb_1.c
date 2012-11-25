#include "SEBS_kb_1.h"

void SEBS_kb_1(
		float fc,
		float LAI,
		float NDVI,
		float LST_K,
		float hc,
		float Zref,
		float Uref,
		float Pref,
		float P0,
		float Ps,
		float Tref_K,
		float qa_ref,
		float * z0m,
		float * d0,
		float * z0h)
{
/*
	% output
	% z0hM Roughness height for heat transfer (m), 0.02
	% d0 Displacement height
	% z0h Roughness height for momentum transfer (m), 0.2
	
	% KB_1M by Massman, 1999, (24) This a surrogate for the full LNF model for describing the combined
	% canopy and soil boundary layer effects on kBH_1. The formulation follows Su (2001), HESS
	% zref: reference height (2hc <= zref <= hst=(0.1 ~0.15)hi, (=za))
	% Uref: wind speed at reference height
	% u_h: wind speed at canopy height
	% hc: canopy height
	% LAI: canopy total leaf area index
	% Cd: foliage drag coefficient, = 0.2
	% Ct: heat transfer coefficient
	% fc : Fractional canopy cover (-)
	% Ta: ambient air temperature (0C)
	% Pa: ambient air pressure (Pa)
	% hs: height of soil roughness obstacles 
	% Nu: Kinematic viscosity of air, (m^2/s)
	% Ta: ambient temperature of the air, (0C)
	% Pa: ambient pressure of the air, (0C)
*/
	float dummy  = 0.0; 
	float fs = 1.0 - fc; /* fractional soil coverage */
	float Nu = 1.327e-5 * (P0/Pref) * pow((Tref_K / T0),1.81); /* Kinematic viscosity of air (Massman 1999b) (10 cm^2/s) */
	
	float c1 = 0.320; /* model constants (Massman 1997) */
	float c2 = 0.264; /* model constants (Massman 1997) */
	float c3 = 15.1; /* model constants (Massman 1997) */
	float ust2u_h = c1 - c2 * exp(-c3 * Cd * LAI); /* Ratio of ustar and u(h) (Su. 2001 Eq. 8) */
	float Cd_bs = 2.0*ust2u_h*ust2u_h; /* Bulk surface drag cofficient (Su 2001) */
	
	/* within-canopy wind speed profile extinction coefficient */
	float n_ec = Cd * LAI / (Cd_bs); /* windspeed profile extinction coefficient (Su. 2002 Eq 7.) */
	bool In_ec = (n_ec!=0.0);
	
	/* Ratio of displacement height and canopy height (derived from Su 2002, eq 9) */
	float d2h=In_ec? 1.0 - 1.0/(2.0*n_ec)*(1.0 - exp(-2.0 * n_ec)) : 0.0;
	
	/* Displacement height */
	*d0 = d2h * hc;
	
	/* Roughness length for momentum */
	*z0m = 0.005; /* roughness height for bare soil */
	if (fc>0.0)
	{
		/* roughness height for vegetation (Eq 10 Su. 2001) */
		*z0m = hc*(1.0-d2h)*exp(-k/ust2u_h); 
	}

	/* KB-1 for Full canopy only (Choudhury and Monteith, 1988) */
	float kB1_c =In_ec? k*Cd/(4.0*Ct*ust2u_h*(1.0-exp(-n_ec*0.5))) : 0.0;
	
	/* KB-1 mixed soil and canopy */
	float u_h = max(0.0,Uref*log((hc-*d0)/ *z0m)/log((Zref-*d0)/ *z0m)); /* (within canopy) Horizontal windspeed at the top of the canopy */
	float ustar_m = ust2u_h * u_h; /* friction velocity */
	/* note that if Zref becomes smaller than d0, ustar becomes on-physical */
	
	float Re_m = (Nu != 0.0)? ustar_m * hc / Nu : 0.0;
	
	float Ct_m = pow(Pr,-2.0/3.0) / sqrt(Re_m); /* heat transfer coefficient for soil */
	float kB1_m = (k * ust2u_h) * (*z0m / hc) / Ct_m; /* KB-1 for Mixed canopy (Choudhury and Monteith, 1988) */
	
	/* displacement height for soil = 0 */
	float d0_s = 0.000; 
	/* momentum roughness parameter (0.009 ~ 0.024)(Su et al., 1997, IJRS, p.2105-2124.) */
	float hs = 0.009; 
	/* Friction velocity in case of soil only. (Brutsaert 2008, Eq 2.41 P43 )[m/2] */
	float Ustar_s = Uref * k / log((Zref - d0_s)/ hs); 
	
	float Re_s=(Nu != 0.0)? Ustar_s* hs / Nu: 0.0;
	
	float kB1_s = 2.46 * pow(Re_s,0.25) - log(7.4); /* KB-1 for Soil only (Brutsaert,1982) */
	
	
	/* KB-1 All */
	float kB_1 = (fc*fc) * kB1_c +
	              2.0*(fc*fs) * kB1_m +
	              (fs*fs) * kB1_s; 
	 
	/* roughness length for Heat */
	*z0h = *z0m / exp(kB_1); /* roughness height for heat (su 2002 eq 8) */
	
	/* extra parameterization for bare soil and snow, according to Yang */
	if (fc==0.0)
	{
		*z0h = SEBS_kb_1_yang(*z0m,Zref,Uref,LST_K,Tref_K,qa_ref, Pref,P0,Ps); 
	}
}


