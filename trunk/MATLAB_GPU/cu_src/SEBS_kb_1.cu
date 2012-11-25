#include "SEBS_kb_1.cuh"

__device__ float Psim(float Zref_m, float z0m,float L)
{
	if (L<0.0f)
	{
		/* unstable */
		float xx2 = sqrtf(sqrtf(( 1.0f-19.0f*Zref_m/L ))); 
		float xx1 = sqrtf(sqrtf(( 1.0f-19.0f*z0m/L )));
		return( 2.0f*logf((1.0f+xx2)/(1.0f+xx1))+logf((1.0f+xx2*xx2)/(1.0f+xx1*xx1))-2.0f*atanf(xx2)+2.0f*atanf(xx1) );
	}
	else
	{
		/* stable */
		return( fmaxf( -5.3f , -5.3f*(Zref_m-z0m)/L ) ); 
	}
}

__device__ float Psih(float z0h,float Zref_h,float L)
{
	if (L<0.0f)
	{
		/* Unstable */
		return( 2.0f*logf((1.0f+sqrtf(1.0f-11.6f*Zref_h/L))/(1.0f+sqrtf(1.0f-11.6f*z0h/L))) ); 
	}
	else
	{
		/* stable */
		return( fmaxf( -8.0f , -8.0f*(Zref_h-z0h)/L ) );
	}
}


__device__ void SEBS_kb_1_yang_bulktransfercoefficients(
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
	
	float log_z_z0m = logf(Zref_m/z0m);
	float log_z_z0h = logf(Zref_h/z0h);
	
	float uprf,ptprf;
	if (L<0.0f)
	{
		/* Unstable case */
		uprf = fmaxf( log_z_z0m-psim , 0.50f*log_z_z0m);
		ptprf = fmaxf( log_z_z0h-psih , 0.33f*log_z_z0h); 
	}
	else
	{
		/* Stable case */
		uprf = fminf( log_z_z0m-psim , 2.0f*log_z_z0m );
		ptprf = fminf( log_z_z0h-psih , 2.0f*log_z_z0h );
	}
	*C_D = k*k/( uprf*uprf ); /* eq 160a p, 1651 */
	
	*C_H = k*k/Pr_u/( uprf*ptprf); /* eq 160b p, 1651 */
}

__device__ float PFunc(
			float alpha, 
			float beta, 
			float VarGamma)
{
	float tmp1=logf(-alpha);
	float tmp2=logf(beta);
	float tmp3=logf(VarGamma);
	return (fmaxf(0.0f,
	            	0.03728f - 
	   				0.093143f*tmp1 + 0.017131f*tmp1*tmp1 -
				   	0.240690f*tmp2 - 0.084598f*tmp2*tmp2 +
	   				0.306160f*tmp3 - 0.125870f*tmp3*tmp3 +
	  	 			0.037666f*tmp1*tmp2 - 0.016498f*tmp1*tmp3 + 0.1828f*tmp2*tmp3
	  	 	   )
	  	   );
}

__device__ float SEBS_kb_1_yang_MOlength(
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
	float zeta=0.0f;
		
	/*  Bulk Richardson number eq 3, p1649 (T0 is replace by Theta_a) */
	float Ri = g*(Zref_m-z0m)/(Uref*2.0f)*(Theta_a-Theta_s)/Theta_a;  
	
	float log_z_z0m  = logf(Zref_m/z0m);
	float log_z_z0h  = logf(Zref_h/z0h);
	
	
	/* Calcualtions */
	
	/*  Unstable case: calculated by an aproximate analytical solution proposed by Yang (2001,JAM) */	
	
	if (Ri<0.0f)
	{
		float Ri_u = fmaxf(Ri,-10.0f);
		float numerator = (Ri_u/Pr_u)*(log_z_z0m*log_z_z0m / log_z_z0h ) * (Zref_m / (Zref_m-z0m )); /* Part of eq 13 */
		
		float p = PFunc( Ri_u/Pr_u, log_z_z0m, log_z_z0h);
		float denominator = 1.0f- ((Ri_u/Pr_u) * 3.890086f*(Zref_m-z0m)/(Zref_h-z0h))*p;
		zeta= numerator/denominator; /*  eq 13, p1650 */
	}
	else
	{
		/*  Stable case: calculated by the exact analytical solution proposed by Yang (2001,JAM) */ 
		float Ri_s2 = fminf( 0.2f , Pr_s*8.0f*(1.0f-z0h/Zref_h)/ 28.09f /(1.0f-z0m/Zref_m) - 0.05f);
		float Ri_s = fminf( Ri , Ri_s2 );
		
		float a=(Ri_s/Pr_s)*28.09f*(1.0f-z0m/Zref_m)*(1.0f-z0m/Zref_m)-8.0f*(1.0f-z0m/Zref_m)*(1.0f-z0h/Zref_h);  /*  eq 9a, p1649 */
		float b=(2.0f*(Ri_s/Pr_s)*5.3f*log_z_z0m-log_z_z0h)*(1.0f-z0m/Zref_m); /*  eq 9b, p1649 */
		float c=(Ri_s/Pr_s)*log_z_z0m*log_z_z0m; /* eq 9c, p1649 */
		zeta = (-b - sqrtf(b*b-4.0f*a*c))/(2.0f*a); /* eq 8, p1649 */
	}
	
	
	/* Monin Obukhov length */
	if (zeta==0.0f)
	{
		zeta=1e-6f;
	}
	float L  = Zref_m/zeta; /* p1648. */
	return(L);
}




__device__ float SEBS_kb_1_yang(
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
	float Theta_a = Tref_K *powf((float)(P0/Pref ),0.286f); /* potential Air temperature [K] */
	float Theta_s = LST_K *powf((float)(P0/Ps ),0.286f); /* potential surface temperature [K] */
	
	float c_u = k/logf(Zref_m/z0m);
	float c_pt = k/logf(Zref_h/z0m);
	float Thetastar = c_pt*(Theta_a - Theta_s);
	float ustar = c_u * Uref; 
	
	float Nu = 1.328e-5f * (P0/Pref) * powf((float)(Theta_a/T0),1.754f); /* viscosity of air */

	/* iteration over z0h */
	int i;
	float z0h,L,C_D,C_H;
	for(i=0;i<3;i++)
	{
		z0h = fminf( Zref_h/10.0f , 
	               fmaxf( 1.0E-10f,
	                    70.0f * Nu / ustar * expf(-7.2f*sqrt(ustar)* powf((float)fabsf(-Thetastar),0.25f))
	                  )
	             );
		L = SEBS_kb_1_yang_MOlength(Zref_m,Zref_h,z0m,z0h,Uref,Theta_s,Theta_a);
		SEBS_kb_1_yang_bulktransfercoefficients(L,z0m,z0h,Zref_m,Zref_h,&C_D,&C_H);
		ustar = sqrtf(C_D*Uref); /* eq 16a, p1651 update ustar */
		Thetastar = C_H*(Theta_a-Theta_s)*rsqrtf(C_D); /* eq 16b, p1651 update Thetastar */
	}	
	return(z0h);
}


__device__ void SEBS_kb_1_G(
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
 	float fs = 1.0f - fc; /* fractional soil coverage */
	float Nu = 1.327e-5f * (P0/Pref) * powf((Tref_K / T0),1.81f); /* Kinematic viscosity of air (Massman 1999b) (10 cm^2/s) */
	
	float ust2u_h = 0.320f -0.264f * expf(-15.1f * Cd * LAI); /* Ratio of ustar and u(h) (Su. 2001 Eq. 8) */
	float Cd_bs = 2.0f*ust2u_h*ust2u_h; /* Bulk surface drag cofficient (Su 2001) */
	
	/* within-canopy wind speed profile extinction coefficient */
	float n_ec = Cd * LAI / (Cd_bs); /* windspeed profile extinction coefficient (Su. 2002 Eq 7.) */
	int In_ec = (n_ec!=0.0f);
	
	/* Ratio of displacement height and canopy height (derived from Su 2002, eq 9) */
	float d2h=In_ec? 1.0f - 1.0f/(2.0f*n_ec)*(1.0f - expf(-2.0f * n_ec)) : 0.0f;
	
	/* Displacement height */
	*d0 = d2h * hc;
	
	/* Roughness length for momentum */
	*z0m = 0.005f; /* roughness height for bare soil */
	if (fc>0.0f)
	{
		/* roughness height for vegetation (Eq 10 Su. 2001) */
		*z0m = hc*(1.0f-d2h)*expf(-k/ust2u_h); 
	}

	/* KB-1 for Full canopy only (Choudhury and Monteith, 1988) */
	float kB1_c =In_ec? k*Cd/(4.0f*Ct*ust2u_h*(1.0f-expf(-n_ec*0.5f))) : 0.0f;
	
	/* KB-1 mixed soil and canopy */
	float u_h = fmaxf(0.0f,Uref*logf((hc-*d0)/ *z0m)/logf((Zref-*d0)/ *z0m)); /* (within canopy) Horizontal windspeed at the top of the canopy */
	float ustar_m = ust2u_h * u_h; /* friction velocity */
	/* note that if Zref becomes smaller than d0, ustar becomes on-physical */
	
	float Re_m = (Nu != 0.0f)? ustar_m * hc / Nu : 0.0f;
	
	float Ct_m = powf(Pr,-2.0f/3.0f) * rsqrtf(Re_m); /* heat transfer coefficient for soil */
	float kB1_m = (k * ust2u_h) * (*z0m / hc) / Ct_m; /* KB-1 for Mixed canopy (Choudhury and Monteith, 1988) */
	
	/* displacement height for soil = 0 */
	float d0_s = 0.000f; 
	/* momentum roughness parameter (0.009 ~ 0.024)(Su et al., 1997, IJRS, p.2105-2124.) */
	float hs = 0.009f; 
	/* Friction velocity in case of soil only. (Brutsaert 2008, Eq 2.41 P43 )[m/2] */
	float Ustar_s = Uref * k / logf((Zref - d0_s)/ hs); 
	
	float Re_s=(Nu != 0.0f)? Ustar_s* hs / Nu: 0.0f;
	
	float kB1_s = 2.46f * sqrtf(sqrtf((float) Re_s)) -2.00148f; /* KB-1 for Soil only (Brutsaert,1982) */
	
	
	/* KB-1 All */
	float kB_1 = (fc*fc) * kB1_c +
	              2.0f*(fc*fs) * kB1_m +
	              (fs*fs) * kB1_s; 
	 
	/* roughness length for Heat */
	*z0h = *z0m / expf(kB_1); /* roughness height for heat (su 2002 eq 8) */
	
	/* extra parameterization for bare soil and snow, according to Yang */
	if (fc==0.0f)
	{
		*z0h = SEBS_kb_1_yang(*z0m,Zref,Uref,LST_K,Tref_K,qa_ref, Pref,P0,Ps); 
	}
}

__global__ void SEBS_kb_1_Kernel(
			float * z0m,
			float * d0,
			float * z0h,
			float * fc,
			float * LAI,
			float * NDVI,
			float * LST_K,
			float * hc,
			float * Zref,
			float * Uref,
			float * Pref,
			float * P0,
			float * Ps,
			float * Tref_K,
			float * qa_ref,
			float * ComptMask,
			int DataSize)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	__shared__ float sfc;
	__shared__ float sLAI;
	__shared__ float sNDVI;
	__shared__ float sLST_K;
	__shared__ float shc;
	__shared__ float sZref;
	__shared__ float sUref;
	__shared__ float sPref;
	__shared__ float sP0;
	__shared__ float sPs;
	__shared__ float sTref_K;
	__shared__ float sqa_ref;

	if (idx<DataSize && ComptMask[idx]==1.0f)
	{
		sfc=fc[idx];
		sLAI=LAI[idx];
		sNDVI=NDVI[idx];
		sLST_K=LST_K[idx];
		shc=hc[idx];
		sZref=Zref[idx];
		sUref=Uref[idx];
		sPref=Pref[idx];
		sP0=P0[idx];
		sPs=Ps[idx];
		sTref_K=Tref_K[idx];
		sqa_ref=qa_ref[idx];
	
		SEBS_kb_1_G(
					sfc,
					sLAI,
					sNDVI,
					sLST_K,
					shc,
					sZref,
					sUref,
					sPref,
					sP0,
					sPs,
					sTref_K,
					sqa_ref,
					( z0m+idx ),
					( d0+idx ),
					( z0h+idx ));
	}
	else
	{
		*( z0m+idx )=-9999.0f;
		*( d0+idx )=-9999.0f;
		*( z0h+idx )=-9999.0f;
		
	}
}


