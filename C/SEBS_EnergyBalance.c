#include "SEBS_EnergyBalance.h"

float PSIm(float zeta)
{
	float Y = -zeta; /* coordinate system change (Brutsaert 2008, p500) */
/*	% Integrated stability function
	% 1a. Stability correction function for momentum, eq.(16) Y=-z/L, z is the height, L the Obukhov length, 
	% both specified in the calling statement. For stable conditions we use the expressions proposed by 
	% Beljaars and Holtslag (1991)and evaluated by Van den Hurk and Holtslag (1995) can be used.
*/

	float b_s = 0.667; /* constants, p. 122, van der Hurk & Holtslag (1997) */
	float c_s = 5.0; /* constants, p. 122, van der Hurk & Holtslag (1997) */
	float d_s = 0.35; /* constants, p. 122, van der Hurk & Holtslag (1997)% QUESTION: In page 24, d_s=1 */
	
	float a_u = 0.33; /* constants, p. 49, Brutsaert(2008) */
	float b_u = 0.41; /* constants, p. 49, Brutsaert(2008) */
	
	float psim = 0.0;
	
	if (zeta > 0.0)
	{
		/* STABLE conditions (According to Beljaars & Holtslag, 1991, eq. 13) */
		psim=-(zeta+b_s*(zeta-c_s/d_s)*exp(-d_s*zeta)+b_s*c_s/d_s); 
	}
	else
	{
/*		% UNSTABLE conditions(% According to Brutsaert 2008, p50)
		% For unstable conditions Kader and Yaglom reasoned that the surface layer should be subdivided into
		% three sublayers. Brutsaert (1992, 1999) combined the functional behaviour of each of these three
		% sublayers and proposed the two regions. */
		
		float y_u= (Y<=14.509365795621073) ? Y : 14.509365795621073;
		
		float x_u = pow(y_u/a_u,1.0/3.0); 

		psim = log(a_u+y_u) -1.2299999892711639*pow(y_u,1.0/3.0) 
		      +0.1416636768559522*log((1.0+x_u)*(1.0+x_u)/(1.0-x_u+x_u*x_u)) 
		      +0.4907373718030569*atan((2.0*x_u-1.0)*0.5773502691896258) 
		      +PSI0;
	}
	
	return(psim);
}

float PSIh(float zeta)
{
/*	% Integrated stability function 1b.
	% Stability correction function for heat, eq.(17)
	% Y=-z/L, z is the height, L the Obukhov length, both specified in the calling statement.
	% For stable conditions
	% we use the expressions proposed by Beljaars and Holtslag (1991)
	% and evaluated by Van den Hurk and Holtslag (1995) can be used. */
	
	float b_s = 0.667; /* constants, p. 122, van der Hurk & Holtslag (1995) */
	float c_s = 5.0; /* constants, p. 122, van der Hurk & Holtslag (1995) */
	float d_s = 0.35; /* constants, p. 122, van der Hurk & Holtslag (1995)% QUESTION: d_s=1 in page 34 of SEBS document of Su */
	
	float c_u = 0.33; /* constants, p. 443, Brutsaert, 2008 */
	float d_u = 0.057; /* constants, p. 443, Brutsaert, 2008 */
	float n = 0.78; /* constants, p. 443, Brutsaert, 2008 */
	
	float psih = 0.0;
	
	if (zeta>0.0)
	{
		/* STABLE conditions (According to Beljaars & Holtslag, 1991 eq. 13 ) */
		psih = -( pow(1.0 + 2.0/3.0 * zeta,1.5) 
		       +b_s * (zeta - c_s / d_s)*exp(-d_s * zeta)
		       +b_s * c_s / d_s - 1.0);
	}
	else
	{
		/* UNSTABLE conditions (According to Brutsaert 2008, p50) */
		psih = ((1.0 - d_u) / n) * log((c_u + pow(-zeta,n)) / c_u); 
	}

	return(psih);
}

float Bw(float hpbl, float L, float z0m, float d0) 
{
/*	% NOTES:
	% Bulk Stability function for momentum, eq.(22), (26)
	% PBL: Height of ABL or PBL
	% L: The Obukhov length
	% z0m: Surface roughnes height for momentum
	% 
	% The ASL height
	% hst = alfa*PBL, with alfa=0.10~0.15, over moderately rough surfaces, or hst=beta*z0m, with beta= 100~150.
	% 
	% Typical values:
	% The equations describe the Free convective conditions in the mixed layer,
	% provided the top of the ABL satisfies the condition -hst > 2L.
	% 
	% For a surface with moderate roughness and PBL=1000m, alfa=0.12, this
	% is -PBL/L > 17 and -L <60 (and similar values over very rough terrain).
	% NOTE: The minus (-) sign infront of B1, B11, B22 are necessary,though not
	% clearly specified by Brutsaert, 1999. This is consistent with the integral
	% form given in (16)&(17) in which y, the variable is defined as -z/L.
	% (z0m LT (alfa/beta)*PBL): Bw = -ALOG(alfa) + PSIm(alfa*PBL/L) - PSIm(z0m/L) ;(22)
	% (z0m GE (alfa/beta)*PBL): Bw = ALOG(PBL/(beta*z0m)) + PSIm(beta*z0m/L)- PSIm(z0m/L) ;(26)
	% B0 = (alfa / beta) * hpbl;
	% B1 = -z0m / L;
	% B11 = -alfa * hpbl / L;
	% B21 = PBL / (beta * z0m);
	% B22 = -beta * z0m / L; */
	
	float alpha = 0.12; /* These are mid values as given by Brutsaert, 1999 */
	float beta = 125.0; /* These are mid values as given by Brutsaert, 1999 */
	
	float bw = 0.0;
	
	if ((-z0m / L) < 0.0)
	{
		/* STABLE conditions (Brutsaert, 1982, Eq. 4.93, p.84) */
		bw = -2.2 * log(1.0 + hpbl/L);
	}
	else
	{
		/* UNSTABLE conditions (Brutsaert 2008, p53) */
		if (z0m<((alpha/ beta) * hpbl))
		{
			/* BAS stability function for unstable conditions for moderately rough terrain */
			bw = PSIm(alpha*(hpbl-d0)/L)
			    -PSIm(z0m/L)
			    +2.1202635585518332;
//			            printf("something");
		}
		else
		{
			/* BAS stability function for unstable conditions for very rough terrain */
			bw = PSIm(beta *(z0m)/L) -
			     PSIm(z0m/L) +
				 log((hpbl-d0)/(beta*z0m));
		}		
	}
	
	return(bw);
}


float Cw(float hpbl, float L, float z0m, float z0h, float d0)
{
/*	% Bulk Stability function for heat tranfer, eq.(23), (27)
	% hpbl: Height of ABL or PBL
	% L: The Obukhov length
	% z0: Surface roughnes height for momentum
	% z0h: Surface roughnes height for height transfer
	% 
	% The ASL height
	% hst = alfa*hpbl, with alfa=0.10~0.15 over moderately rough surfaces, or hst=beta*z0, with beta= 100~150. */
	
	float alfa = 0.12; /* These are mid values as given by Brutsaert, 1999 */
	float beta = 125.0; /* These are mid values as given by Brutsaert, 1999 */
	
	float cw = 0.0;
	
	
	if ((-z0h / L) < 0.0)
	{
		/* STABLE conditions (Brutsaert, 1982, Eq. 4.93,p.84) */
		cw = -7.6 * log(1.0 + hpbl/L);
	}
	else
	{
		/* UNSTABLE conditions (Brutsaert 2008, p53) */
		if (z0m < ((alfa / beta) * hpbl))
		{
			/* BAS stability functions for stable atmosphere (mod. rough terrain) */
			cw = PSIh(alfa * (hpbl-d0)/L) 
			    -PSIh(z0h/L)
			    +2.1202635585518332;
		}
		else
		{
			/* BAS stability functions for stable atmosphere (very rough terrain) */
			cw = PSIh(beta * z0m/L) -
				 PSIh(z0h/L) +
				 log((hpbl-d0)/(beta*z0m)); 
		}
	}
	
	return(cw);
}


void SEBS_EnergyBalance(
		float d0,
		float z0m,
		float z0h,
		float fc,
		float LAI,
		float hc,
		float albedo,
		float emissivity,
		float LST_K,
		float NDVI,
		float SWd,
		float LWd,
		float hpbl,
		float Zref,
		float Tref_K,
		float Uref,
		float qaref,
		float Pref,
		float Ps,
		float P0,
		float * Rn,
		float * G0,
		float * H_i,
		float * LE,
		float * evap_fr,
		float * re_i,
		float * ustar,
		float * H_DL,
		float * H_WL)
{
/*	%% NOTES
	%---------------------------------------------------------------------
	% [outputs]=EnergyBalance(inputs)
	% syntax: outputs [ ustar, LE, LE0, G0, Rn, H_DL, H_WL, H_i, evap_fr, re_i]= EnergyBalance( d0, z0m, z0h, fc, ..., 
	% inputs:[ albedo, emissivity, LST_K, SWd, LWd,hpbl, Zref, Tref_K, Uref, Earef, qaref, Pref, Ps, ri_i)
	%---------------------------------------------------------------------
	% Description of parameters
	% Zref -> Reference height (m)
	% hpbl -> Height of the PBL (m)
	% (if not available, use 1000m)
	% d0 -> Zero plane displacement height (m)
	% z0m -> Roughness height for momentum tranfer (m)
	% z0h -> Roughness height for heat tranfer (m)
	% fc -> Fractional vegetaion cover (-)
	% Uref -> Wind speed at reference height (m/s)
	% Tref_K -> Air temperature at reference height (K)
	% Pref -> Pressure at reference height (Pa)
	% qaref -> Specific humidity at reference height (kg/kg)
	% LST_K -> Surface temperature (K)
	% Ps -> Surafce pressure (Pa)
	% SWd -> Downward Solar Radiation (Watt/m^2)
	% LWd -> Downward long wave radiation (Watt/m^2)",
	% albedo -> Albedo (-)
	% emissivity-> Emissivity of the surface (-)
	%---------------------------------------------------------------------
	% Solving of 3 equations
	% Nots: Here we start to solve the system of three equations
	% i.e. the equation of momentum transfer, the equation of heat transfer, and the equation for the stability length.
	% We use ASL functions of Brutsaert, 1999, if Zref < hst, the ASL height
	% We use BAS functions of Brutsaert, 1999, if Zref > hst, Zref <= hpbl, the PBL height.
	
	% Note: We will determine the Monin-Obukhov length using the definition
	% given by Brutsaert (1982), P.65, eqn. 5.25.
	% i.e. L = - ustar^3*rhoa/(k*g*((H/Ta*Cp)+0.61*E))
	% By using the quantity, ef=Le*E/(Rn-G0), we can write
	% H = (1-ef)*(Rn-G0) and E = ef/Le*(Rn-G0)
	% So that L =-ustar^3*rhoam/(k*g*((1-ef)/(Ta*CP)+0.61*ef/Le)*(Rn-G0))
	% From this eqn., it is obvious that L=f(ustar^3) and L=f(ef^-1)
	
	
	% LIMITING CASES: A good idea could be to take L_DL and L_WL respectively as
	% the max and min stability length.
	% This has the advantage that the feedback effects between land
	% surface and the hpbl are considered.
	% For theoretical limits, H=Rn-G0, and E=Rn-G0 respectively.
	% (Note: For wet surfaces in a dry climate, E may be bigger than
	% Rn-G0 due to negative H,
	% i.e. the Oasis effect. In this case, a small modification of
	% L_WL can be expected. This is ignored for the time being)
	% Previously, we interpolated mu_i between mu_i(0) and
	% mu_i(-1500), by means of temperature difference
	% Though other factors may aslo have influences as seen in the
	% calculation of T0Ta_l and T0Ta_u,
	% the uncertainties associated to those factors do not warrant
	% their adequate applications for this case.
	% This is consistant with the definition of resistences for the
	% limiting cases, that is, that we ignore the stable case
	% over the whole region. */
	
	float dummy = 0.0;

	/* Classification */
	bool I_snow = (NDVI<0.0) && (albedo>0.47);
	bool I_water = (NDVI==-1.0);
	
	float alpha_m = 0.11; 	/* transfer coefficient () */
	float alpha_h = 0.40; 	/* transfer coefficient (Brutsaert 1982) */
	float alpha_Ch = 0.62; /* Charnok coefficient */
	
	/* Meteorological Parameters */
	
	float Nu = 1.327e-5 * (P0 / Pref) * pow(Tref_K / T0,1.81); /* Kinematic viscosity of air (Massman 1999b) (10 cm^2/s) */
	float Earef = Pref * qaref * (Rv / Rd); /* actual vapour pressure (based on Pressure at reference height) */
	
	/* Temperatures (Brutsaert 2008, P32, Eq. 2.23 ) */
	float Theta_a = Tref_K *pow(P0/Pref ,0.286); /* potential Air temperature [K] */
	float Theta_s = LST_K *pow(P0/Ps,0.286); /* potential surface temperature [K] */
	float Theta_av = Theta_a*(1.0 + 0.61 * qaref); /* virtual potential air temperature [K] */
	
	/* air densities (Brutsaert 2008, p25 , Eq 2.4,2.6 ) */
	float rhoa_m = (Pref - 0.378 * Earef)/(Rd * Tref_K); /* density of moist air [kg/m3] */
	float rhoa_WL = (Pref - 1.000 * Earef)/(Rd * LST_K); /* density of dry air. [kg/m3] */
/*	% NOTE: rhoa_WL is only used for the wet-limit. To get a true upperlimit for the sensible heat
	% the Landsurface Temperature is used as a proxy instead of air temperature. */
	
	/* moist air density (Brutsaert 2008, p29) */
	float Cp = qaref*Cpw + (1.0-qaref)*Cpd; /* Specific heat for constant pressure [kg / m3] */
	float rhoa_m_Cp = rhoa_m * Cp; /* specific air heat capacity [J K-1 m-3] */
	float rhoa_WL_Cp = rhoa_WL * Cp; /* specific air heat capacity [J K-1 m-3] */
	
	/* saturated vapor pressure at wet limit (WL) (Campbell & Norman, 1998, p41, eq 3.8 and 3.9) */
	float LST_C = LST_K - 273.15; /* Landsurface temperature[C]. */
	float A = 611.0; /* constant [Pa] */
	float B = 17.502; /* constant [- ] */
	float C = 240.97; /* constant [C ] */
	float esat_WL = A * exp(B * LST_C/(LST_C + C)); /* saturated vapor pressure [Pa] */
	float slope_WL = B * C * esat_WL / pow(C + LST_C,2.0) ; /* Slope of saturation vapor pressure [Pa C-1] */
	
/*	% NOTE: esat_WL is only used for the wet-limit. To get a true upperlimit for the sensible heat
	% the Landsurface Temperature is used as a proxy instead of air temperature.
	%% Net Radiation */
	float SWnet = (1.0 - albedo) * SWd; /* Shortwave Net Radiation [W/m2] */
	float LWnet = emissivity*LWd - emissivity*Sigma_SB*pow(LST_K,4.0); /* Longwave Net Radiation [W/m2] */
	*Rn = SWnet+LWnet; /* Total Net Radiation [W/m2] */
	
	/* Ground Heat Flux */
	/* Kustas et al 1993 */
	C = 0.34; /* amplitude of LAI, obtained through fit with SCOPE */
	float beta = 0.46; /* extinction coefficient, obtained through fit with fit with SCOPE */
	/* Kustas, W.P., Daughtry, C.S.T. van Oevelen P.J., 
	Analatytical Treatment of Relationships between Soil heat flux/net radiation and Vegetation Indices, 
	Remote sensing of environment,46:319-330 (1993) */
	*G0 = *Rn*C*exp(-beta * LAI); 

	if (I_water || I_snow)
	{
		*G0= 0.5* *Rn; 
	}
	
	/* ASL height */
	/* hst= alfa*hpbl, with alfa=0.10~0.15 over moderately rough surfaces, or hst=beta*z0, with beta= 100~150. */
	float alfa = 0.12; /* These are mid values as given by Brutsaert,1999 */
	beta = 125.0; /* These are mid values as given by Brutsaert,1999 */
	float hst = max(alfa * hpbl, beta * z0m);
	
	
	/* U* and L (Brutsaert 2008, P47, Eq. 2.46, 2.54 and 2.55 and 2.56)
	% Initial guess: (Brutsaert 2008, p46 and p57, Eq. 2.54, 2.55 and 2.46)
	% MOS (Brutsaert 2008, p46 and p57, Eq. 2.54, 2.55 and 2.46)
	% BAS (Brutsaert 2008, p46 and p52, Eq. 2.67, 2.68 and 2.46) */
	float dTheta = Theta_s - Theta_a;
	float CH = (dTheta) * k * rhoa_m_Cp;
	float CL = -rhoa_m_Cp * Theta_av/ (k * g); 
	
	float z_d0 = Zref - d0;
	float ku = k * Uref; 
	float log_z_d0_z0m = log(z_d0 / z0m);
	float log_z_d0_z0h = log(z_d0 / z0h);
	
	/* Initial guess for u*, H and L assuming neutral stability */
	float L = dummy; /* initial L is zero for neutral condition */
	*ustar = ku / log_z_d0_z0m; /* U* in neutral condition when stability factors are zero */
    float H = CH * *ustar / log_z_d0_z0h; /* H in neutral condition when stability factors are zero */
	
	float errorH = 10.0;
	float H0 = H; /* H0 is H in neutral condition */
	int steps = 0;
	
	bool IMOS = (Zref <= hst);
	bool IBAS = Zref > hst;
	float bw = dummy;
	float cw = dummy;
		
	/* tall vegetation parameterization */
	bool Itall = ( (hc>1.0) && (fc>=0.7) );
	float l = 0.027;
	float o = 0.69;
	
/*	while ( (errorH>0.01) && (steps<100) ) */
while	(steps<14)
/*	while ( (errorH>0.01) && (steps<1) ) */
	{
		/* Stability Function */
		L = CL * (*ustar * *ustar * *ustar)/ H; /* Obukhov Stability length */
				
		/* Friction Velocity */
		if (IBAS)
		{
			bw = Bw(hpbl,L,z0m,d0);
			cw = Cw(hpbl,L,z0m,z0h,d0);
		}

		*ustar = ku/(log_z_d0_z0m - bw);
		
		/* update z0h and z0m depending on land cover */
		/* tall vegeation */
		if (Itall)
		{
			z0h = z0m / exp(52.0*sqrt(l* *ustar)/LAI - o);
			log_z_d0_z0h = log(z_d0 / z0h );
		}
		
		/* over sea - water bodies, beljaars 1994 */
		if (I_water)
		{
			z0m = alpha_m*(Nu / *ustar ) + alpha_Ch*( (*ustar * *ustar) /g); 
			z0h = alpha_h*(Nu / *ustar);
			log_z_d0_z0m = log(z_d0 / z0m);
			log_z_d0_z0h = log(z_d0 / z0h );
		}
		
		
		/* Sensible Heat Flux */
		if (IMOS)
		{
			bw = PSIm(z_d0/L) + PSIm(z0m/L); /* MOST stability functions */
			cw = PSIh(z_d0/L) + PSIh(z0h/L); 
		}

		H = CH * *ustar / (log_z_d0_z0h - cw); 
		
		/* Error */
		errorH = fabs(H0 - H);
		H0 = H;
		steps++;
	}

	/* Post iteration */
	float C_i1 = dummy;
	float C_i2 = dummy;
	
	if (IMOS)
	{
		C_i1 = PSIh(Zref/L); /* 1st Stability correction term for heat (BAS condition) */
		C_i2 = PSIh(z0h/L); /* 2nd Stability correction term for heat (BAS condition) */
	}
	
	if (IBAS)
	{
		C_i1 = Cw(Zref,L,z0m,z0h,d0); 
	}
	
	/* Sensible heat Flux
	%resistances (Su 2002, eq 17) */
	*re_i=0.0;
	if ((log_z_d0_z0h + C_i2) > C_i1)
	{
		*re_i=(log_z_d0_z0h-C_i1+C_i2)/(k * *ustar); /* Actual resistance to heat transfer [s/m] */
	}
	else
	{
		*re_i=(log_z_d0_z0h)/(k * *ustar); /* Actual resistance to heat transfer [s/m] */
	}

	*H_i = rhoa_m_Cp* (Theta_s-Theta_a)/ *re_i; /* Sensible heat flux */
	
	/* Sensible heat flux at theoretical Dry limit */
	/* Dry limit */
	*H_DL = *Rn - *G0; /* Sensible heat at Dry Limit */
	
	/* Sensible heat flux at theoretical wet limit */
	/* Dry air is assumed.. eact=0, and we need to take the density of dry air */
	float L_WL = -(*ustar * *ustar * *ustar)* rhoa_WL/(k*g*(0.61* (*Rn - *G0)/L_e)); /* Obukhov stability length at Wet Limit */
	
	/* Bulk Stability Corrections */
	float C_WL = 0.0;
	if (Zref < hst)
	{
		C_WL = PSIh(-Zref/L_WL); /* Stability correction term for heat (MOS condition) */
	}
	else
	{
		C_WL = Cw(Zref,L_WL,z0m,z0h,d0); /* Stability correction term for heat (BAS condition) */
	}
	
	/* Calculating Resistances (Su 2002, p 88, eq 18) */
	float re_WL = 0.0;
	if (log_z_d0_z0h>C_WL)
	{
		re_WL = (log_z_d0_z0h-C_WL)/(k * *ustar); /* Actual resistance to heat transfer at wet limit(s/m) */
	}
	else
	{
		re_WL = (log_z_d0_z0h)/(k * *ustar); /* Actual resistance to heat transfer at wet limit(s/m) */
	}
	
	*H_WL = ((*Rn - *G0) - (rhoa_WL_Cp/re_WL)*((esat_WL)/ gammaConst))/(1.0 + slope_WL/gammaConst); /* Sensible heat at Wet Limit [W/m2](Su 2002, p88, eq 16) */
	
	/* Evaporative fraction */
	*H_i = min(*H_i, *H_DL); /* set lower limit for sensible heat [W/m2] */
	*H_i = max(*H_i, *H_WL); /* %set upper limit for sensible heat [W/m2] */
	
	/* Relative evaporation */
	float evap_re = 0.0;
	if (*H_DL <= *H_WL)
	{
		evap_re = 1.0; /* relative evaporation for water & wet surfaces [] */
	}
	else
	{
		evap_re = 1.0-(*H_i-*H_WL)/(*H_DL-*H_WL); /* relative evaporation for land surface [] */
	}
	
	
	/* Evaporative fraction */
	*evap_fr = 0.0;
	if ((*Rn - *G0) != 0.0)
	{
		*evap_fr = evap_re*(*Rn-*G0-*H_WL)/(*Rn-*G0); /* evaporative fraction [] */
	}
	else
	{
		*evap_fr = 1.0; /* evaporative fraction upper limit [] (for negative available energy) */
	}

	/* Postprocessing for arid conditions 
	% I_a = NDVI>=0.0 & NDVI<0.2 & (Theta_s-Theta_a)>2;
	% evap_fr(I_a) = 5.0511 * NDVI(I_a)./(Theta_s(I_a)-Theta_a(I_a)) - 0.006;
	% evap_fr(I_a) = 3.1537 * NDVI(I_a)./(Theta_s(I_a)-Theta_a(I_a)) + 0.2695;  */
	
	*evap_fr = max(*evap_fr,0.0);
	*evap_fr = min(*evap_fr,1.0);
	
	/* Latent heat flux [W/m2] */
	*LE = *evap_fr * (*Rn - *G0); 
	
}





