#include "SEBS_Daily_Evapotranspiration.cuh"

__device__ void SEBS_Daily_Evapotranspiration_G(
				float day_angle,
				float lat_rad,
				float albedo,
				float SWd24_Wm,
				float Ta_av_K,
				float evapfr,
				float Ns,
				float * Rndaily,
				float * Edaily)
{  
	/* const float DT = 60.0*60.0*24.0; */ /* seconds per day */
	
	/* Solar irradiance parameters (FAO 1998, p46/p48) (eq24, eq25 eq31, eq34) */
	float solar_decl = 0.4093f*sinf(day_angle-1.405f); /*solar declination P4.31 */
	
	float omegas = acosf(-tanf(solar_decl)*tanf(lat_rad)); /*sunset hour angle in radiance */
	float Nt  =  24.0f*omegas/pi; /* total day length */
/*	% e_es  =  1.0+0.033* cos(day_angle); % Correction factor for eccentricity of earth orbit
	% cosbeta =  cos(lat_rad).* cos(solar_decl).* sin(omegas)+... % Daily zenith angles
	%   sin(lat_rad).* sin(solar_decl).*omegas;
	% omega = pi/12; % solar time angle at midpoint of period [rad/h]
	
	%% Short wave net radiation 
	% TOA Incoming Shortwave Radiation (Brutsaert 2005, p61, p62 eq 2.78) 
	% Rin24_S0 = 440.0255*DT/1e6;  % MJ m-2 day-1 should be based on Rso
	% Rin24_Sse = Rin24_S0*e_es.*cosbeta; % MJ m-2 day-1
	% 
	% % Bottom of atmosphere downwelling shortwave radiation (Hydrology Handbook 4.7, 4.2.6)
	% as = 0.25; % cloudiness/transmissivity coefficients 
	% bs = 0.5; % cloudiness/transmissivity coefficients 
	% Rin24_St = (as+bs*Ns./Nt).*Rin24_Sse_TOA;
*/	
	float Rin24_St = SWd24_Wm*DT*1.0e-6f; /* MJ m-2 day-1 (calculated by integration of daily ECMWF) */
	
	/* Net solar radation bottom of atmosphere */
	float Rn24_Sse = Rin24_St*(1.0f-albedo); /* MJ m-2 day-1 (Hydrology Handbook, 4.7, eq 4.2.5) */
	
	/* Thermal net radiation 
	% Net BOA Long wave radiation (Hydrology Handbook 4.7 eq4.2.7, eq4.2.12 and eq4.2.9): */
	float fcl = (0.1f + 0.9f * Ns/Nt); /*Cloudiness factor */
	/* float Ns_Over_Nt=min(max((1.0-0.85*tcdc)/1.1,0.0),1.0); */
	/* float fcl = (0.1 + 0.9 * (1.0-0.85*tcdc)/1.1 );*/
	
	float Ta_av_C = Ta_av_K +kel2deg;
	float enet = -0.02f+0.261f* expf(-7.77e-4f*Ta_av_C*Ta_av_C); /* Net (daily) emissivity */
	float Rn24_L=-fcl*enet*(Sigma_SB*Ta_av_K*Ta_av_K*Ta_av_K*Ta_av_K)*(DT*1.0e-6f); /* Net Longwave radiation (in MJ m-2 day-1) */

	/* Net Radiation 
	%(Hydrology Handbook 4.2.13) */
	*Rndaily = Rn24_Sse + Rn24_L;  /* Net Radiation (MJ m-2 day-1) */
	
	/* Daily Evapotranspiration */
	*Edaily = evapfr * fmaxf(*Rndaily,0.0f) / ((2.501f-0.002361f*(Ta_av_K-273.15f))*rho_w);  /* Daily Evapotranspiration [mm day-1] */
}

__global__ void SEBS_Daily_Evapotranspiration_Kernel(
				int * ComptMask,
				float * day_angle,
				float * lat_rad,
				float * albedo,
				float * SWd24_Wm,
				float * Ta_av_K,
				float * evapfr,
				float * Ns,
				float * Rndaily,
				float * Edaily,
				int DataSize)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	__shared__ float sTa_av_K;
	if (idx<DataSize && ComptMask[idx]==1)
	{
		sTa_av_K=Ta_av_K[idx];
		SEBS_Daily_Evapotranspiration_G(
					*day_angle,
					lat_rad[idx],
					albedo[idx],
					SWd24_Wm[idx],
					sTa_av_K,
					evapfr[idx],
					Ns[idx],
					(Rndaily+idx),
					(Edaily+idx));
	}
	else
	{
		*(Rndaily+idx)=-9999.0f;
		*(Edaily+idx)=-9999.0f;
	}
}










