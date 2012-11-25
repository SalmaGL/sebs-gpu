/* SEBS_EnergyBalance Header */

#include "Constants.h"

__device__ float PSIm(float zeta);

__device__ float PSIh(float zeta);

__device__ float Bw(float hpbl, float L, float z0m, float d0);

__device__ float Cw(float hpbl, float L, float z0m, float z0h, float d0);

__device__ void SEBS_EnergyBalance(
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
			float * H_WL);
		
__global__ void SEBS_EnergyBalance_Kernel(
			int * ComptMask,
			float * d0,
			float * z0m,
			float * z0h,
			float * fc,
			float * LAI,
			float * hc,
			float * albedo,
			float * emissivity,
			float * LST_K,
			float * NDVI,
			float * SWd,
			float * LWd,
			float * hpbl,
			float * Zref,
			float * Tref_K,
			float * Uref,
			float * qaref,
			float * Pref,
			float * Ps,
			float * P0,
			float * Rn,
			float * G0,
			float * H_i,
			float * LE,
			float * evap_fr,
			float * re_i,
			float * ustar,
			float * H_DL,
			float * H_WL,
			int DataSize);		

	