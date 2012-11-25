/* SEBS_EnergyBalance Header */

#include "Constants.h"
#include "AUXFunc.h"
#include <math.h>
#include <stdbool.h>

float PSIm(float zeta);

float PSIh(float zeta);

float Bw(float hpbl, float L, float z0m, float d0) ;

float Cw(float hpbl, float L, float z0m, float z0h, float d0);

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
		float * H_WL);
		
		

	