/* SEBS_Daily_Evapotranspiration Header */

#include "Constants.h"

__device__ void SEBS_Daily_Evapotranspiration_G(
				float day_angle,
				float lat_rad,
				float albedo,
				float SWd24_Wm,
				float Ta_av_K,
				float evapfr,
				float Ns,
				float * Rndaily,
				float * Edaily);

__global__ void SEBS_Daily_Evapotranspiration_Kernel(
				float * Rndaily,
				float * Edaily,
				float * ComptMask,
				float * day_angle,
				float * lat_rad,
				float * albedo,
				float * SWd24_Wm,
				float * Ta_av_K,
				float * evapfr,
				float * Ns,
				int DataSize);		




