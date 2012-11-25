/* SEBS_kb_1 Header */

#include "Constants.h"


__device__ float Psim(float Zref_m, float z0m,float L);

__device__ float Psih(float z0h,float Zref_h,float L);

__device__ void SEBS_kb_1_yang_bulktransfercoefficients(
			float L,
			float z0m,
			float z0h,
			float Zref_m,
			float Zref_h,
			float * C_D,
			float * C_H);

__device__ float PFunc(
			float alpha, 
			float beta, 
			float VarGamma);

__device__ float SEBS_kb_1_yang_MOlength(
			float Zref_m,
			float Zref_h,
			float z0m,
			float z0h,
			float Uref,
			float Theta_s,
			float Theta_a);


__device__ float SEBS_kb_1_yang(
			float z0m,
			float Zref,
			float Uref,
			float LST_K,
			float Tref_K,
			float qa_ref,
			float Pref,
			float P0,
			float Ps);

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
			float * z0h);

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
			int DataSize);		
