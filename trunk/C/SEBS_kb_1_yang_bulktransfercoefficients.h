/* SEBS_kb_1_yang_bulktransfercoefficients Header */

#include "Constants.h"
#include "AUXFunc.h"
#include <math.h>

float Psim(float Zref_m, float z0m,float L);

float Psih(float z0h,float Zref_h,float L);

void SEBS_kb_1_yang_bulktransfercoefficients(
		float L,
		float z0m,
		float z0h,
		float Zref_m,
		float Zref_h,
		float * C_D,
		float * C_H);

