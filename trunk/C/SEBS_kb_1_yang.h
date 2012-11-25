/* SEBS_kb_1_yang Header */

#include "Constants.h"
#include "AUXFunc.h"
#include <math.h>

#include "SEBS_kb_1_yang_MOlength.h"
#include "SEBS_kb_1_yang_bulktransfercoefficients.h"

float z0mz0h(float Zref,float Nu,float ustar,float Thetastar);

float SEBS_kb_1_yang(
		float z0m,
		float Zref,
		float Uref,
		float LST_K,
		float Tref_K,
		float qa_ref,
		float Pref,
		float P0,
		float Ps);
		

