/* SEBS_kb_1 Header */

#include "Constants.h"
#include "SEBS_kb_1_yang.h"
#include "AUXFunc.h"
#include <math.h>
#include <stdbool.h>

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
		float * z0h);
		
