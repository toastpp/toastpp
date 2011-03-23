#include <mathlib.h>
#include <felib.h>
#include "phaseFunc.h"

using namespace toast;
toast::complex phaseFunc(const double g, const double costheta){
	double num = 1-g*g;
	double denom = pow(1 + g*g - 2*g*costheta, 1.5);
	return(toast::complex(num/(4*M_PI*denom),0));
}
