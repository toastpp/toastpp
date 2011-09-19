#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

RVector BoundLogScaler::Unscale (const RVector &scaled) const
{
    RVector res = xmin*xmin + xmax*xmax - 2.0*xmin*xmax - 4.0*exp(scaled);
    for (int i = 0; i < res.Dim(); i++)
	if (res[i] >= 0) res[i] = xmin[i]+xmax[i] + sqrt(res[i]);
	else xERROR ("Invalid argument for sqrt");
    return res;
}

RVector BoundLogScaler::JScale (const RVector &unscaled) const
{
    RVector s = xmin+xmax - 2.0*unscaled;
    for (int i = 0; i < s.Dim(); i++)
	if (s[i]) s[i] = (xmax[i]-unscaled[i])*(unscaled[i]-xmin[i])/s[i];
	else xERROR ("Division by zero");
    return s;
}
