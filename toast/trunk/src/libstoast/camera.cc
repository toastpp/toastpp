#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

#include "camera.h"
#include <math.h>
#include <iostream>

using namespace std;

double Camera::getAspect() const 
{ 
    return (double)(w)/(double)(h);
}

void PinholeCamera::getRayVector(const double ix, const double iy, RVector & rayStart, RVector & rayDir) const
{
    double dx = ( (double)(ix) - (double)(w-1)*0.5 );
    double dy = ( (double)(iy) - (double)(h-1)*0.5 );
    rayDir = z*f + dx*pixelSize*x + dy*pixelSize*y;
    rayDir /= l2norm(rayDir);
    rayStart = pos;
}

double PinholeCamera::getFoVy() const 
{ 
    double fovy = atan2((double)(h)*pixelSize*0.5, f) * 360.0 / Pi;
    return fovy;
}

void PinholeCamera::getPixelCoords(const RVector & p, double & ix, double & iy) const
{
    xERROR("Not implemented");
}

void OrthoCamera::getRayVector(const double ix, const double iy, RVector & rayStart, RVector & rayDir) const
{
    double dx = ( (double)(ix) - (double)(w)*0.5 + 0.5 );
    double dy = ( (double)(iy) - (double)(h)*0.5 + 0.5);
    rayStart = pos + dx*pixelSize*x + dy*pixelSize*y;
    rayDir = z;
}

void OrthoCamera::getPixelCoords(const RVector & p, double & ix, double & iy) const
{
    ix = (double)(w-1)*0.5 + dot((p-pos), x) / pixelSize;
    iy = (double)(h-1)*0.5 + dot((p-pos), y) / pixelSize;
}
