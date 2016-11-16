#ifndef CAMERA_H
#define CAMERA_H

#include "mathlib.h"

class STOASTLIB Camera
{
    public:
	Camera()
	    {
		pixelSize = 1.0;
	    }
	Camera(	int w_, int h_,
		const RVector & pos_, 
		const RVector & x_, const RVector & y_, const RVector & z_,
		double pixSize = 1.0)
	    : w(w_), h(h_), pos(pos_), x(x_), y(y_), z(z_), pixelSize(pixSize)
	    {
	    }
	void init(int w_, int h_, double f_, 
		const RVector & pos_, 
		const RVector & x_, const RVector & y_, const RVector & z_)
	{
	    w=w_; h=h_;
	    pos=pos_; x=x_; y=y_; z=z_;
	}
	virtual void getRayVector(const double ix, const double iy, RVector & rayStart, RVector & rayDir) const = 0;
	virtual void getPixelCoords(const RVector & p, double & ix, double & iy) const = 0;
	RVector getPos() const {return pos;}
	RVector getUp() const {return y;}
	RVector getViewDir() const {return z;}
	int getImageWidth() const {return w;}
	int getImageHeight() const {return h;}
	double getAspect() const;
	virtual double getFoVy() const { return 0; }
	double getPixelSize() {return pixelSize;}

    protected:
	int w, h;
	RVector pos, x, y, z;
	double pixelSize;
};

class STOASTLIB PinholeCamera : public Camera
{
    public:
	PinholeCamera(){}
	PinholeCamera(	int w_, int h_, double f_, 
			const RVector & pos_, 
			const RVector & x_, const RVector & y_, 
			const RVector & z_,
			double pixSize = 1.0)
	    : Camera(w_, h_, pos_, x_, y_, z_, pixSize), f(f_)
	    {
	    }
	void getRayVector(const double ix, const double iy, RVector & rayStart, RVector & rayDir) const;
	void getPixelCoords(const RVector & p, double & ix, double & iy) const;
	double getFoVy() const; 
	double getFocalLength() const { return f; }
	
    protected:
	double f;	
};

class STOASTLIB OrthoCamera : public Camera
{
    public:
	OrthoCamera(){}
	OrthoCamera(int w_, int h_, double pixSize,
		const RVector & pos_, 
		const RVector & x_, const RVector & y_, const RVector & z_)
	    : Camera(w_, h_, pos_, x_, y_, z_, pixSize)
	{
	}
	void getRayVector(const double ix, const double iy, RVector & rayStart, RVector & rayDir) const;
	void getPixelCoords(const RVector & p, double & ix, double & iy) const;
};

#endif
