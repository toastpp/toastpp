// -*-C++-*-
// ============================================================================
// TOAST v.15
// Library: libfe   File: surface.h
//
// Declaration of class Surface and derived classes
// 2-D and 3-D surface parametrisations
// ============================================================================

#ifndef __SURFACE_H
#define __SURFACE_H

class Mesh;
class Surface;

Surface *ReadSurface (std::istream &is);

// ============================================================================
// class Surface

class Surface {
public:
    Surface() {}
    virtual ~Surface() {}

    virtual Surface *Clone () const = 0;
    // return a copy of *this

    virtual int Dimension() const = 0; 
    // surface dimension

    virtual int ParamDim() const = 0;
    // dimension of parameter vectors

    virtual void Scale (double scale) = 0;
    // Scales the surface by factor scale

    virtual void Scale (const RVector &scale) = 0;
    // Anisotropic scaling. 'scale' must have the same dimension as the surface

    virtual void Point2Param (const Point &p, RVector &prm) const = 0;
    // returns parametrisation of point p in prm
    // p should be located on the surface, otherwise the result is
    // undefined.

    virtual Point Param2Point (const RVector &param) const = 0;
    // converts surface parameter param into a point

    virtual RVector DParam (const RVector &param1, const RVector &param2)
      const { return param2 - param1; }
    // Parameter difference. Derived classes take into account wrapping etc.

    virtual double ChordDist (const RVector &param1, const RVector &param2)
	const = 0;
    // returns chord length along surface between two points given by
    // parametrisations param1 and param2

    virtual RVector Normal (const RVector &param) const = 0;
    // Outward normal of surface at point defined by param

    virtual void SetCamera (const Point &cam);
    // Set camera position for subsequent calls to RayIntersect

    virtual bool RayIntersect (const RVector &dir, Point &pi) const;
    // Calculate the intersection of the surface with a ray, given by
    // reference pos set in 'SetCamera' and direction 'dir'
    // Returns intersection closest to 'p0' in positive 'dir' in 'pi'
    // Return value is false if no intersection in positive 'dir' exists
    // 'SetCamera' must have been called

    virtual bool Inside (const Point &p) const;
    // returns true if point p is inside the volume bounded by the surface

    friend Surface *ReadSurface (std::istream &is);
    // Read a surface from stream is and return pointer to it (polymorphic)
    // Return value is NULL on error

    virtual std::ostream &put (std::ostream &os) const = 0;
    // Output the surface to a stream

    friend std::ostream &operator<< (std::ostream &os, const Surface &surf);

    friend int BndNodeParams (const Mesh &mesh, RDenseMatrix &pvecs);
    // returns in the rows of pvecs the parametrisations of all boundary
    // nodes. Return value is #rows of pvecs, or 0 if mesh doesn't have
    // a surface

protected:

    virtual std::istream &get (std::istream &is) = 0;
    // Auxiliary routine for ReadSurface to retrieve surface data once
    // an instance of the specialised surface has been created
};

// ============================================================================
// class Surface2D

class Surface2D: public Surface {
public:
    Surface2D(): Surface() {}
    int Dimension() const { return 2; }
    virtual double Circumference() const = 0;

    virtual double ChordDiff (const RVector &param1, const RVector &param2)
      const = 0;
    // signed chord distance between two points. Takes into account wrapping
    // return value > 0 if param1 'right of' (>) param2
};

// ============================================================================
// class Surface3D

class Surface3D: public Surface {
public:
    Surface3D(): Surface() {}
    int Dimension() const { return 3; }
};

// ============================================================================
// class Surface_Circle
// * Parametrisation is angle from x-axis (-Pi..Pi)

class Surface_Circle: public Surface2D {
public:
    Surface_Circle ();
    Surface_Circle (double _radius, const Point &_centre);
    Surface_Circle (const Surface_Circle &surf);

    int ParamDim() const { return 1; }
    Surface *Clone() const { return new Surface_Circle (*this); }

    double Circumference() const { return Pi2*radius; }

    void Scale (double scale);
    virtual void Scale (const RVector &scale);
    void Point2Param (const Point &p, RVector &prm) const;
    Point Param2Point (const RVector &param) const;
    RVector DParam (const RVector &param1, const RVector &param2) const;
    double ChordDiff (const RVector &param1, const RVector &param2) const;
    double ChordDist (const RVector &param1, const RVector &param2) const
    { return fabs (ChordDiff (param1, param2)); }
    RVector Normal (const RVector &param) const;

    std::ostream &put (std::ostream &os) const { return os << *this; }
    friend std::istream &operator>> (std::istream &is, Surface_Circle &surf);
    friend std::ostream &operator<< (std::ostream &os,
        const Surface_Circle &surf);

protected:
    std::istream &get (std::istream &is) { return is >> *this; }

private:
    double radius;
    Point  centre;
};

// ============================================================================
// class Surface_Sphere
// Parametrisation: polar coordinates
//     param[0] = azimuth angle (-Pi..Pi)
//     param[1] = polar angle (0..Pi) where 0 is north pole (0,0,r)

class Surface_Sphere: public Surface3D {
public:
    Surface_Sphere ();
    Surface_Sphere (double _radius, const Point &_centre);
    Surface_Sphere (const Surface_Sphere &surf);

    int ParamDim() const { return 2; }
    Surface *Clone () const { return new Surface_Sphere (*this); }

    void Scale (double scale);
    void Scale (const RVector &scale);
    void Point2Param (const Point &p, RVector &prm) const;
    Point Param2Point (const RVector &param) const;
    double ChordDist (const RVector &param1, const RVector &param2) const;
    RVector Normal (const RVector &param) const;
    void SetCamera (const Point &_cam);
    bool RayIntersect (const RVector &dir, Point &pi) const;
    bool Inside (const Point &p) const;

    std::ostream &put (std::ostream &os) const { return os << *this; }
    friend std::istream &operator>> (std::istream &is, Surface_Sphere &surf);
    friend std::ostream &operator<< (std::ostream &os,
        const Surface_Sphere &surf);

protected:
    std::istream &get (std::istream &is) { return is >> *this; }

private:
    double radius;
    Point  centre;
    RVector cam, rcam;
    double camfac;
};

// ============================================================================
// class Surface_Cone
// Parametrisation: cylinder coordinates
//     param[0] = azimuth angle (-Pi..Pi)
//     param[1] = z (distance from tip, base = 1)
// Parameters of the cone:
//     tip:    coordinates of the tip
//     dir:    direction into which the cone opens
//     sap:    semi-aperture [rad]
//     height: height from tip to base
// Note that the base plane of the cone is not included in the parametrisation

class Surface_Cone: public Surface3D {
public:
    Surface_Cone ();
    Surface_Cone (const Point &_tip, const RVector &_dir, double _sap,
        double _height);
    Surface_Cone (const Surface_Cone &surf);
  
    int ParamDim() const { return 2; }
    Surface *Clone () const { return new Surface_Cone (*this); }

    void Scale (double scale);
    virtual void Scale (const RVector &scale);
    void Point2Param (const Point &p, RVector &prm) const;
    Point Param2Point (const RVector &param) const;
    double ChordDist (const RVector &param1, const RVector &param2) const;
    RVector Normal (const RVector &param) const;

    std::ostream &put (std::ostream &os) const { return os << *this; }
    friend std::istream &operator>> (std::istream &is, Surface_Cone &surf);
    friend std::ostream &operator<< (std::ostream &os,
        const Surface_Cone &surf);

protected:
    std::istream &get (std::istream &is) { return is >> *this; }

private:
    Point tip;
    RVector dir;
    double sap;
    double height;

    void ComputeRotation (const RVector &d);
    RDenseMatrix Rot, IRot; // rotation matrices local <-> global
};

#endif // !__SURFACE_H
