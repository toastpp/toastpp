#define FELIB_IMPLEMENTATION

#include "felib.h"

using namespace std;

// ============================================================================

ostream &operator<< (ostream &os, const Surface &surf)
{
    return surf.put (os);
}

Surface *ReadSurface (istream &is)
{
    char cbuf[256], surftype[256];
    Surface *surf;

    do {
        if (!is.getline (cbuf, 256)) { is.clear(); return 0; }
    } while (strncasecmp (cbuf, "Surface", 7));
    sscanf (cbuf+7, "%s", surftype);

    // This is ugly since it requires the base class to know about
    // derived classes. I don't know how else to do this

    if (!strcasecmp (surftype, "Circle")) {
        surf = new Surface_Circle ();
    } else if (!strcasecmp (surftype, "Sphere")) {
        surf = new Surface_Sphere ();
    } else if (!strcasecmp (surftype, "Cone")) {
        surf = new Surface_Cone ();
    } else {
        surf = 0; // not a known surface type
    }
    // add more surface types here

    if (surf && !surf->get(is)) { // read specs from file
        delete surf;
	surf = 0;
    }

    return surf;
}

int BndNodeParams (const Mesh &mesh, RDenseMatrix &pvecs)
{
    const Surface *surf = mesh.Boundary();
    if (!surf) return 0;
    int nbnd = mesh.nbnd();
    int pdim = surf->ParamDim();
    int n, nb;
    pvecs.New (nbnd, pdim);
    RVector prm(pdim);
    for (n = nb = 0; n < mesh.nlen(); n++) {
        if (!mesh.nlist[n].isBnd()) continue;
	surf->Point2Param (mesh.nlist[n], prm);
	pvecs.SetRow (nb++, prm);
    }
    return nbnd;
}

void Surface::SetCamera (const Point &cam)
{
    ERROR_UNDEF;
}

bool Surface::RayIntersect (const RVector &dir, Point &pi) const
{
    ERROR_UNDEF;
    return false;
}

bool Surface::Inside (const Point &p) const
{
    ERROR_UNDEF;
    return false;
}

// ============================================================================

Surface_Circle::Surface_Circle (): Surface2D ()
{
    centre.New(2); // (0,0): default
    radius = 1.0;  // default
}

Surface_Circle::Surface_Circle (double _radius, const Point &_centre)
  : Surface2D ()
{
    dASSERT(_centre.Dim() == 2, "Circle centre must be 2D point");
    centre.New (2);
    centre = _centre;
    radius = _radius;
}

Surface_Circle::Surface_Circle (const Surface_Circle &surf)
{
    centre.New (2);
    centre = surf.centre;
    radius = surf.radius;
}

void Surface_Circle::Scale (double scale)
{
    // This scales w.r.t circle centre, not coordinate origin
    radius *= scale;
}

void Surface_Circle::Scale (const RVector &scale)
{
    double EPS = 1e-8;
    dASSERT(scale.Dim() == 2, "Wrong scale dimension");
    if (fabs (scale[0]-scale[1]) < EPS) Scale (scale[0]);
    else xERROR("Anisotropic scaling for circular surfaces is not defined");
}

void Surface_Circle::Point2Param (const Point &p, RVector &prm) const
{
    dASSERT(p.Dim() == 2, "Arg 1 dimension 2 required");
    dASSERT(prm.Dim() == 1, "Arg 2 dimension 1 required");
    prm[0] = atan2 (p[1]-centre[1], p[0]-centre[0]);
}

Point Surface_Circle::Param2Point (const RVector &param) const
{
    dASSERT(param.Dim() == 1, "Parameter dimension 1 expected");
    Point p(2);
    p[0] = radius * cos(param[0]) + centre[0];
    p[1] = radius * sin(param[0]) + centre[1];
    return p;
}

RVector Surface_Circle::DParam (const RVector &param1, const RVector &param2)
    const
{
    dASSERT(param1.Dim() == 1 && param2.Dim() == 1,
	    "Invalid parameter dimensions");
    double dphi = param1[0] - param2[0];
    if      (dphi >  Pi) dphi -= 2.0*Pi;
    else if (dphi < -Pi) dphi += 2.0*Pi;
    RVector dp(1);
    dp[0] = dphi;
    return dp;
}

double Surface_Circle::ChordDiff (const RVector &param1, const RVector &param2)
    const
{
    RVector dp = DParam (param1, param2);
    return radius * dp[0];
}

RVector Surface_Circle::Normal (const RVector &param) const
{
    dASSERT(param.Dim() == 1, "Arg 1 invalid vector dimension");
    RVector nml(2);
    nml[0] = cos(param[0]);
    nml[1] = sin(param[0]);
    return nml;
}

istream &operator>> (istream &is, Surface_Circle &surf)
{
    is >> surf.centre >> surf.radius;
    return is;
}

ostream &operator<< (ostream &os, const Surface_Circle &surf)
{
    os << "Surface Circle" << endl;
    os << surf.centre << ' ' << surf.radius;
    return os;
}

// ============================================================================

Surface_Sphere::Surface_Sphere (): Surface3D ()
{
    centre.New(3); // (0,0,0): default
    radius = 1.0;  // default
}

Surface_Sphere::Surface_Sphere (double _radius, const Point &_centre)
  : Surface3D ()
{
    dASSERT(_centre.Dim() == 3, "Arg 2 invalid vector dimension");
    centre.New (3);
    centre = _centre;
    radius = _radius;
}

Surface_Sphere::Surface_Sphere (const Surface_Sphere &surf)
  : Surface3D (surf)
{
    centre.New (3);
    centre = surf.centre;
    radius = surf.radius;
}

void Surface_Sphere::Scale (double scale)
{
    // This scales w.r.t circle centre, not coordinate origin
    radius *= scale;
}

void Surface_Sphere::Scale (const RVector &scale)
{
    double EPS = 1e-8;
    dASSERT(scale.Dim() == 3, "Arg 1 wrong vector dimension");
    if (fabs (scale[0]-scale[1]) < EPS &&
	fabs (scale[0]-scale[2]) < EPS) Scale (scale[0]);
    else xERROR("Anisotropic scaling for sphere surfaces is not defined");
}

void Surface_Sphere::Point2Param (const Point &p, RVector &prm) const
{
    dASSERT(p.Dim() == 3, "Arg 1 dimension 3 required");
    dASSERT(prm.Dim() == 2, "Arg 2 dimension 2 required");
    prm[0] = atan2 (p[1]-centre[1], p[0]-centre[0]); // azimuth
    prm[1] = acos ((p[2]-centre[2])/radius);  // polar
}

Point Surface_Sphere::Param2Point (const RVector &param) const
{
    dASSERT(param.Dim() == 2, "Arg 1 wrong vector dimension");
    Point p(3);
    double rsintheta = radius * sin(param[1]);
    p[0] = rsintheta * cos(param[0]);
    p[1] = rsintheta * sin(param[0]);
    p[2] = radius * cos(param[1]);
    return p;
}

double Surface_Sphere::ChordDist (const RVector &param1, const RVector &param2)
    const
{
    dASSERT(param1.Dim() == 2, "Arg 1 wrong vector dimension");
    dASSERT(param2.Dim() == 2, "Arg 2 wrong vector dimension");

    RVector p1(3), p2(3);
    double sint1 = sin(param1[1]);
    p1[0] = sint1 * cos(param1[0]);
    p1[1] = sint1 * sin(param1[0]);
    p1[2] = cos(param1[1]);
    double sint2 = sin(param2[1]);
    p2[0] = sint2 * cos(param2[0]);
    p2[1] = sint2 * sin(param2[0]);
    p2[2] = cos(param2[1]);
    double alpha = acos (p1 & p2);
    return radius * alpha;
}

RVector Surface_Sphere::Normal (const RVector &param) const
{
    dASSERT(param.Dim() == 2, "Arg 1 invalid vector dimension");
    RVector nml(3);

    double sinth = sin(param[1]);
    nml[0] = sinth * cos(param[0]);
    nml[1] = sinth * sin(param[0]);
    nml[2] = cos(param[1]);
    return nml;
}

void Surface_Sphere::SetCamera (const Point &_cam)
{
    dASSERT(_cam.Dim() == 3, "Arg 1 wrong vector dimension");
    cam.New(3), rcam.New(3);
    cam = _cam;
    rcam = cam-centre; // cam pos in relative coords
    camfac = rcam[0]*rcam[0] + rcam[1]*rcam[1] + rcam[2]*rcam[2] -
        radius*radius;
    camfac *= 4.0;
}

bool Surface_Sphere::RayIntersect (const RVector &dir, Point &pi) const
{
    dASSERT(dir.Dim() == 3, "Arg 2 wrong vector dimension");
    RVector ndir = dir/length(dir);
    double b = 2.0 * (rcam[0]*ndir[0] + rcam[1]*ndir[1] + rcam[2]*ndir[2]);
    
    double r1, r2, arg = b*b - camfac;
    if (arg < 0.0) return false;
    arg = sqrt(arg);
    r1 = (-b + arg) * 0.5;
    if (r1 < 0.0) return false;
    r2 = (-b - arg) * 0.5;
    if (r2 < 0.0) pi = cam + ndir*r1;
    else          pi = cam + ndir*r2;
    return true;
}

bool Surface_Sphere::Inside (const Point &p) const
{
    double dx = p[0]-centre[0];
    double dy = p[1]-centre[1];
    double dz = p[2]-centre[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    return (r2 <= radius*radius);
}

istream &operator>> (istream &is, Surface_Sphere &surf)
{
    is >> surf.centre >> surf.radius;
    return is;
}

ostream &operator<< (ostream &os, const Surface_Sphere &surf)
{
    os << "Surface Sphere" << endl;
    os << surf.centre << ' ' << surf.radius;
    return os;
}

// ============================================================================

Surface_Cone::Surface_Cone (): Surface3D ()
{
    tip.New(3);               // (0,0,0): default
    dir.New(3); dir[2] = 1.0; // cone opening in z-direction (default)
    sap = Pi*0.25;            // 45 deg semi-aperture (default)
    height = 1.0;             // default
    ComputeRotation(dir);
}

Surface_Cone::Surface_Cone (const Point &_tip, const RVector &_dir,
    double _sap, double _height) : Surface3D ()
{
    dASSERT(_tip.Dim() == 3, "Arg 1 invalid vector dimension");
    dASSERT(_dir.Dim() == 3, "Arg 2 invalid vector dimension");
    dASSERT(_sap > 0.0 && _sap < 0.5*Pi, "Arg 3 invalid value");
    dASSERT(_height > 0.0, "Arg 4 invalid value");

    tip.New (3); tip = _tip;
    dir.New (3); dir = _dir / length(_dir);
    sap = _sap;
    height = _height;
    ComputeRotation(dir);
}

Surface_Cone::Surface_Cone (const Surface_Cone &surf)
  : Surface3D (surf)
{
    tip.New (3); tip = surf.tip;
    dir.New (3); dir = surf.dir;
    sap = surf.sap;
    height = surf.height;
    Rot = surf.Rot;
    IRot = surf.IRot;
}

void Surface_Cone::Scale (double scale)
{
    // This scales w.r.t cone tip, not coordinate origin
    height *= scale;
}

void Surface_Cone::Scale (const RVector &scale)
{
    // need to think about this
}

void Surface_Cone::Point2Param (const Point &p, RVector &prm) const
{
    dASSERT(p.Dim() == 3, "Arg 1 dimension 3 required");
    dASSERT(prm.Dim() == 2, "Arg 2 dimension 2 required");
    RVector param(2);
    Point ps = IRot * (p-tip); // global -> local
    prm[0] = atan2 (ps[1], ps[0]); // azimuth
    prm[1] = ps[2]/height;         // elevation
}

Point Surface_Cone::Param2Point (const RVector &param) const
{
    dASSERT(param.Dim() == 2, "Arg 1 wrong vector dimension");
    Point p(3);
    p[2] = param[1]*height;
    p[0] = p[2] * cos (param[0]);
    p[1] = p[2] * sin (param[0]);
    return (Rot * p) + tip;
}

double Surface_Cone::ChordDist (const RVector &param1, const RVector &param2)
    const
{
    // Algorithm: we cut the cone from a point at the base to the tip such
    // that one of the points is on the cutting line. Then we unroll the
    // mantle onto a plane and calculate the distance between the points in
    // the plane

    // azimuth angular distance between points
    double dphi = param1[0] - param2[0];
    if (dphi < 0.0) dphi = -dphi;
    dphi = fmod (dphi, 2.0*Pi);
    if (dphi >= Pi) dphi = 2.0*Pi-dphi; // 0 <= dphi < Pi

    // map into planar pie segment
    double dphi_plane = dphi * sin(sap);
    
    // distance of points from origin in plane
    double scale = height/cos(sap);
    double r1 = param1[1]*scale;
    double r2 = param2[1]*scale;

    // distance between points in plane
    double dx = r2*cos(dphi_plane) - r1;
    double dy = r2*sin(dphi_plane);
    return hypot (dx, dy);
}

RVector Surface_Cone::Normal (const RVector &param) const
{
    RVector nml(3);
    nml[2] = -asin(sap);
    double xy = acos(sap);
    nml[0] = xy * cos(param[0]);
    nml[1] = xy * sin(param[0]);
    return Rot*nml;
}

istream &operator>> (istream &is, Surface_Cone &surf)
{
    is >> surf.tip >> surf.dir >> surf.sap >> surf.height;
    surf.ComputeRotation (surf.dir);
    return is;
}

ostream &operator<< (ostream &os, const Surface_Cone &surf)
{
    os << "Surface Cone" << endl;
    os << surf.tip << ' ' << surf.dir << ' ' << surf.sap << ' ' << surf.height;
    return os;
}

void Surface_Cone::ComputeRotation (const RVector &d)
{
    double phi = atan2 (d[1], d[0]); // azimuth angle
    double tht = acos (d[2]);        // polar angle, assuming |d|=1
    double sphi = sin(phi), cphi = cos(phi);
    double stht = sin(tht), ctht = cos(tht);

    Rot.Zero(3,3);   // glob = Rot*loc
    IRot.Zero(3,3);  // loc = IRot*glob

    Rot(0,0) = cphi*ctht;
    Rot(0,1) = -sphi;
    Rot(0,2) = cphi*stht;
    Rot(1,0) = sphi*ctht;
    Rot(1,1) = cphi;
    Rot(1,2) = sphi*stht;
    Rot(2,0) = -stht;
    Rot(2,1) = 0.0;
    Rot(2,2) = ctht;

    IRot = inverse(Rot);
}
