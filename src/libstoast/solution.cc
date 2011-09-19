#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

using namespace std;

Solution::Solution (int nparam, int n)
{
    nprm   = nparam;
    param  = new RVector[nprm];
    active = new bool[nprm];

    for (int i = 0; i < nprm; i++) {
        param[i].New (n);
	active[i] = false;
    }
    nactive = 0;
}

Solution::Solution (const Solution &sol)
{
    nprm   = sol.nprm;
    param  = new RVector[nprm];
    active = new bool[nprm];

    for (int i = 0; i < nprm; i++) {
        param[i] = sol.param[i];
	active[i] = sol.active[i];
    }
    nactive = sol.nactive;
}

Solution::~Solution()
{
    delete []param;
    delete []active;
}

void Solution::SetParam (int which, const RVector &prm)
{ 
    param[which] = prm;
}

void Solution::Set (const RVector &prm)
{
    // this is a raw vector copy
    // assumes that val is already converted to log if required

    int i, j, idx;

#ifdef FEM_DEBUG
    int len = 0;
    for (i = 0; i < nprm; i++) if (active[i]) len += param[i].Dim();
    xASSERT (len == prm.Dim(), "Invalid parameter vector length");
#endif

    for (i = idx = 0; i < nprm; i++) {
        if (active[i])
	    for (j = 0; j < param[i].Dim(); j++)
	        param[i][j] = prm[idx++];
    }
}

void Solution::SetActive (int which, bool act)
{
    active[which] = act;
    nactive = 0;
    for (int i = 0; i < nprm; i++)
        if (active[i]) nactive++;
}

const RVector Solution::GetParam (int which) const
{
    return param[which];
}

const RVector Solution::GetActiveParams () const
{
    int i, n;
    for (i = n = 0; i < nprm; i++)
        if (active[i]) n += param[i].Dim();
    RVector p(n);
    for (i = n = 0; i < nprm; i++) {
        if (active[i]) {
	    p.Copy (param[i], n, 0, param[i].Dim());
	    n += param[i].Dim();
	}
    }
    return p;
}

void Solution::SetActiveParams (const RVector &prm)
{
    int i, n = 0;
    for (i = 0; i < nprm; i++) {
        if (active[i]) {
	    param[i].Copy (prm, 0, n, param[i].Dim());
	    n += param[i].Dim();
	}
    }
}

int Solution::ActiveDim () const
{
    int dim = 0;
    for (int i = 0; i < nprm; i++)
	if (active[i]) dim += param[i].Dim();
    return dim;
}

void Solution::Extents (int which, double &vmin, double &vmax) const
{
    vmin = vmax = param[which][0];
    for (int i = 1; i < param[which].Dim(); i++) {
	if      (param[which][i] < vmin) vmin = param[which][i];
	else if (param[which][i] > vmax) vmax = param[which][i];
    }
}

bool Solution::Bounded (int which, double vmin, double vmax)
{
    for (int i = 0; i < param[which].Dim(); i++)
        if (param[which][i] < vmin || param[which][i] > vmax)
	    return false;
    return true;
}

bool Solution::Valid ()
{
    // Note: this only works for default parameter configuration
    // (0=cmua, 1=ckappa)

    // using default parameter ranges
    if (!Bounded (0, 0, 0.05)) return false;
    if (!Bounded (1, 0.007, 1.0)) return false;
    return true;
}

void Solution::Scale (int which, double factor)
{
    param[which] *= factor;
}

void Solution::Add (const Solution &sol)
{
    for (int i = 0; i < nprm; i++)
        if (active[i]) param[i] += sol.param[i];
}
void Solution::vmul (const Solution &sol)
{
    for (int i = 0; i < nprm; i++)
        if (active[i]) param[i] *= sol.param[i];
}
void Solution::vmul (const RVector &v) const
{
    int i, j, k;
    for (i = k = 0; i < nprm; i++) {
        if (active[i])
	    for (j = 0; j < param[i].Dim(); j++)
	        param[i][j] *= v[k++];
    }
}
Solution &Solution::operator= (const Solution &sol)
{
    for (int i = 0; i < nprm; i++) {
        param[i] = sol.param[i];
	active[i] = sol.active[i];
    }
    return *this;
}

Solution Solution::operator* (double s) const
{
    Solution res (*this);
    for (int i = 0; i < nprm; i++)
	if (active[i]) res.Scale (i, s);
    return res;
}

Solution Solution::operator+ (const Solution &sol) const
{
    Solution res (*this);
    res.Add (sol);
    return res;
}

Solution Solution::operator+ (const RVector &v) const
{
    int i, j, k;
    Solution res (*this);
    for (i = k = 0; i < nprm; i++) {
        if (active[i])
	    for (j = 0; j < param[i].Dim(); j++)
	        res.param[i][j] += v[k++];
    }
    return res;
}

Solution &Solution::operator+= (const RVector &v)
{
    int i, j, k;
    for (i = k = 0; i < nprm; i++) {
        if (active[i])
	    for (j = 0; j < param[i].Dim(); j++)
	        param[i][j] += v[k++];
    }
    return *this;
}
Solution &Solution::operator*= (const RVector &v) 
{
    int i, j, k;
    for (i = k = 0; i < nprm; i++) {
        if (active[i])
	    for (j = 0; j < param[i].Dim(); j++)
	        param[i][j] *= v[k++];
    }
    return *this;
}

Solution &Solution::operator*= (double s)
{
    for (int i = 0; i < nprm; i++)
	if (active[i]) Scale (i, s);
    return *this;
}

Solution &Solution::operator+= (const Solution &sol)
{
    Add (sol);
    return *this;
}

void Solution::WriteImgGeneric (int imgno, const char *filename,
    const RVector &img, bool append)
{
    int i, n = img.Dim();
    ofstream ofs;
    if (append) ofs.open (filename, ios::app);
    else        ofs.open (filename);
    
    ofs << "Image " << imgno << endl;
    for (i = 0; i < n; i++) ofs << img[i] << ' ';
    ofs << endl;    
}

void Solution::WriteImgGeneric (int imgno, const char *filename, 
				int prmind, bool append)
{
    WriteImgGeneric (imgno, filename, param[prmind], append);
}

void Solution::WriteImg_mua (int imgno, const char *nimname, bool append)
{
    int i, n = param[0].Dim();
    ofstream ofs;
    if (append) ofs.open (nimname, ios::app);
    else        ofs.open (nimname);
    ofs << "Image " << imgno << endl;
    for (i = 0; i < n; i++)
        ofs << param[OT_CMUA][i] * (param[OT_N][i]/0.3)<< ' ';
    ofs << endl;
}

void Solution::WriteImg_mus (int imgno, const char *nimname, bool append)
{
    int i, n = param[0].Dim();
    ofstream ofs;
    if (append) ofs.open (nimname, ios::app);
    else        ofs.open (nimname);
    ofs << "Image " << imgno << endl;
    for (i = 0; i < n; i++) {
	double c = 0.3/param[OT_N][i];
        if (param[1][i]) ofs << c / (3.0*param[1][i]) - param[0][i]/c << ' ';
        else ofs << 0 << ' ';
    }
    ofs << endl;
}

double l2norm (const Solution &sol)
{
    double sum = 0.0;
    for (int i = 0; i < sol.nprm; i++)
        if (sol.active[i]) sum += l2normsq (sol.param[i]);
    return sqrt (sum);
}

Solution exp (const Solution &sol)
{
    Solution res(sol);
    for (int i = 0; i < sol.nprm; i++)
        if (res.active[i]) res.param[i] = exp(sol.param[i]);
    return res;
}
