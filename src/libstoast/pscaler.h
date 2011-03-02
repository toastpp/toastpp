// -*-C++-*-
// ==========================================================================
// Parameter scaling
// ==========================================================================

#ifndef __PSCALER_H
#define __PSCALER_H

#include "stoastlib.h"

// =========================================================================

class STOASTLIB Scaler {
public:
    Scaler () {}
    virtual const char *ScalingMethod() const = 0;
    virtual RVector Scale (const RVector &unscaled) const = 0;
    virtual RVector Unscale (const RVector &scaled) const = 0;
    virtual RVector JScale (const RVector &unscaled) const = 0;
    virtual void ScaleJacobian (const RVector &unscaled, RMatrix &J) const {
	return J.ColScale (JScale (unscaled));
    }
    virtual void ScaleGradient (const RVector &unscaled, RVector &grad) const {
	grad *= JScale (unscaled);
    }
};

class STOASTLIB NullScaler: public Scaler {
public:
    NullScaler (): Scaler () {}
    const char *ScalingMethod () const { return "NULL"; }
    RVector Scale (const RVector &unscaled) const { return unscaled; }
    RVector Unscale (const RVector &scaled) const { return scaled; }
    RVector JScale (const RVector &unscaled) const {
	RVector s(unscaled.Dim()); s = 1.0; return s;
    }
    void ScaleJacobian (const RVector &unscaled, RMatrix &J) const {}
    void ScaleGradient (const RVector &unscaled, RVector &grad) const {}
};

class STOASTLIB ConstScaler: public Scaler {
public:
    ConstScaler (const RVector &scale): Scaler () { sc = scale; }
    const char *ScalingMethod () const { return "CONST"; }
    RVector Scale (const RVector &unscaled) const { return unscaled*sc; }
    RVector Unscale (const RVector &scaled) const { return scaled/sc; }
    RVector JScale (const RVector &unscaled) const { return inv(sc); }
private:
    RVector sc;
};

class STOASTLIB LogScaler: public Scaler {
public:
    LogScaler (): Scaler () {}
    const char *ScalingMethod () const { return "LOG"; }
    RVector Scale (const RVector &unscaled) const { return log(unscaled); }
    RVector Unscale (const RVector &scaled) const { return exp(scaled); }
    RVector JScale (const RVector &unscaled) const { return unscaled; }
};

class STOASTLIB LinLogScaler: public Scaler {
public:
    LinLogScaler (const RVector &scale): Scaler () { sc = scale; }
    const char *ScalingMethod () const { return "LINLOG"; }
    RVector Scale (const RVector &unscaled) const { return log(unscaled*sc); }
    RVector Unscale (const RVector &scaled) const { return exp(scaled)/sc; }
    RVector JScale (const RVector &unscaled) const { return unscaled; }
private:
    RVector sc;
};

class STOASTLIB BoundLogScaler: public Scaler {
public:
    BoundLogScaler (const RVector &_xmin, const RVector &_xmax): Scaler()
    { xmin = _xmin, xmax = _xmax; }
    const char *ScalingMethod () const { return "BOUNDED_LOG"; }
    RVector Scale (const RVector &unscaled) const
    { return log ((unscaled-xmin)*(xmax-unscaled)); }
    RVector Unscale (const RVector &scaled) const;
    RVector JScale (const RVector &unscaled) const;
private:
    RVector xmin, xmax;
};

#endif // !__PSCALER_H
