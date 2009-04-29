#ifndef __TOASTTYPE_H
#define __TOASTTYPE_H

// =========================================================================
// type definitions for parameter sets

typedef enum {
    RESET_MESH, RESET_HOMOG, RESET_NIM
} PRM_RESET_MODE;

typedef enum {
    PROF_GAUSSIAN, PROF_COSINE, PROF_COMPLETETRIG, PROF_POINT
} SRC_PROFILE;

typedef enum {
    SRCMODE_NEUMANN,
    SRCMODE_ISOTROPIC
} SourceMode;

typedef enum {          // measurement types
    MEAS_INTENSITY,         //   CW only
    MEAS_FRE_FIM,           //   complex intensity (re+im)
    MEAS_FMOD_FARG          //   complex intensity (mod+arg)
} Measurement;

typedef struct {
    char prmname[256];
    char logname[256];
    char meshname[256];
    char qmname[256];
    Measurement dtype;
    double freq;
    SourceMode qtype;
    SRC_PROFILE qprof;
    double qwidth;
    SRC_PROFILE mprof;
    double mwidth;
    struct {
	PRM_RESET_MODE resettp;
	union {
	    double homog;
	    char nimf[256];
	};
    } initprm[3];
} PARAMS;

#endif // !__TOASTTYPE_H
