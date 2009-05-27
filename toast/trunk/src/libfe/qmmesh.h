// ==========================================================================
// Module libfe
// File qmmesh.h
// Declaration of class QMMesh
// ==========================================================================

#ifndef __QMMESH_H
#define __QMMESH_H

// measurement profile type

typedef enum {
    PROFILE_POINT,
    PROFILE_GAUSSIAN,
    PROFILE_COSINE,
    PROFILE_TOPHAT
} MProfileType;

// ==========================================================================
// class QMMesh

class FELIB QMMesh : public Mesh {
public:
    int nQ, nM;		// number of source/measurement points
    int nQM;		// total number of measurements (not always nQ*nM !)
    Point *Q, *M;	// lists of source/measurement points

    int *nQMref;	// nQMref[q]: number of detectors connected to q
    int **QMref;	// QMref[q][i] is the absolute detector index
                        // for the i-th detector connected to q
    int *Qofs;          // offset of source blocks in data array
    int **QMofs;        // QMofs[q][m] is the absolute offset into data
                        // array for measurement from source q and detector
                        // m (or -1 if combination is not used)

    int *Mel;		// list of measurement elements
    RDenseMatrix *Mcosingder; // list of cosin*gder matrices for measurement
                              // elements
    RVector *source_profile;
    RVector *meas_profile;
    // Arrays of dimension nQ and nM of boundary source and measurement
    // profile vectors, respectively.
    // Each vector is of dimension nbnd() (i.e. number of boundary nodes)
    // and contains the weights of measurement i at each boundary node
    // Valid after InitM

    bool fixed_q_pos;   // TRUE if Q contains final source positions
    bool fixed_m_pos;   // TRUE if M contains final measurement positions

    QMMesh ();
    ~QMMesh ();

    void SetupQM (Point *q, int nq, Point *m, int nm);
    void SetupRegularQM (int nqm);
    void LoadQM (std::istream &is);

    void ScaleMesh (double scale);
    // rescale the whole mesh

    void ScaleMesh (const RVector &scale);
    // Anisotropic scaling

    bool Connected (int q, int m) const;
    bool Connected (int qm) const;
    // returns TRUE if measurement with source q and detector m exists

    int Meas (int q, int m) const;
    // returns the index number of absolute measurement 'm' for source 'q'
    // (i.e. inverse of QMref: QMref[q][Meas(q,m)] = m)
    // returns -1 if q and m not connected

    void GetMesh (std::istream& i);	// read mesh definition from stream
    void PutMesh (std::ostream& o);	// write mesh definition to stream
    // IO functions

    MProfileType *mptype;  // list of measurement profiles
    double *mwidth;        // measurement width parameters
    double *msup;          // measurement support parameters

    MProfileType *qptype;  // list of source profiles
    double *qwidth;        // source width parameters
    double *qsup;          // source support parameters

private:
    void InitM ();	// precalculation of measurement matrices
};

#endif // !__QMMESH_H
