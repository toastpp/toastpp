// -*-C++-*-

#ifndef __GMSH_INTERFACE_H
#define __GMSH_INTERFACE_H

#include "Gmsh.h"
#include "GModel.h"
#include "felib.h"

#ifdef TOASTGMSHLIB_IMPLEMENTATION
#define TOASTGMSHLIB DLLEXPORT
#else
#define TOASTGMSHLIB DLLIMPORT
#endif

class TOASTGMSHLIB GmshInterface {
public:
    GmshInterface();
    ~GmshInterface();

    GModel *Create_Model();
    GModel *Load_Model(const char *fname);
    void LoadModel(GModel *model, const char *fname);

    void MeshModel(GModel *model, int dim,
	double minlength = -1.0, double maxlength = -1.0);
    
    void RemeshModel(GModel *model, double minlength, double maxlength = -1.0);

    void WriteMesh(GModel *model, const char *fname);
    GModel *Load_gMesh(const char *fname);
    GModel *gMesh_From_Model(GModel *model);
    QMMesh *Mesh_From_gMesh(GModel *model);
};

#endif // !__GMSH_INTERFACE_H
