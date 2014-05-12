// ============================================================================
// TOAST module                                       (c) Martin Schweiger 2001
// lem: linear elasticity module (Tcl/Tk application)
//
// Requires:
//    * Tcl/Tk 8.0
//    * Tix 4.1
//    * BLT
//    * Tcl script $TOAST_SCRIPT_PATH/lem.tcl
// ============================================================================

#include <togl.h>
#include <GL/glut.h>
#include <tcl.h>
#include <tk.h>
#include <tix.h>
#include <iostream.h>
#include "toasttcl.h"
#include "mathlib.h"
#include "felib.h"
#include "gmres.h"

// ============================================================================
// Prototypes

void create_cb   (struct Togl *togl);
void display_cb  (struct Togl *togl);
void reshape_cb  (struct Togl *togl);
int  setXrot_cb  (struct Togl *togl, int argc, char *argv[]);
int  setYrot_cb  (struct Togl *togl, int argc, char *argv[]);
int  setZrot_cb  (struct Togl *togl, int argc, char *argv[]);
int  setDist_cb  (struct Togl *togl, int argc, char *argv[]);
int  showDisp_cb (struct Togl *togl, int argc, char *argv[]);
int  Solve_cb    (struct Togl *togl, int argc, char *argv[]);
int  getXrot_cb  (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[]);
int  getYrot_cb  (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[]);
int  getDist_cb  (ClientData clientData, Tcl_Interp *interp,
		  int argc, char *argv[]);
void display_mesh (void);
void init_display (void);
void InitRGB ();
void UpdateMesh (RVector &sol);
void AllocVtxColour (RVector &disp,
    IVector &dspx, IVector &dspy, IVector &dspz);

// ============================================================================
// Global variables

typedef struct {
    unsigned char r, g, b;
} RGB;

const char *PalettePath = "/home/ucacmas/toast/scales/";
const char *Palette = "fire2.pal";

char *tclscriptpath;
Tcl_Interp *g_interp;
static GLuint FontBase;
static float xAngle = 0.0, yAngle = 0.0, zAngle = 0.0;
static GLfloat CornerX, CornerY, CornerZ;  /* where to print strings */

double MAX_SIZE =20000;
double MIN_SIZE =-20000;
double xmin = -1000.0;
double ymin = -1000.0;
double zmin = -1000.0;
double xmax =  1000.0;
double ymax =  1000.0;
double zmax =  1000.0;
double value=0.0;
double range=0.5;
double RotX = 0.0;
double RotY = 0.0;
double Xpoint = 100.0;
int Xres =600;
int Yres = 600;

bool boundary = true;
bool voidsurf = false;
bool wireframe = false;
//int   nvtx;
//float *vtxptr;
//float *nmlptr;
//float *col1, *col2, *col3, *colptr;
IVector xidx, yidx, zidx;
RGB rgb[256];

static void print_string (const char *s)
{ glCallLists (strlen (s), GL_UNSIGNED_BYTE, s); }

class RunParam {
public:
    RunParam();
    ~RunParam();
    bool LoadMesh (char *fname);
    bool LoadMaterialFile (char *fname);
    bool LoadDispFile (char *fname);
    void Run (struct Togl *togl);
    void AllocSurfBufs ();
    Mesh mesh;
    RVector x;
    int *dispnodes;
  
    float *vtx, *nml;
    float *col1, *col2, *col3, *col;
    int ntri;

private:
    char meshfile[256];
    char dispfile[256];
    char bforcefile[256];
    RCompRowMatrix F;
    RVector rhs;
    RVector nu, E;
    double tol;
} runprm;

RunParam::RunParam()
{
    meshfile[0] = '\0';
    dispfile[0] = '\0';
    bforcefile[0] = '\0';
    dispnodes = 0;
    ntri = 0;
    tol = 1e-13;
}

RunParam::~RunParam()
{
    if (dispnodes) delete []dispnodes;
    if (ntri) {
        delete []vtx;
	delete []nml;
	delete []col1;
	delete []col2;
	delete []col3;
    }
}

bool RunParam::LoadMesh (char *fname)
{
    cerr << "Loading mesh ... " << flush;
    if (strcmp (fname, meshfile)) {
        strcpy (meshfile, fname);
	ifstream ifs (meshfile);
	if (!ifs) return false;
	ifs >> mesh;
	if (!ifs.good()) return false;
	mesh.Setup();
    }
    cerr << "finished!" << endl;

    // Allocate system stiffness matrix
    cerr << "Allocating system matrix ... " << flush;
    int n   = mesh.nlen();
    int dim = mesh.Dimension();
    int *rowptr, *colidx, *drowptr, *dcolidx, nzero, dn;
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    BlockExpand (rowptr, colidx, n, drowptr, dcolidx, dn, dim, dim);
    F.New (dn, dn);
    F.Initialise (drowptr, dcolidx);
    delete []rowptr;
    delete []colidx;
    delete []drowptr;
    delete []dcolidx;

    // Allocate rhs and solution
    rhs.New (dn);
    x.New (dn);
    cerr << "Finished!" << endl;

    // Allocate parameter arrays
    E.New (mesh.elen());
    nu.New (mesh.elen());

    if (dispnodes) delete []dispnodes;
    dispnodes = new int[n];

    // Allocate render-related
    // data buffers
    AllocSurfBufs ();

    return true;
}

void RunParam::AllocSurfBufs ()
{
    int i, j, k, sn, nsn, nbndsd, *node;
    Element *pel;
    Mesh &mesh = runprm.mesh;

    // find the number of boundary sides
    for (i = nbndsd = 0; i < mesh.elen(); i++) {
        pel = mesh.elist[i];
	node = pel->Node;
	for (j = 0; j < pel->nSide(); j++) {
	    nsn = pel->nSideNode(j);
	    for (k = 0; k < nsn; k++) {
	        sn = node[pel->SideNode(j,k)];
		if (!mesh.nlist[sn].isBnd()) break;
	    }
	    if (k == nsn) nbndsd++;
	}
    }
    cout << "Found " << nbndsd << " boundary sides" << endl;

    if (ntri) {
        delete []vtx;
	delete []nml;
	delete []col1;
	delete []col2;
	delete []col3;
    }
    // create vertex, normal and index buffers
    // 3 vertices per triangle, 3 coordinates per vertex
    vtx = new float[nbndsd*3*3];
    nml = new float[nbndsd*3*3];
    //int   *vtxidx = new int[nbndsd*3];
    // create buffers for displacement vertex colours
    col1 = new float[nbndsd*3*3];
    col2 = new float[nbndsd*3*3];
    col3 = new float[nbndsd*3*3];
    col  = col1;

    ntri = nbndsd;
}

bool RunParam::LoadMaterialFile (char *fname)
{
    int i, nmat, nidx, idx;

    cerr << "Loading material file ... " << flush;
    ifstream ifs (fname);
    if (!ifs) return false;
    ifs >> nmat;
    struct MatList { double E, nu, dns; } *matlist = new MatList[nmat];
    for (i = 0; i < nmat; i++)
        ifs >> matlist[i].E >> matlist[i].nu >> matlist[i].dns;
    ifs >> nidx;
    if (nidx != mesh.elen()) return false;
    for (i = 0; i < nidx; i++) {
        ifs >> idx;
	if (idx < 0 || idx >= nmat) return false;
	nu[i] = matlist[idx].nu;
	E[i]  = matlist[idx].E;
    }
    cerr << "finished!" << endl;

    // assemble system matrix
    cerr << "Assembling system matrix ... " << flush;
    AddToSysMatrix_elasticity (mesh, F, E, nu);
    cerr << "finished!" << endl;

    return true;
}

bool RunParam::LoadDispFile (char *fname)
{
    const double BIGSPRING = 1e7;
    const char *WS = " \t";
    char cbuf[256], *s;
    int i, nd, res;
    double disp;

    if (!strcmp (fname, "[none]")) return true; // no displacements

    cerr << "Assigning explicit displacements ... " << flush;
    ifstream ifs (fname);
    if (!ifs) return false;
    while (ifs.getline (cbuf, 256)) {
        s = strtok (cbuf, WS);
	res = sscanf (s, "%d", &nd);
	if (!res) continue;
	for (i = 0; i < 3; i++) {
	    s = strtok (NULL, WS);
	    res = sscanf (s, "%lf", &disp);
	    if (res) {
	        F(nd*3+i,nd*3+i) *= BIGSPRING;
		rhs[nd*3+i] = F(nd*3+i,nd*3+i)*disp;
	    }
	}  
    }
    cerr << "finished!" << endl;

    return true;
}

struct Togl *togl_glob;

void RenderIntermediate (void *data)
{
    Tcl_Interp *interp = Togl_Interp(togl_glob);
    RVector *x = (RVector*)data;
    UpdateMesh (*x);
    //Togl_PostRedisplay (togl_glob);
    display_cb(togl_glob);
    Tcl_Eval (interp, "update");
}

void RunParam::Run (struct Togl *togl)
{
    togl_glob = togl;  // temporary

    cerr << "Starting GMRES solver" << endl;
    RPreconditioner *precon = new RPrecon_Diag;
    precon->Reset (&F);
    int niter;
    //niter = BiCGSTAB (F, rhs, x, tol, precon);
    //niter = GMRES (F, rhs, x, tol, precon, &RenderIntermediate);
    niter = gmres (10, F, rhs, x, precon, tol, &RenderIntermediate);
    delete precon;
    UpdateMesh (x);
    cerr << "converged after " << niter << " iterations" << endl;
    Togl_PostRedisplay (togl);
}

// ============================================================================
// Prototypes

int Tk_AppInit (Tcl_Interp *interp);
void ErrorHandler (char *msg);

// ============================================================================
// main

int main (int argc, char *argv[])
{
    CHECK_EXPIRED();

    SetErrorhandler (ErrorHandler);

    if (!(tclscriptpath = getenv ("TOAST_SCRIPT_PATH"))) {
        cerr << "Environment variable TOAST_SCRIPT_PATH not defined.\n";
	cerr << "Exiting.\n";
	exit (1);
    }
    Tk_Main (argc, argv, Tk_AppInit);
    exit (0);
}

int Tk_AppInit (Tcl_Interp *interp)
{
    char cbuf[256];

    if (Tcl_Init  (interp) == TCL_ERROR) return TCL_ERROR;
    if (Tk_Init   (interp) == TCL_ERROR) return TCL_ERROR;
    if (Tix_Init  (interp) == TCL_ERROR) return TCL_ERROR;
    if (Togl_Init (interp) == TCL_ERROR) return TCL_ERROR;

    Togl_CreateFunc  (create_cb);
    Togl_DisplayFunc (display_cb);
    Togl_ReshapeFunc (reshape_cb);
    Togl_CreateCommand ("setXrot", setXrot_cb);
    Togl_CreateCommand ("setYrot", setYrot_cb);
    Togl_CreateCommand ("setZrot", setZrot_cb);
    Togl_CreateCommand ("setDist", setDist_cb);
    Togl_CreateCommand ("showDisp", showDisp_cb);
    Togl_CreateCommand ("Solve",   Solve_cb);

    Tcl_CreateCommand (interp, "getXrot", getXrot_cb, (ClientData)NULL,
		       (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand (interp, "getYrot", getYrot_cb, (ClientData)NULL,
		       (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand (interp, "getDist", getDist_cb, (ClientData)NULL,
		       (Tcl_CmdDeleteProc *)NULL);
    
    g_interp = interp;

    InitRGB();

    strcpy (cbuf, tclscriptpath);
    strcat (cbuf, "/lem.tcl");
    if (Tcl_EvalFile (interp, cbuf) == TCL_ERROR)
        return TCL_ERROR;

    return TCL_OK;
}

void ErrorHandler (char *msg)
{
    cerr << endl << msg << endl;
    Tcl_Eval (g_interp, "ErrorBox");
}

int Solve_cb (struct Togl *togl, int argc, char *argv[])
{
    cerr << "argc=" << argc << endl;
    for (int i = 0; i < argc; i++)
      cerr << argv[i] << endl;

    // Run the linear elasticity solver with provided parameters

    Tcl_Interp *interp = Togl_Interp(togl);
    cerr << "label2" << endl;
    char *varname = argv[2], *bforcefile;
    if (!runprm.LoadMesh (Tcl_GetVar2 (interp, varname, "meshfile", 0)))
        ErrorHandler ("Error loading mesh");
    cerr << "label3" << endl;

    if (!runprm.LoadMaterialFile (Tcl_GetVar2 (interp, varname, "matfile", 0)))
        ErrorHandler ("Error loading material file");
    cerr << "label4" << endl;

    if (!runprm.LoadDispFile (Tcl_GetVar2 (interp, varname, "dispfile", 0)))
        ErrorHandler ("Error loading displacement file");
    cerr << "label5" << endl;

    bforcefile = Tcl_GetVar2 (interp, varname, "bforcefile", 0);

    runprm.Run (togl);

    return TCL_OK;
}
    
void create_cb (struct Togl *togl)
{
    FontBase = Togl_LoadBitmapFont (togl, TOGL_BITMAP_8_BY_13);
    if (!FontBase) {
        printf("Couldn't load font!\n");
	exit(1);
    }
    init_display();
}

/*
 * Togl widget display callback.  This is called by Tcl/Tk when the widget's
 * contents have to be redrawn.  Typically, we clear the color and depth
 * buffers, render our objects, then swap the front/back color buffers.
 */
void display_cb (struct Togl *togl)
{ 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.0,0.0,0.0,0.0);

    glLoadIdentity();	/* Reset modelview matrix to the identity matrix */
    glTranslatef(0.0, 0.0, -Xpoint);   /* Move the camera back three units */
    glRotatef(xAngle, 1.0, 0.0, 0.0);  /* Rotate by X, Y, and Z angles */
    glRotatef(yAngle, 0.0, 1.0, 0.0);
    glRotatef(zAngle, 0.0, 0.0, 1.0);

    if (runprm.mesh.nlen() > 0) // have mesh?
	display_mesh();
    Togl_SwapBuffers (togl);
}

/*
 * Togl widget reshape callback.  This is called by Tcl/Tk when the widget
 * has been resized.  Typically, we call glViewport and perhaps setup the
 * projection matrix.
 */
void reshape_cb (struct Togl *togl)
{
   int width = Togl_Width( togl );
   int height = Togl_Height( togl );
   float aspect = (float) width / (float) height;

   glViewport( 0, 0, width, height );

   /* Set up projection transform */
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-aspect, aspect, -1.0, 1.0, 1.0, 1000.0);

   CornerX = -aspect;
   CornerY = -1.0;
   CornerZ = -1.1;

   /* Change back to model view transform for rendering */
   glMatrixMode(GL_MODELVIEW);
}

int setXrot_cb (struct Togl *togl, int argc, char *argv[])
{
   Tcl_Interp *interp = Togl_Interp(togl);

   /* error checking */
   if (argc != 3) {
      Tcl_SetResult( interp,
                     "wrong # args: should be \"pathName setXrot ?angle?\"",
                     TCL_STATIC );
      return TCL_ERROR;
   }

   xAngle = atof( argv[2] );
   
/* printf( "before %f ", xAngle ); */

   if ( xAngle < 0.0 ) {
     xAngle += 360.0;
   } else if ( xAngle > 360.0 ) {
     xAngle -= 360.0;
   }

/* printf( "after %f \n", xAngle ); */

   Togl_PostRedisplay (togl);

   /* Let result string equal value */
   strcpy( interp->result, argv[2] );
   return TCL_OK;
}

int setYrot_cb (struct Togl *togl, int argc, char *argv[])
{
   Tcl_Interp *interp = Togl_Interp(togl);

   /* error checking */
   if (argc != 3) {
      Tcl_SetResult( interp,
                     "wrong # args: should be \"pathName setYrot ?angle?\"",
                     TCL_STATIC );
      return TCL_ERROR;
   }

   yAngle = atof( argv[2] );
   
   if ( yAngle < 0.0 ) {
     yAngle += 360.0;
   } else if ( yAngle > 360.0 ) {
     yAngle -= 360.0;
   }

   Togl_PostRedisplay(togl);

   /* Let result string equal value */
   strcpy( interp->result, argv[2] );
   return TCL_OK;
}

int setZrot_cb (struct Togl *togl, int argc, char *argv[])
{
   Tcl_Interp *interp = Togl_Interp(togl);

   /* error checking */
   if (argc != 3) {
      Tcl_SetResult( interp,
                     "wrong # args: should be \"pathName setZrot ?angle?\"",
                     TCL_STATIC );
      return TCL_ERROR;
   }

   zAngle = atof( argv[2] );
   
   if ( zAngle < 0.0 ) {
     zAngle += 360.0;
   } else if ( zAngle > 360.0 ) {
     zAngle -= 360.0;
   }

   Togl_PostRedisplay(togl);

   /* Let result string equal value */
   strcpy( interp->result, argv[2] );
   return TCL_OK;
}

int setDist_cb (struct Togl *togl, int argc, char *argv[])
{
    Xpoint = atof (argv[2]);
    Togl_PostRedisplay (togl);
    return TCL_OK;
}

int showDisp_cb (struct Togl *togl, int argc, char *argv[])
{
    char axis = argv[2][0];
    switch (axis) {
    case 'x':
        runprm.col = runprm.col1;
	break;
    case 'y':
        runprm.col = runprm.col2;
	break;
    case 'z':
        runprm.col = runprm.col3;
	break;
    }
    Togl_PostRedisplay (togl);
    return TCL_OK;
}

int getXrot_cb (ClientData clientData, Tcl_Interp *interp,
		int argc, char *argv[])
{
    sprintf (interp->result, "%d", (int)xAngle);
    return TCL_OK;
}

int getYrot_cb (ClientData clientData, Tcl_Interp *interp,
		int argc, char *argv[])
{
    sprintf (interp->result, "%d", (int)yAngle);
    return TCL_OK;
}

int getDist_cb (ClientData clientData, Tcl_Interp *interp,
		int argc, char *argv[])
{
    sprintf (interp->result, "%d", (int)Xpoint);
    return TCL_OK;
}

void set_nodes (void)
{
    Mesh &mesh = runprm.mesh;
    for (int i = 0; i < mesh.nlen(); i++) {
        runprm.dispnodes[i] = 1;
	if (boundary && voidsurf && !mesh.nlist[i].isBnd())
	    runprm.dispnodes[i] = 0;
#ifdef TOAST_VOID_VERSION
	if (voidsurf && !mesh.nlist[i].isVBnd())
	    runprm.dispnodes[i] = 0;
#endif
	if (boundary && !mesh.nlist[i].isBnd())
	    runprm.dispnodes[i]=0;
#ifdef TOAST_VOID_VERSION
	if (voidsurf && boundary  &&mesh.nlist[i].isBnd())
	    runprm.dispnodes[i]=1;
#endif
    }
}

void surf_node_ptr ()
{
    int i, j, k, l, idx, nsn, sn, nd, nbndnd, ntri, *node;
    Element *pel;
    Mesh &mesh = runprm.mesh;

    int    nbndsd  = runprm.ntri;
    float *vtxptr  = runprm.vtx;
    float *nmlptr  = runprm.nml;
    float *colptr1 = runprm.col1;
    float *colptr2 = runprm.col2;
    float *colptr3 = runprm.col3;
    int   *vtxidx  = new int[nbndsd*3];

    for (i = ntri = 0; i < mesh.elen(); i++) {
        pel = mesh.elist[i];
	node = pel->Node;

	for (j = 0; j < pel->nSide(); j++) {
	    nsn = pel->nSideNode(j);
	    Point lcnt = pel->SideCentre(j);
	    Point gcnt = pel->Global(mesh.nlist, lcnt);
	    RDenseMatrix lder = pel->LocalShapeD (lcnt);
	    RDenseMatrix jacin = inverse (lder *mesh.ElGeom(i));
	    for (k = 0; k < nsn; k++) {
	        sn = node[pel->SideNode(j,k)];
		if (!mesh.nlist[sn].isBnd()) break;
	    }
	    if (k == nsn) {
	        // generate vertices
	        for (k = 0; k < 3; k++) {
  		    // note for higher order shape functions we skip over
		    // midpoints here
		    sn = node[pel->SideNode(j,k)];
		    for (l = 0; l < 3; l++)
		        vtxptr[ntri*9+k*3+l] = mesh.nlist[sn][l];
		    vtxidx[ntri*3+k] = sn;
		}
		// generate side normal
		RVector dcos = pel->DirectionCosine (j, jacin);
		for (k = 0; k < 3; k++) {
		    for (l = 0; l < 3; l++)
		        nmlptr[ntri*9+k*3+l] = dcos[l];
		}
		// generate colour
		for (k = 0; k < 3; k++) {
		    sn = node[pel->SideNode(j,k)];
		    idx = xidx[sn];
		    colptr1[ntri*9+k*3+0] = (float)rgb[idx].r / 256.0;
		    colptr1[ntri*9+k*3+1] = (float)rgb[idx].g / 256.0;
		    colptr1[ntri*9+k*3+2] = (float)rgb[idx].b / 256.0;
		    idx = yidx[sn];
		    colptr2[ntri*9+k*3+0] = (float)rgb[idx].r / 256.0;
		    colptr2[ntri*9+k*3+1] = (float)rgb[idx].g / 256.0;
		    colptr2[ntri*9+k*3+2] = (float)rgb[idx].b / 256.0;
		    idx = zidx[sn];
		    colptr3[ntri*9+k*3+0] = (float)rgb[idx].r / 256.0;
		    colptr3[ntri*9+k*3+1] = (float)rgb[idx].g / 256.0;
		    colptr3[ntri*9+k*3+2] = (float)rgb[idx].b / 256.0;
		}
		ntri++;
	    }
	}
    }
    if (ntri != nbndsd) cout << "Problems!" << endl;
    
#define SMOOTH_SHADING
#ifdef SMOOTH_SHADING
    // for smooth surfaces we need to average normals
    float *avgnml = new float[mesh.nlen()*3];
    int *nnml = new int[mesh.nlen()];
    for (i = 0; i < mesh.nlen()*3; i++) avgnml[i] = 0;
    for (i = 0; i < mesh.nlen(); i++) nnml[i] = 0;
    for (i = 0; i < nbndsd*3; i++) {
        nd = vtxidx[i];
	for (j = 0; j < 3; j++) avgnml[nd*3+j] += nmlptr[i*3+j];
	nnml[nd]++;
    }
    for (i = 0; i < mesh.nlen(); i++) {
        double x = avgnml[i*3];
	double y = avgnml[i*3+1];
	double z = avgnml[i*3+2];
	double len = sqrt (x*x + y*y + z*z);
	for (j = 0; j < 3; j++) avgnml[i*3+j] /= len;
    }
    for (i = 0; i < nbndsd*3; i++) {
        nd = vtxidx[i];
	for (j = 0; j < 3; j++) nmlptr[i*3+j] = avgnml[nd*3+j];
    }
    delete []avgnml;
    delete []nnml;
#endif

    delete []vtxidx;
}

void AllocVtxColour (RVector &disp,
    IVector &dspx, IVector &dspy, IVector &dspz)
{
    Mesh    &mesh = runprm.mesh;

    int i, val, n = mesh.plist.Len();
    double dx, dy, dz, irange;
    double dispmin = 1e100, dispmax = -1e100;

    // find parameter range
    for (i = 0; i < n; i++) {
        dx = disp[i*3+0];
	if (dx < dispmin) dispmin = dx;
	if (dx > dispmax) dispmax = dx;
	dy = disp[i*3+1];
	if (dy < dispmin) dispmin = dy;
	if (dy > dispmax) dispmax = dy;
	dz = disp[i*3+2];
	if (dz < dispmin) dispmin = dz;
	if (dz > dispmax) dispmax = dz;
    }
    irange = 256.0/(dispmax-dispmin);

    dspx.New(n);
    dspy.New(n);
    dspz.New(n);

    for (i = 0; i < n; i++) {
        dspx[i] = (int)((disp[i*3+0]-dispmin)*irange - 0.5);
	if      (dspx[i] < 0)    dspx[i] = 0;
	else if (dspx[i] >= 256) dspx[i] = 255;
        dspy[i] = (int)((disp[i*3+1]-dispmin)*irange - 0.5);
	if      (dspy[i] < 0)    dspy[i] = 0;
	else if (dspy[i] >= 256) dspy[i] = 255;
        dspz[i] = (int)((disp[i*3+2]-dispmin)*irange - 0.5);
	if      (dspz[i] < 0)    dspz[i] = 0;
	else if (dspz[i] >= 256) dspz[i] = 255;
    }
}

void InitRGB ()
{
    int i, r, g, b;
    char cbuf[256];
    strcpy (cbuf, PalettePath);
    strcat (cbuf, Palette);
    ifstream ifs (cbuf);
    for (i = 0; i < 256; i++) {
        ifs >> r >> g >> b;
	rgb[i].r = (unsigned char)r;
	rgb[i].g = (unsigned char)g;
	rgb[i].b = (unsigned char)b;
    }
}

void display_mesh (void)
{
    int i, j, k, nsn, nd;
    Element *pel;
    Mesh &mesh = runprm.mesh;

    if (!wireframe) {
	glEnableClientState (GL_VERTEX_ARRAY);
	glEnableClientState (GL_NORMAL_ARRAY);
	glEnableClientState (GL_COLOR_ARRAY);
        glVertexPointer (3, GL_FLOAT, 3*sizeof(float), runprm.vtx);
	glNormalPointer (   GL_FLOAT, 3*sizeof(float), runprm.nml);
	glColorPointer  (3, GL_FLOAT, 3*sizeof(float), runprm.col);
	glDrawArrays    (GL_TRIANGLES, 0, runprm.ntri*3);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	return;
    }

    for (i = 0; i < mesh.elen(); i++) {
       pel = mesh.elist[i];
       if(pel->Size() > MAX_SIZE) continue;
       if(pel->Size() < MIN_SIZE) continue;
       for (j = 0; j < pel->nSide(); j++) {
	   nsn = pel->nSideNode(j);
	   for (k = 0; k < nsn; k++) {
	       nd = pel->Node[pel->SideNode(j,k)];
	       if (runprm.dispnodes[nd] == 0)       break;
	       if (mesh.nlist[nd][0] < xmin) break;
	       if (mesh.nlist[nd][1] < ymin) break;
	       if (mesh.nlist[nd][2] < zmin) break;
	       if (mesh.nlist[nd][0] > xmax) break;
	       if (mesh.nlist[nd][1] > ymax) break;
	       if (mesh.nlist[nd][2] > zmax) break;
	   }
	   if(k == nsn) {
	       glBegin(GL_LINE_LOOP);
	       //glColor3f(0.0,0.0,1.0);
	       for (k = 0; k < nsn; k++) {
		   nd = pel->Node[pel->SideNode(j,k)];
		   glVertex3f(mesh.nlist[nd][0], mesh.nlist[nd][1], mesh.nlist[nd][2]);
	       } 
	       glEnd();
	       if (!wireframe) {
		   glBegin(GL_POLYGON);
		   //glColor3f(1.0,0.0,0.0);
		   for(k = 0; k < nsn; k++) {
		       nd = pel->Node[pel->SideNode(j,k)];
		       glVertex3f(mesh.nlist[nd][0], mesh.nlist[nd][1], mesh.nlist[nd][2]);
		       glNormal3f(mesh.nlist[nd][0], mesh.nlist[nd][1], mesh.nlist[nd][2]);
		   }
		   glEnd();
	       }
	   }  
       }
    }
    glFlush();  /* Single buffered, so needs a flush. */
}

void init_display (void) {
    //tbInit(GLUT_MIDDLE_BUTTON);
    //tbAnimate(GL_TRUE);
    /* Use depth buffering for hidden surface elimination. */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);

    /* Setup the view of the cube. */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( /* field of view in degree */ 90,
		    /* aspect ratio */ 1.0,
		    /* z near */ 1.0, /* Z far */ 10000.0);
    glViewport(0,0,Xres,Yres);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(Xpoint, 0.0, 0.0,  /* eye is at (0,0,5) */
	      0.0, 0.0, 0.0,      /* center is at (0,0,0) */
	      0.0, 1.0, 0.);      /* up is in positive Y direction */

    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };

    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    //glDepthFunc(GL_LESS);
    //glEnable(GL_DEPTH_TEST);
}

void UpdateMesh (RVector &sol)
{
   AllocVtxColour (sol, xidx, yidx, zidx);
   surf_node_ptr();
}
