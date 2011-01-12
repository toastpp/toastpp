
typedef enum{
  MUSHOMOG_G0,
  MUSINHOMOG_G0,
  MUSHOMOG_GCONST,    
  MUSINHOMOG_GENERAL,
} ScatKernType;


void genmat_angint_3D(RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Anvec, const Mesh& S2mesh);

void genmat_angint_sdm_3D(RCompRowMatrix& Aintscsc,  RCompRowMatrix& Aintscss,  RCompRowMatrix& Aintscc,  RCompRowMatrix& Aintssss,  RCompRowMatrix& Aintssc,  RCompRowMatrix& Aintcc, RCompRowMatrix& Anvec_sc, RCompRowMatrix& Anvec_ss, RCompRowMatrix& Anvec_c, const Mesh& S2mesh);

void genmat_spatint_nobf_3D(const Mesh& mesh, RCompRowMatrix& Sint, RCompRowMatrix& Sgrad, RCompRowMatrix& Sx, RCompRowMatrix& Sy, RCompRowMatrix& Sz, RCompRowMatrix* &SP);

void genmat_spatint_sdm_nobf_3D(const Mesh& mesh,  const RVector& delta, RCompRowMatrix& Sdx, RCompRowMatrix& Sdy,  RCompRowMatrix& Sdz, RCompRowMatrix& Sdxx, RCompRowMatrix& Sdxy, RCompRowMatrix& Sdyx, RCompRowMatrix& Sdyy, RCompRowMatrix& Sdxz, RCompRowMatrix& Sdzx, RCompRowMatrix& Sdyz, RCompRowMatrix& Sdzy, RCompRowMatrix& Sdzz, RCompRowMatrix* &SPdx, RCompRowMatrix* &SPdy, RCompRowMatrix* &SPdz);

void genmat_boundint_3D(RCompRowMatrix& Sint, const Mesh& mesh, const Mesh& S2mesh);

void genmat_sourceint_3D(RCompRowMatrix& Sint, const Mesh& mesh, const Mesh& S2mesh);

void precompute_3DBman_FEM(const Mesh& mesh,  const RVector& delta, const Mesh& S2mesh, const double w, const double c, CCompRowMatrix& A0, RCompRowMatrix& A1, RCompRowMatrix& A2,RCompRowMatrix& b1, RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix* &SP, RCompRowMatrix* &SPdx, RCompRowMatrix* &SPdy, RCompRowMatrix* &SPdz, RCompRowMatrix& Anvec, RCompRowMatrix& Anvec_sc, RCompRowMatrix& Anvec_ss, RCompRowMatrix& Anvec_c);

void genmat_3DBman_newF(RCompRowMatrix& A3, RCompRowMatrix& A4, const Mesh& mesh, const RDenseMatrix* sigma, const RVector& sigmatot, const int dirc, const RCompRowMatrix& Aint, const RCompRowMatrix& Aintsc, const RCompRowMatrix& Aintss, const RCompRowMatrix& Aintc, const RCompRowMatrix* SP, const RCompRowMatrix* SPdx, const RCompRowMatrix* SPdy,const RCompRowMatrix* SPdz, const RCompRowMatrix& Anvec, const RCompRowMatrix& Anvec_sc, const RCompRowMatrix& Anvec_ss, const RCompRowMatrix& Anvec_c, const ScatKernType sktyp =  MUSINHOMOG_GENERAL);

void calc_paramdistr_nobf_new_3D(RVector& sigma, RVector& sigmatot, RVector& intst, const Mesh& mesh, const RVector& muscat, const RVector& muabs, const RVector& g,const Mesh& S2mesh);

void genmat_sourceint_3D(RCompRowMatrix& Sint, const Mesh& mesh, const Mesh& S2mesh);
