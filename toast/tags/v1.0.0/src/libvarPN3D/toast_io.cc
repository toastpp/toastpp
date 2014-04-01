#include "toast_io.h"
#include "pparse.h"
#include "source.h"
#include "detector.h"

using namespace std;
void ReadMesh(const char* fname, QMMesh &qmmesh)
{
    cout << "MESH: " << fname << endl;
    ifstream ifs;
    ifs.open (fname);
    xASSERT (ifs.is_open(), "Mesh file not found.");
    ifs >> qmmesh;
    xASSERT (ifs.good(), "Problem reading mesh.");
    ifs.close ();
    cout << "* " << qmmesh.elen() << " elements, " << qmmesh.nlen()
	 << " nodes\n";
    int dimension = qmmesh.nlist[0].Dim();
    for (int i = 1; i < qmmesh.nlist.Len(); i++)
	xASSERT(qmmesh.nlist[i].Dim() == dimension, "Inconsistent node dimensions.");
    qmmesh.Setup();

}

void ReadQM(const char* fname, QMMesh &qmmesh, int &nQ, int &nM,  RCompRowMatrix &qvec, RCompRowMatrix &mvec)
{
    ParamParser pp;
    int qprof, mprof;   // source/measurement profile (0=Gaussian, 1=Cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]
    SourceMode srctp;

    cout << "QM: " <<  fname<< endl;
    ifstream ifs;
    ifs.open (fname);
    xASSERT (ifs.is_open(), "QM file not found.");
    qmmesh.LoadQM (ifs);
    xASSERT (ifs.good(), "Problem reading QM.");
    ifs.close ();
    nQ = qmmesh.nQ;
    nM = qmmesh.nM;
    cout<<"Number of sources: "<<nQ<<endl;
    cout<<"Number of detectors: "<<nM<<endl;
    SelectSourceProfile (pp, qprof, qwidth, srctp);
    qvec.New (nQ, qmmesh.nlen());
    for (int i = 0; i < nQ; i++) {
	RVector q(qmmesh.nlen());
	switch (qprof) {
	case 0:
	    q = QVec_Point (qmmesh, qmmesh.Q[i], srctp);
	    break;
	case 1:
	    q = QVec_Gaussian (qmmesh, qmmesh.Q[i], qwidth, srctp);
	    break;
	case 2:
	    q = QVec_Cosine (qmmesh, qmmesh.Q[i], qwidth, srctp);
	    break;
	}
	qvec.SetRow (i, q);
    }
    SelectMeasurementProfile (pp, mprof, mwidth);
    mvec.New (nM, qmmesh.nlen());
    for (int i = 0; i < nM; i++) {
	RVector m(qmmesh.nlen());
	switch (mprof) {
	case 0:
	    m = QVec_Point (qmmesh, qmmesh.M[i], SRCMODE_NEUMANN);
	    break;
	case 1:
	    m = QVec_Gaussian (qmmesh, qmmesh.M[i], mwidth, SRCMODE_NEUMANN);
	    break;
	case 2:
	    m = QVec_Cosine (qmmesh, qmmesh.M[i], mwidth, SRCMODE_NEUMANN);
	    break;
	}
	mvec.SetRow (i, m);
    }

}

void ReadParams(const char* fname, RVector &muabs, RVector &muscat, RVector &ref)
{
    int size = muabs.Dim();
    cin.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
    ifstream ifs;
    ifs.open (fname);
    xASSERT (ifs.is_open(), "Optical parameter file not found.");

    for(int i =0; i < size; i++)
    {
	try{
		ifs>>muabs[i];
	}
	catch(std::ifstream::failure e){
                cerr << e.what();
        }
	try{
		ifs>>muscat[i];
	}
	catch(std::ifstream::failure e){
                cerr << e.what();
        }
	try{
		ifs>>ref[i];
	}
	catch(std::ifstream::failure e){
                cerr << e.what();
        }
    }
    ifs.close ();
	
}

void ReadSphOrder(const char* fname, IVector &sphorder)
{
    int size = sphorder.Dim();
    cin.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
    ifstream ifs;
    ifs.open (fname);
    xASSERT (ifs.is_open(), "Nodal spherical harmonic order specification file not found.");
    for(int i =0; i < size; i++)
    {
	try{
		ifs>>sphorder[i];
	}
	catch(std::ifstream::failure e){
                cerr << e.what();
        }
    }
    ifs.close ();
	
	
}


void ReadDirections(const char*fname, const int nQ, RVector* &dirVec)
{
    cin.exceptions( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );
    ifstream ifs;
    ifs.open (fname);
    xASSERT (ifs.is_open(), "Direction vectors file not found.");
    for(int i=0; i < nQ; i++)
    {
	for(int j=0; j < dirVec[i].Dim(); j++)
	{
		try{
			ifs>>dirVec[i][j];
		}
		catch(std::ifstream::failure e){
                	cerr << e.what();
        	}
	
	}
    } 
    ifs.close();	
}
void WriteData (const RVector &data, char *fname)
{
    ofstream ofs (fname);
    ofs << setprecision(14);
    ofs << data << endl;
}

void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname)
{
    int q, m, i;
    ofstream ofs (fname);
    for (m = i = 0; m < mesh.nM; m++) {
	for (q = 0; q < mesh.nQ; q++) {
	    if (mesh.Connected (q,m)) ofs << data[i++];
	    else                      ofs << '-';
	    ofs << (q == mesh.nQ-1 ? '\n' : '\t');
	}
    }   
}

void OpenNIM (const char *nimname, const char *meshname, int size)
{
    ofstream ofs (nimname);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = N/A" << endl;
    ofs << "ImageSize = " << size << endl;
    ofs << "EndHeader" << endl;
}

void WriteNIM (const char *nimname, const RVector &img, int size, int no)
{
    ofstream ofs (nimname, ios::app);
    ofs << "Image " << no << endl;
    for (int i = 0; i < size; i++)
        ofs << img[i] << ' ';
    ofs << endl;
}

bool ReadNim (char *nimname, RVector &img)
{
    char cbuf[256];
    int i, imgsize = 0;

    ifstream ifs (nimname);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    do {
        ifs.getline (cbuf, 256);
    } while (strncasecmp (cbuf, "Image", 5));
    for (i = 0; i < imgsize; i++)
        ifs >> img[i];
    return true;
}


