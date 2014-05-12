// mesh2nim.cc
// Convert mua and mus parameters stored in mesh to nim files

#include <mathlib.h>
#include <felib.h>
#include <fstream>

using namespace std;

void PrintUsage();
void WriteNIM (const char *nimfile, const char *meshfile, const RVector &vec);

int main (int argc, char *argv[])
{
    char *meshfile;
    Mesh mesh;
    int n;

    if (argc < 2) {
	PrintUsage ();
	exit (1);
    }
    meshfile = argv[1];

    cout << "Reading " << meshfile << " ..." << endl;
    ifstream ifs (meshfile);
    if (!ifs) {
	cerr << " FAILED!" << endl;
	exit (1);
    }
    ifs >> mesh;
    mesh.Setup();
    n = mesh.nlen();

    RVector mua(n);
    RVector mus(n);
    mua = mesh.plist.Mua();
    mus = mesh.plist.Mus();

    WriteNIM ("mua.nim", meshfile, mua);
    cout << "MUA image written to mua.nim" << endl;
    WriteNIM ("mus.nim", meshfile, mus);
    cout << "MUS image written to mus.nim" << endl;
    return 0;
}

void PrintUsage()
{
    cerr << "Usage: meshinfo meshfile" << endl;
}

void WriteNIM (const char *nimfile, const char *meshfile, const RVector &vec)
{
    ofstream ofs (nimfile);
    ofs << "NIM\n";
    ofs << "Mesh = " << meshfile << endl;
    ofs << "SolutionType = N/A" << endl;
    ofs << "ImageSize = " << vec.Dim() << endl;
    ofs << "EndHeader" << endl;
    ofs << "Image 0" << endl;
    for (int i = 0; i < vec.Dim(); i++)
        ofs << vec[i] << ' ';
    ofs << endl;
}
