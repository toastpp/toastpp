#include "mathlib.h"
#include <string.h>

using namespace std;

void DisplayInfo();

int main (int argc, char *argv[])
{
    int i, img;
    char cbuf[256];
    char **hdrline;
    int nhdrline;
    //char meshstr[256], soltpstr[256];
    int imgno = -1; // default: last image
    int imsize = 0;
    double *imgbuf;
    bool isnim;
    bool omit_header = false;

    // command line parser
    for (i = 1; i < argc; i++) {
	xASSERT(argv[i][0] == '-', "Error parsing command line");
	switch (argv[i][1]) {
	case 'H':
	case 'h':
	case '?':
	    DisplayInfo();
	    exit(0);
	case 'n':
	    if (sscanf (argv[i+1], "%d", &imgno)) {
		i++;
		break;
	    } else {
		DisplayInfo();
		exit(1);
	    }
	case 'x':
	    omit_header = true;
	    break;
	}
    }

    // check image header
    cin.getline(cbuf, 256);
    xASSERT(!strcmp(cbuf,"NIM") || !strcmp(cbuf,"RIM"), "Unknown image format");
    isnim = (!strcmp (cbuf, "NIM"));

    cerr << "Searching for: ";
    if (imgno < 0) cerr << "last image in stream" << endl;
    else           cerr << "image " << imgno << endl;

    nhdrline = 0;
    do {
	cin.getline(cbuf,256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imsize);
	if (strcasecmp (cbuf, "EndHeader")) {
	    char **tmp = new char*[nhdrline+1];
	    if (nhdrline) {
		memcpy (tmp, hdrline, nhdrline*sizeof(char*));
		delete []hdrline;
	    }
	    hdrline = tmp;
	    hdrline[nhdrline] = new char[strlen(cbuf)+1];
	    strcpy (hdrline[nhdrline], cbuf);
	    nhdrline++;
	}
    } while (!cin.fail() && strcasecmp(cbuf, "EndHeader"));

    xASSERT(!strcasecmp (cbuf, "EndHeader"), "Unexpected end of file");
    xASSERT(imsize > 0, "Could not determine image size");

    imgbuf = new double[imsize];

    for (img = 0;; img++) {
	cin.getline(cbuf,256);
	if (cin.fail() || strncasecmp (cbuf, "Image", 5)) {
	    if (imgno < 0 && img > 0) imgno = --img;
	    else xERROR("Requested image not found in stream");
	} else {
	    for (i = 0; i < imsize; i++)
		cin >> imgbuf[i];
	    cin.getline(cbuf,256);
	}

	if (img == imgno) {
	    if (!omit_header) {
		cout << (isnim ? "NIM" : "RIM") << endl;
		for (i = 0; i < nhdrline; i++)
		    cout << hdrline[i] << endl;
		cout << "EndHeader" << endl;
		cout << "Image 0" << endl;
	    }
	    for (i = 0; i < imsize; i++)
		cout << imgbuf[i] << ' ';
	    cout << endl;
	    cerr << "Extracted image " << img << endl;
	    break;
	}
    }

    if (nhdrline) {
	for (i = 0; i < nhdrline; i++)
	    delete []hdrline[i];
	delete []hdrline;
    }

    delete []imgbuf;
    return 0;
}


void DisplayInfo()
{
    cout << "getnim: Extract a single image from a NIM or RIM image file.\n";
    cout << "\nSyntax: getnim [-n <imgno>] [-x] < infile > outfile\n";
    cout << "-n <imgno>: Image index (0-based). Default: last image in stream.\n";
    cout << "-x: omit header: only write out image data.\n";
    cout << "\ngetnim is a filter which reads a stream of images in NIM\n"
	 << "(nodal image) or RIM (raw image) format from stdin and\n"
	 << "writes the extracted single image to stdout.\n";
}
