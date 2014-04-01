#include "stoastlib.h"
#include "util.h"
#include "supertoast_util.h"

using namespace std;

// ============================================================================

void ReadDataFile (char *fname, RVector &data)
{
    int i, n = data.Dim();
    char c;
    ifstream ifs (fname);
    do {
        ifs >> c;
    } while (ifs.good() && c != '[');
    if (!ifs.good())
	xERROR("Data file not found or invalid");
    for (i = 0; i < n; i++)
        ifs >> data[i];
}

// ============================================================================

void WriteJacobian (const RMatrix *J, const Raster &raster,
    const QMMesh &mesh)
{
    int q, m, i, j, idx, xofs, yofs, dim = raster.Dim();
    double vmin, vmax, scal;

    if (dim != 2) { // for now, just dump the whole thing
	for (i = 0; i < J->nRows(); i++) {
	    char cbuf[256];
	    sprintf (cbuf, "pmdf_%03d.dat", i);
	    ofstream ofs (cbuf);
	    ofs << J->Row(i) << endl;
	}
	return; // this function only handles 2D problems
    }

    IVector bdim = raster.BDim();
    int blen = raster.BLen();
    int slen = raster.SLen();
    int imgw = mesh.nM*bdim[0];
    int imgh = mesh.nQ*bdim[1];
    IVector idim(2);
    idim[0] = imgw, idim[1] = imgh;
    RVector img(imgw*imgh);
    RVector Ji(J->nCols());

    // 1: real, mua
    for (q = idx = 0; q < mesh.nQ; q++) {
        yofs = q*bdim[1];
        for (m = 0; m < mesh.nM; m++) {
	    if (!mesh.Connected (q, m)) continue;
	    xofs = m*bdim[0];
	    Ji = J->Row(idx);
	    RVector sol(Ji, 0, slen);
	    RVector p(blen);
	    raster.Map_SolToBasis (sol, p);
	    ImageScale (p, vmin, vmax);
	    scal = 1.0/(vmax-vmin);
	    for (j = 0; j < bdim[1]; j++)
	        for (i = 0; i < bdim[0]; i++)
		    img[(yofs+j)*imgw + (xofs+i)] = (p[j*bdim[0]+i]-vmin)*scal;
	    idx++;
	}
    }
    vmin = 0.0, vmax = 1.0;
    WritePPM (img, idim, &vmin, &vmax, "j_re_mua.ppm");
    cout << "Jacobian (real, mua) written to j_re_mua.ppm" << endl;

    // 2: real, kappa
    if (J->nCols() >= slen*2) { // J contains kappa
        for (q = idx = 0; q < mesh.nQ; q++) {
	    yofs = q*bdim[1];
	    for (m = 0; m < mesh.nM; m++) {
	        if (!mesh.Connected (q, m)) continue;
		xofs = m*bdim[0];
		Ji = J->Row(idx);
		RVector sol(Ji, slen, slen);
		RVector p(blen);
		raster.Map_SolToBasis (sol, p);
		ImageScale (p, vmin, vmax);
		scal = 1.0/(vmax-vmin);
		for (j = 0; j < bdim[1]; j++)
		    for (i = 0; i < bdim[0]; i++)
		        img[(yofs+j)*imgw + (xofs+i)] =
			    (p[j*bdim[0]+i]-vmin)*scal;
		idx++;
	    }
	}
	vmin = 0.0, vmax = 1.0;
	WritePPM (img, idim, &vmin, &vmax, "j_re_kappa.ppm");
	cout << "Jacobian (real, kappa) written to j_re_kappa.ppm" << endl;
    }

    if (J->nRows() >= mesh.nQM*2) {
        // 3: imag, mua
        for (q = 0, idx = mesh.nQM; q < mesh.nQ; q++) {
	    yofs = q*bdim[1];
	    for (m = 0; m < mesh.nM; m++) {
	        if (!mesh.Connected (q, m)) continue;
		xofs = m*bdim[0];
		Ji = J->Row(idx);
		RVector sol(Ji, 0, slen);
		RVector p(blen);
		raster.Map_SolToBasis (sol, p);
		ImageScale (p, vmin, vmax);
		scal = 1.0/(vmax-vmin);
		for (j = 0; j < bdim[1]; j++)
		    for (i = 0; i < bdim[0]; i++)
		        img[(yofs+j)*imgw + (xofs+i)] =
			    (p[j*bdim[0]+i]-vmin)*scal;
		idx++;
	    }
	}
	vmin = 0.0, vmax = 1.0;
	WritePPM (img, idim, &vmin, &vmax, "j_im_mua.ppm");
	cout << "Jacobian (imag, mua) written to j_im_mua.ppm" << endl;

	// 4: real, kappa
	if (J->nCols() >= slen*2) { // J contains kappa
	    for (q = 0, idx = mesh.nQM; q < mesh.nQ; q++) {
	        yofs = q*bdim[1];
		for (m = 0; m < mesh.nM; m++) {
		    if (!mesh.Connected (q, m)) continue;
		    xofs = m*bdim[0];
		    Ji = J->Row(idx);
		    RVector sol(Ji, slen, slen);
		    RVector p(blen);
		    raster.Map_SolToBasis (sol, p);
		    ImageScale (p, vmin, vmax);
		    scal = 1.0/(vmax-vmin);
		    for (j = 0; j < bdim[1]; j++)
		        for (i = 0; i < bdim[0]; i++)
			    img[(yofs+j)*imgw + (xofs+i)] =
			        (p[j*bdim[0]+i]-vmin)*scal;
		    idx++;
		}
	    }
	    vmin = 0.0, vmax = 1.0;
	    WritePPM (img, idim, &vmin, &vmax, "j_im_kappa.ppm");
	    cout << "Jacobian (imag, kappa) written to j_im_kappa.ppm" << endl;
	}
    }
}

