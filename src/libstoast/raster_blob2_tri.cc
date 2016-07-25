// Implementation of those routines in Raster_Blob2 specific to triangular
// meshes

#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"
#include "raster_blob2.h"
#include "tri_qr.h"
#include "sh.h"

const int nprog = 50; // length of progress bars

// MC integral of two blob basis functions, with distance basis_dst between
// their centres (basis_dst in units of grid spacing)
double Raster_Blob2::MC_integral_2D(double basis_dst) const
{
    if (basis_dst >= 2.0*sup) return 0.0;

    double px, py, r1, r2, v1, v2, sum = 0.0;
    //const int nsample = 100000000;
    //for (int i = 0; i < nsample; i++) {
    const int nsample = 10000;
    for (int i = 0; i < nsample; i++) {
	px = -sup + (i+0.5)/(double)nsample * 2.0*sup;
	for (int j = 0; j < nsample; j++) {
	    py = -sup + (j+0.5)/(double)nsample * 2.0*sup;
	    //px = (drand48()-0.5)*2.0*sup;
	    //py = (drand48()-0.5)*2.0*sup;
	    r1 = sqrt(px*px+py*py);
	    v1 = RadValue(r1);
	    if (basis_dst) {
		px = basis_dst-px;
		r2 = sqrt(px*px+py*py);
		v2 = RadValue(r2);
	    } else v2 = v1;
	    sum += v1*v2;
	}
    }
    return sum/((double)nsample*(double)nsample)*4.0*sup*sup;
}

int QRule_tri_dummy (const double **wght, const Point **absc)
{
    const double x[20][20] = {
	{2.8516e-05,0.00014952,0.00036423,0.00066767,0.0010527,0.0015104,0.0020299,0.0025991,0.0032047,0.0038324,
	 0.0044676 ,0.0050953 ,0.0057009 ,0.0062701 ,0.0067897,0.0072473,0.0076324,0.0079358,0.0081505,0.0082715},
	{9.4973e-05,0.00049796,0.0012131 ,0.0022236 ,0.0035061,0.0050303,0.0067606,0.0086563,0.010673 ,0.012764,
	 0.014879  ,0.01697   ,0.018987  ,0.020882  ,0.022613 ,0.024137 ,0.025419 ,0.02643  ,0.027145 ,0.027548},
	{0.00019767,0.0010364 ,0.0025248 ,0.0046282 ,0.0072973,0.01047  ,0.014071 ,0.018017 ,0.022214 ,0.026566,
	 0.030969  ,0.03532   ,0.039518  ,0.043463  ,0.047065 ,0.050237 ,0.052906 ,0.05501  ,0.056498 ,0.057337},
	{0.00033431,0.0017528 ,0.00427   ,0.0078273 ,0.012341 ,0.017707 ,0.023797 ,0.03047  ,0.03757  ,0.044929,
         0.052375  ,0.059734  ,0.066834  ,0.073507  ,0.079597 ,0.084963 ,0.089477 ,0.093034 ,0.095551 ,0.09697},
	{0.00050183,0.0026312 ,0.0064097 ,0.01175   ,0.018526 ,0.02658  ,0.035722 ,0.045739 ,0.056396 ,0.067443,
	 0.07862   ,0.089667  ,0.10032   ,0.11034   ,0.11948  ,0.12754  ,0.13431  ,0.13965  ,0.14343  ,0.14556},
	{0.00069649,0.0036519 ,0.008896  ,0.016307  ,0.025712 ,0.03689  ,0.049579 ,0.063482 ,0.078273 ,0.093604,
	 0.10912   ,0.12445   ,0.13924   ,0.15314   ,0.16583  ,0.17701  ,0.18642  ,0.19383  ,0.19907  ,0.20203},
	{0.00091395,0.004792  ,0.011674  ,0.021399  ,0.03374  ,0.048408 ,0.065059 ,0.083302 ,0.10271  ,0.12283,
	 0.14319   ,0.16331   ,0.18271   ,0.20096   ,0.21761  ,0.23228  ,0.24462  ,0.25434  ,0.26122  ,0.2651},
	{0.0011493 ,0.0060262 ,0.01468   ,0.02691   ,0.04243  ,0.060876 ,0.081815 ,0.10476  ,0.12916  ,0.15446,
	 0.18007   ,0.20537   ,0.22977   ,0.25272   ,0.27365  ,0.2921   ,0.30762  ,0.31985  ,0.3285   ,0.33338},
	{0.0013974 ,0.0073269 ,0.017849  ,0.032718  ,0.051588 ,0.074015 ,0.099474 ,0.12737  ,0.15704  ,0.1878,
	 0.21893   ,0.24969   ,0.27937   ,0.30726   ,0.33272  ,0.35515  ,0.37402  ,0.38889  ,0.39941  ,0.40534},
	{0.0016526 ,0.008665  ,0.021108  ,0.038694  ,0.061009 ,0.087532 ,0.11764  ,0.15063  ,0.18572  ,0.2221,
	 0.25891   ,0.29529   ,0.33039   ,0.36338   ,0.39348  ,0.42001  ,0.44232  ,0.45991  ,0.47235  ,0.47936},
	{0.0019093 ,0.010011  ,0.024386  ,0.044703  ,0.070484 ,0.10113  ,0.13591  ,0.17402  ,0.21457  ,0.25659,
	 0.29912   ,0.34115   ,0.38169   ,0.41981   ,0.45459  ,0.48523  ,0.51101  ,0.53133  ,0.5457   ,0.55381},
	{0.0021616 ,0.011334  ,0.027609  ,0.050611  ,0.079799 ,0.11449  ,0.15387  ,0.19702  ,0.24292  ,0.29051,
	 0.33866   ,0.38624   ,0.43214   ,0.47529   ,0.51467  ,0.54936  ,0.57855  ,0.60155  ,0.61783  ,0.627},
	{0.002404  ,0.012605  ,0.030706  ,0.056286  ,0.088748 ,0.12733  ,0.17113  ,0.21911  ,0.27017  ,0.32309,
	 0.37663   ,0.42955   ,0.4806    ,0.52859   ,0.57239  ,0.61097  ,0.64343  ,0.66901  ,0.68711  ,0.69732},
	{0.0026311 ,0.013795  ,0.033606  ,0.061603  ,0.097131 ,0.13936  ,0.18729  ,0.23981  ,0.29568  ,0.3536,
	 0.41221   ,0.47012   ,0.526     ,0.57852   ,0.62645  ,0.66868  ,0.70421  ,0.7322   ,0.75201  ,0.76318},
	{0.0028377 ,0.014879  ,0.036245  ,0.066441  ,0.10476  ,0.1503   ,0.202    ,0.25864  ,0.31891  ,0.38137,
	 0.44458   ,0.50705   ,0.56731   ,0.62395   ,0.67565  ,0.72119  ,0.75951  ,0.78971  ,0.81107  ,0.82312},
	{0.0030193 ,0.015831  ,0.038565  ,0.070693  ,0.11146  ,0.15992  ,0.21493  ,0.2752   ,0.33931  ,0.40578,
	 0.47303   ,0.5395    ,0.60361   ,0.66388   ,0.71889  ,0.76735  ,0.80812  ,0.84025  ,0.86298  ,0.87579},
	{0.0031718 ,0.016631  ,0.040513  ,0.074264  ,0.11709  ,0.168    ,0.22578  ,0.2891   ,0.35645  ,0.42627,
	 0.49692   ,0.56675   ,0.6341    ,0.69742   ,0.7552   ,0.80611  ,0.84894  ,0.88269  ,0.90657  ,0.92003},
	{0.0032918 ,0.01726   ,0.042045  ,0.077073  ,0.12152  ,0.17435  ,0.23433  ,0.30004  ,0.36994  ,0.4424,
	 0.51573   ,0.58819   ,0.65809   ,0.7238    ,0.78377  ,0.83661  ,0.88106  ,0.91608  ,0.94087  ,0.95484},
	{0.0033767 ,0.017705  ,0.043129  ,0.079059  ,0.12465  ,0.17885  ,0.24036  ,0.30777  ,0.37947  ,0.4538,
	 0.52902   ,0.60335   ,0.67505   ,0.74245   ,0.80397  ,0.85816  ,0.90376  ,0.93969  ,0.96511  ,0.97944},
	{0.0034244 ,0.017955  ,0.043739  ,0.080178  ,0.12642  ,0.18138  ,0.24377  ,0.31212  ,0.38484  ,0.46022,
	 0.5365    ,0.61188   ,0.6846    ,0.75296   ,0.81535  ,0.87031  ,0.91655  ,0.95298  ,0.97877  ,0.9933}};

    const double y[20][20] = {
	{0.0082715,0.0081505,0.0079358,0.0076324,0.0072473,0.0067897,0.0062701 ,0.0057009 ,0.0050953 ,0.0044676,
	 0.0038324,0.0032047,0.0025991,0.0020299,0.0015104,0.0010527,0.00066767,0.00036423,0.00014952,2.8516e-05},
	{0.027548 ,0.027145 ,0.02643  ,0.025419 ,0.024137 ,0.022613 ,0.020882  ,0.018987  ,0.01697   ,0.014879,
	 0.012764 ,0.010673 ,0.0086563,0.0067606,0.0050303,0.0035061,0.0022236 ,0.0012131 ,0.00049796,9.4973e-05},
	{0.057337 ,0.056498 ,0.05501  ,0.052906 ,0.050237 ,0.047065 ,0.043463  ,0.039518  ,0.03532   ,0.030969,
	 0.026566 ,0.022214 ,0.018017 ,0.014071 ,0.01047  ,0.0072973,0.0046282 ,0.0025248 ,0.0010364 ,0.00019767},
	{0.09697  ,0.095551 ,0.093034 ,0.089477 ,0.084963 ,0.079597 ,0.073507  ,0.066834  ,0.059734  ,0.052375,
	 0.044929 ,0.03757  ,0.03047  ,0.023797 ,0.017707 ,0.012341 ,0.0078273 ,0.00427   ,0.0017528 ,0.00033431},
	{0.14556  ,0.14343  ,0.13965  ,0.13431  ,0.12754  ,0.11948  ,0.11034   ,0.10032   ,0.089667  ,0.07862,
	 0.067443 ,0.056396 ,0.045739 ,0.035722 ,0.02658  ,0.018526 ,0.01175   ,0.0064097 ,0.0026312 ,0.00050183},
	{0.20203  ,0.19907  ,0.19383  ,0.18642  ,0.17701  ,0.16583  ,0.15314   ,0.13924   ,0.12445   ,0.10912,
	 0.093604 ,0.078273 ,0.063482 ,0.049579 ,0.03689  ,0.025712 ,0.016307  ,0.008896  ,0.0036519 ,0.00069649},
	{0.2651   ,0.26122  ,0.25434  ,0.24462  ,0.23228  ,0.21761  ,0.20096   ,0.18271   ,0.16331   ,0.14319,
	 0.12283  ,0.10271  ,0.083302 ,0.065059 ,0.048408 ,0.03374  ,0.021399  ,0.011674  ,0.004792  ,0.00091395},
	{0.33338  ,0.3285   ,0.31985  ,0.30762  ,0.2921   ,0.27365  ,0.25272   ,0.22977   ,0.20537   ,0.18007,
	 0.15446  ,0.12916  ,0.10476  ,0.081815 ,0.060876 ,0.04243  ,0.02691   ,0.01468   ,0.0060262 ,0.0011493},
	{0.40534  ,0.39941  ,0.38889  ,0.37402  ,0.35515  ,0.33272  ,0.30726   ,0.27937   ,0.24969   ,0.21893,
	 0.1878   ,0.15704  ,0.12737  ,0.099474 ,0.074015 ,0.051588 ,0.032718  ,0.017849  ,0.0073269 ,0.0013974},
	{0.47936  ,0.47235  ,0.45991  ,0.44232  ,0.42001  ,0.39348  ,0.36338   ,0.33039   ,0.29529   ,0.25891,
	 0.2221   ,0.18572  ,0.15063  ,0.11764  ,0.087532 ,0.061009 ,0.038694  ,0.021108  ,0.008665  ,0.0016526},
	{0.55381  ,0.5457   ,0.53133  ,0.51101  ,0.48523  ,0.45459  ,0.41981   ,0.38169   ,0.34115   ,0.29912,
	 0.25659  ,0.21457  ,0.17402  ,0.13591  ,0.10113  ,0.070484 ,0.044703  ,0.024386  ,0.010011  ,0.0019093},
	{0.627    ,0.61783  ,0.60155  ,0.57855  ,0.54936  ,0.51467  ,0.47529   ,0.43214   ,0.38624   ,0.33866,
	 0.29051  ,0.24292  ,0.19702  ,0.15387  ,0.11449  ,0.079799 ,0.050611  ,0.027609  ,0.011334  ,0.0021616},
	{0.69732  ,0.68711  ,0.66901  ,0.64343  ,0.61097  ,0.57239  ,0.52859   ,0.4806    ,0.42955   ,0.37663,
	 0.32309  ,0.27017  ,0.21911  ,0.17113  ,0.12733  ,0.088748 ,0.056286  ,0.030706  ,0.012605  ,0.002404},
	{0.76318  ,0.75201  ,0.7322   ,0.70421  ,0.66868  ,0.62645  ,0.57852   ,0.526     ,0.47012   ,0.41221,
	 0.3536   ,0.29568  ,0.23981  ,0.18729  ,0.13936  ,0.097131 ,0.061603  ,0.033606  ,0.013795  ,0.0026311},
	{0.82312  ,0.81107  ,0.78971  ,0.75951  ,0.72119  ,0.67565  ,0.62395   ,0.56731   ,0.50705   ,0.44458,
	 0.38137  ,0.31891  ,0.25864  ,0.202    ,0.1503   ,0.10476  ,0.066441  ,0.036245  ,0.014879  ,0.0028377},
	{0.87579  ,0.86298  ,0.84025  ,0.80812  ,0.76735  ,0.71889  ,0.66388   ,0.60361   ,0.5395    ,0.47303,
	 0.40578  ,0.33931  ,0.2752   ,0.21493  ,0.15992  ,0.11146  ,0.070693  ,0.038565  ,0.015831  ,0.0030193},
	{0.92003  ,0.90657  ,0.88269  ,0.84894  ,0.80611  ,0.7552   ,0.69742   ,0.6341    ,0.56675   ,0.49692,
	 0.42627  ,0.35645  ,0.2891   ,0.22578  ,0.168    ,0.11709  ,0.074264  ,0.040513  ,0.016631  ,0.0031718},
	{0.95484  ,0.94087  ,0.91608  ,0.88106  ,0.83661  ,0.78377  ,0.7238    ,0.65809   ,0.58819   ,0.51573,
	 0.4424   ,0.36994  ,0.30004  ,0.23433  ,0.17435  ,0.12152  ,0.077073  ,0.042045  ,0.01726   ,0.0032918},
	{0.97944  ,0.96511  ,0.93969  ,0.90376  ,0.85816  ,0.80397  ,0.74245   ,0.67505   ,0.60335   ,0.52902,
	 0.4538   ,0.37947  ,0.30777  ,0.24036  ,0.17885  ,0.12465  ,0.079059  ,0.043129  ,0.017705  ,0.0033767},
	{0.9933   ,0.97877  ,0.95298  ,0.91655  ,0.87031  ,0.81535  ,0.75296   ,0.6846    ,0.61188   ,0.5365,
	 0.46022  ,0.38484  ,0.31212  ,0.24377  ,0.18138  ,0.12642  ,0.080178  ,0.043739  ,0.017955  ,0.0034244}};

    const double wx[20] = {
	0.00011538, 0.00068306, 0.0020115, 0.0043232, 0.0077277, 0.012204, 0.017597, 0.023625, 0.029902, 0.035965,
	0.041318,   0.045471,   0.047985,  0.048516,  0.046842,  0.04289,  0.03675,  0.028668, 0.019035, 0.0083708};
    const double wy[20] = {
	0.008807,   0.020301,   0.031336,  0.041638,  0.050965,  0.059097, 0.065844, 0.071048, 0.074586, 0.076377,
	0.076377,   0.074586,   0.071048,  0.065844,  0.059097,  0.050965, 0.041638, 0.031336, 0.020301, 0.008807};

    static Point a[400];
    static double w[400];
    static bool need_setup = true;

    if (need_setup) {
	for (int j = 0; j < 20; j++) {
	    for (int i = 0; i < 20; i++) {
		a[i+j*20].New(2);
		a[i+j*20][0] = x[i][j];
		a[i+j*20][1] = y[i][j];
		w[i+j*20] = wx[i]*wy[j];
	    }
	}
	need_setup = false;
    }
    *absc = a;
    *wght = w;
	    
    return 400;
}

// ===========================================================================

RCompRowMatrix *Raster_Blob2::CreateBvv_tri () const
{
    const double eps = 1e-10;
    const int subgrd = 200;
    int i, i0, j0, i1, j1, bidx0, bidx1, prog, count;
    double xcnt0, ycnt0, xcnt1, ycnt1, dstx, dsty;
    double xmin, ymin, xmax, ymax, x, y, r, r0, r1, v0, v1, bvvi;
    double dx0_2, dy0_2, dx1_2, dy1_2, area;
    double xstep = 1.0/((bdim[0]-1.0)*subgrd);
    double ystep = 1.0/((bdim[1]-1.0)*subgrd);
    double dx = grid[0];
    double dy = grid[1];
    double radlimit2 = sup*sup;
    double radlimit22 = radlimit2*4.0;
    IVector npx(blen_pad);
    struct OVERLAP {
	double dst;
	double v;
    } *overlap;
    int noverlap = 0;
    
    std::cout << "Bvv: pass 1 [" << std::flush; prog = 0;
    for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	ycnt0 = bbmin[1]+j0*dy;
	for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
	    xcnt0 = bbmin[0]+i0*dx;
	    bidx0 = (i0+npad) + (j0+npad)*bdim_pad[0];
	    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
		ycnt1 = bbmin[1]+j1*dy;
		for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
		    xcnt1 = bbmin[0]+i1*dx;
		    dstx = xcnt0-xcnt1;
		    dsty = ycnt0-ycnt1;
		    if (dstx*dstx + dsty*dsty < radlimit22)
			npx[bidx0]++;
		}
	    }
	}
	while (((j0+1)*nprog)/bdim_pad[1] > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    int *rowptr = new int[blen_pad+1];
    rowptr[0] = 0;
    for (i = 0; i < blen_pad; i++)
	rowptr[i+1] = rowptr[i] + npx[i];
    int nz = rowptr[blen_pad];
    int *colidx = new int[nz];
    npx.Clear();

    std::cout << "Bvv: pass 2 [" << std::flush; prog = 0;
    for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	ycnt0 = bbmin[1]+j0*dy;
	for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
	    xcnt0 = bbmin[0]+i0*dx;
	    bidx0 = (i0+npad) + (j0+npad)*bdim_pad[0];
	    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
		ycnt1 = bbmin[1]+j1*dy;
		for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
		    xcnt1 = bbmin[0]+i1*dx;
		    bidx1 = (i1+npad) + (j1+npad)*bdim_pad[0];
		    dstx = xcnt0-xcnt1;
		    dsty = ycnt0-ycnt1;
		    if (dstx*dstx + dsty*dsty < radlimit22)
			colidx[rowptr[bidx0]+npx[bidx0]++] = bidx1;
		}
	    }
	}
	while (((j0+1)*nprog)/bdim_pad[1] > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    RCompRowMatrix *bvv = new RCompRowMatrix (blen_pad, blen_pad,
					      rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    std::cout << "Bvv: pass 3 [" << std::flush; prog = 0;
    for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	ycnt0 = bbmin[1]+j0*dy;
	for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
	    xcnt0 = bbmin[0]+i0*dx;
	    bidx0 = (i0+npad) + (j0+npad)*bdim_pad[0];
	    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
		ycnt1 = bbmin[1]+j1*dy;
		for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
		    xcnt1 = bbmin[0]+i1*dx;
		    bidx1 = (i1+npad) + (j1+npad)*bdim_pad[0];
		    dstx = xcnt0-xcnt1;
		    dsty = ycnt0-ycnt1;
		    if (dstx*dstx + dsty*dsty < radlimit22) {
			r = hypot(dstx,dsty);
			// check if we already have an entry
			for (i = 0; i < noverlap; i++) {
			    if (fabs(overlap[i].dst-r) < eps) {
				bvvi = overlap[i].v;
				break;
			    }
			}
			if (i == noverlap) { // not present - calculate
			    xmin = max(xcnt0,xcnt1)-sup;
			    xmax = min(xcnt0,xcnt1)+sup;
			    ymin = max(ycnt0,ycnt1)-sup;
			    ymax = min(ycnt0,ycnt1)+sup;
			    bvvi = 0.0;
			    count = 0;
			    for (x = xmin; x <= xmax; x += xstep) {
				dstx = x-xcnt0;
				dx0_2 = dstx*dstx;
				dstx = x-xcnt1;
				dx1_2 = dstx*dstx;
				for (y = ymin; y <= ymax; y += ystep) {
				    dsty = y-ycnt0;
				    dy0_2 = dsty*dsty;
				    if (dx0_2 + dy0_2 >= radlimit2) continue;
				    dsty = y-ycnt1;
				    dy1_2 = dsty*dsty;
				    if (dx1_2 + dy1_2 >= radlimit2) continue;
				    r0 = sqrt(dx0_2+dy0_2);
				    r1 = sqrt(dx1_2+dy1_2);
				    v0 = RadValue(r0);
				    v1 = RadValue(r1);
				    bvvi += v0*v1;
				    count++;
				}
			    }
			    // area of overlap
			    if (r < eps) {
				area = sup*sup*Pi;
			    } else {
				double theta = 2.0*acos(r/(2.0*sup));
				area = sup*sup*0.5*(theta - sin(theta));
			    }
			    if (count)
				bvvi *= area/(double)count;
			    OVERLAP *tmp = new OVERLAP[noverlap+1];
			    if (noverlap) {
				memcpy(tmp,overlap,noverlap*sizeof(OVERLAP));
				delete []overlap;
			    }
			    overlap = tmp;
			    overlap[noverlap].dst = r;
			    overlap[noverlap].v = bvvi;
			    noverlap++;
			}
			(*bvv)(bidx0,bidx1) = bvvi;
		    }
		}
	    }
	}
	while (((j0+1)*nprog)/bdim_pad[1] > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;
    //(*bvv) *= xstep*ystep;
    
    // diagonal conditioning
    RVector d = bvv->Diag();
    double dmean = mean(d);
    for (i = 0; i < blen; i++)
	(*bvv)(i,i) += dmean*dgscale;

    bvv->Shrink();
    return bvv;
}

#ifdef UNDEF
RCompRowMatrix *Raster_Blob2::CreateBvv_tri () const
{
    int i, j, i0, j0, i1, j1, k, ii, jj, nv, nv0, nv1, m, el, nd, idx_i, idx_j;
    int prog, nel = meshptr->elen(), n = meshptr->nlen();
    bool intersect;
    double rx, ry, djac, v, dstx, dsty;
    RVector fun;

    // grid range and spacing (assumes regular arrangement of blob
    // basis functions)
    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double radlimit2 = sup*sup;
    double radlimit22 = radlimit2*4.0;

    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    //int np = QRule_tri_9_19 (&wght, &absc);
    int np = QRule_tri_dummy (&wght, &absc);

    int *npx = new int[blen_pad];
    for (i = 0; i < blen_pad; i++) npx[i] = 0;

    bool *have_diag = new bool[blen_pad];
    for (i = 0; i < blen_pad; i++) have_diag[i] = false;

    // pass 1: determine matrix fill structure
    double px0, py0, px1, py1;
    std::cout << "Bvv: pass 1 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	    py0 = bbmin[1] + j0*dy;
	    for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
		px0 = bbmin[0] + i0*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px0-meshptr->nlist[nd][0];
		    ry = py0-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
			py1 = bbmin[1] + j1*dy;
			for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
			    px1 = bbmin[0] + i1*dx;
			    dstx = px1-px0;
			    dsty = py1-py0;
			    if (dstx*dstx + dsty*dsty >= radlimit22)
				continue;
			    intersect = false;
			    for (k = 0; k < pel->nNode(); k++) {
				nd = pel->Node[k];
				rx = px1-meshptr->nlist[nd][0];
				ry = py1-meshptr->nlist[nd][1];
				if (rx*rx + ry*ry < radlimit2) {
				    intersect = true;
				    break;
				}
			    }
			    if (!intersect) continue;
			    idx_i = (i0+npad) + (j0+npad)*bdim_pad[0];
			    if (idx_i == (i1+npad) + (j1+npad)*bdim_pad[0])
				have_diag[idx_i] = true;
			    npx[idx_i]++;
			}
		    }
		}
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    for (i = 0; i < blen_pad; i++) {
	if (!have_diag[i]) // reserve space for diag element
	npx[i]++;
	have_diag[i] = false;
    }

    int *rowptr = new int[blen_pad+1];
    rowptr[0] = 0;
    for (i = 0; i < blen_pad; i++)
	rowptr[i+1] = rowptr[i]+npx[i];
    int nz = rowptr[blen_pad];
    int *colidx = new int[nz];
    for (i = 0; i < blen_pad; i++)
	npx[i] = 0;

    // pass 2
    std::cout << "Bvv: pass 2 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	    py0 = bbmin[1] + j0*dy;
	    for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
		px0 = bbmin[0] + i0*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px0-meshptr->nlist[nd][0];
		    ry = py0-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		idx_i = (i0+npad) + (j0+npad)*bdim_pad[0];
		if (intersect) {
		    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
			py1 = bbmin[1] + j1*dy;
			for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
			    px1 = bbmin[0] + i1*dx;
			    dstx = px1-px0;
			    dsty = py1-py0;
			    if (dstx*dstx + dsty*dsty >= radlimit22)
				continue;
			    intersect = false;
			    for (k = 0; k < pel->nNode(); k++) {
				nd = pel->Node[k];
				rx = px1-meshptr->nlist[nd][0];
				ry = py1-meshptr->nlist[nd][1];
				if (rx*rx + ry*ry < radlimit2) {
				    intersect = true;
				    break;
				}
			    }
			    if (!intersect) continue;
			    idx_j = (i1+npad) + (j1+npad)*bdim_pad[0];
			    if (idx_i == idx_j)
				have_diag[idx_i] = true;
			    colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
			}
		    }
		}
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    for (i = 0; i < blen_pad; i++) {
	if (!have_diag[i])
	    colidx[rowptr[i] + npx[i]++] = i;
    }
    delete []have_diag;

    RCompRowMatrix *bvv = new RCompRowMatrix (
        blen_pad, blen_pad, rowptr, colidx);

    // pass 3: fill the matrix
    std::cout << "Bvv: pass 3 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	const int npoly = 10;
	Point poly0[npoly];
	Point poly1[npoly];
	Point poly[npoly];
	Triangle3 t3;
	NodeList n3(3);
	for (i = 0; i < 3; i++) {
	    t3.Node[i] = i;
	    n3[i].New(2);
	}

	for (idx_i = 0; idx_i < blen_pad; idx_i++) {
	    i0 = idx_i % bdim_pad[0] - npad;
	    j0 = idx_i / bdim_pad[0] - npad;
	    px0 = bbmin[0] + i0*dx;
	    py0 = bbmin[1] + j0*dy;
	    intersect = false;
	    for (k = 0; k < pel->nNode(); k++) {
		nd = pel->Node[k];
		rx = px0-meshptr->nlist[nd][0];
		ry = py0-meshptr->nlist[nd][1];
		if (rx*rx + ry*ry < radlimit2) {
		    intersect = true;
		    break;
		}
	    }
	    if (!intersect) continue;
	    nv0 = SutherlandHodgman (el, px0, py0, poly0, npoly);
	    for (idx_j = 0; idx_j <= idx_i; idx_j++) {
		i1 = idx_j % bdim_pad[0] - npad;
		j1 = idx_j / bdim_pad[0] - npad;
		px1 = bbmin[0] + i1*dx;
		py1 = bbmin[1] + j1*dy;
		dstx = px1-px0;
		dsty = py1-py0;
		if (dstx*dstx + dsty*dsty >= radlimit22)
		    continue;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px1-meshptr->nlist[nd][0];
		    ry = py1-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (!intersect) continue;
		nv1 = SutherlandHodgman (el, px1, py1, poly1, npoly);
		nv  = SutherlandHodgman (poly0, nv0, poly1, nv1,
					 poly, npoly);

		// split into nv-2 triangles
		for (k = 0; k < nv-2; k++) {
		    n3[0][0] = poly[0][0];
		    n3[0][1] = poly[0][1];
		    n3[1][0] = poly[k+1][0];
		    n3[1][1] = poly[k+1][1];
		    n3[2][0] = poly[k+2][0];
		    n3[2][1] = poly[k+2][1];
		    t3.Initialise(n3);
		    
		    for (m = 0; m < np; m++) {
			djac = t3.DetJ(absc[m], &n3);
			Point glob = t3.Global(n3, absc[m]);
			v = wght[m] * djac *
			    Value_nomask(glob,idx_i,false) *
			    Value_nomask(glob,idx_j,false);
			(*bvv)(idx_i,idx_j) += v;
			if (idx_i != idx_j)
			    (*bvv)(idx_j,idx_i) += v;
		    }
		}
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    delete []rowptr;
    delete []colidx;
    delete []npx;

    // diagonal conditioning
    RVector d = bvv->Diag();
    double dmean = mean(d);
    for (i = 0; i < blen_pad; i++)
	(*bvv)(i,i) += dmean*dgscale;

    bvv->Shrink();
    return bvv;
}
#endif

// ===========================================================================

RCompRowMatrix *Raster_Blob2::CreateBuv_tri () const
{
    int i, j, k, ii, nv, m, el, nd, idx_i, idx_j, prog;
    int nel = meshptr->elen(), n = meshptr->nlen();
    bool intersect;
    double rx, ry, djac, v;
    RVector fun;

    // grid range and spacing (assumes regular arrangement of blob
    // basis functions)
    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double radlimit2 = sup*sup;
    
    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    int np = QRule_tri_4_6 (&wght, &absc);

    int *npx = new int[n];
    for (i = 0; i < n; i++) npx[i] = 0;

    // pass 1: determine matrix fill structure
    double px, py;
    std::cout << "Buv: pass 1 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j = -npad; j < bdim[1]+npad; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = -npad; i < bdim[0]+npad; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    for (k = 0; k < pel->nNode(); k++)
			npx[pel->Node[k]]++;
		}
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    int *rowptr = new int[n+1];
    rowptr[0] = 0;
    for (i = 0; i < n; i++)
	rowptr[i+1] = rowptr[i]+npx[i];
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (i = 0; i < n; i++)
	npx[i] = 0;

    // pass 2
    std::cout << "Buv: pass 2 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j = -npad; j < bdim[1]+npad; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = -npad; i < bdim[0]+npad; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    for (k = 0; k < pel->nNode(); k++) {
			nd = pel->Node[k];
			colidx[rowptr[nd]+npx[nd]++] =
			    (i+npad) + (j+npad)*bdim_pad[0];
		    }
		}
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    RCompRowMatrix *buv = new RCompRowMatrix (n, blen_pad, rowptr, colidx);

    // pass 3: fill the matrix
    std::cout << "Buv: pass 3 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	const int npoly = 10;
	Point poly[npoly];
	Triangle3 t3;
	NodeList n3(3);
	for (i = 0; i < 3; i++) {
	    t3.Node[i] = i;
	    n3[i].New(2);
	}

	for (j = -npad; j < bdim[1]+npad; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = -npad; i < bdim[0]+npad; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    nv = SutherlandHodgman (el, px, py, poly, npoly);
		    idx_j = (i+npad) + (j+npad)*bdim_pad[0];

		    // split into nv-2 triangles
		    for (k = 0; k < nv-2; k++) {
			n3[0][0] = poly[0][0];
			n3[0][1] = poly[0][1];
			n3[1][0] = poly[k+1][0];
			n3[1][1] = poly[k+1][1];
			n3[2][0] = poly[k+2][0];
			n3[2][1] = poly[k+2][1];
			t3.Initialise(n3);

			for (m = 0; m < np; m++) {
			    djac = t3.DetJ(absc[m], &n3);
			    Point glob = t3.Global(n3, absc[m]);
			    Point loc = pel->Local(meshptr->nlist, glob);
			    fun = pel->LocalShapeF (loc);
			    v = wght[m] * djac * Value_nomask(glob,idx_j,false);
			    for (ii = 0; ii < fun.Dim(); ii++) {
				idx_i = pel->Node[ii];
				(*buv)(idx_i, idx_j) += v*fun[ii];
			    }
			}
		    }
		}    
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    delete []rowptr;
    delete []colidx;
    delete []npx;

    buv->Shrink();
    return buv;
}

// ===========================================================================

RCompRowMatrix *Raster_Blob2::CreateDvv_tri () const
{
    const double eps = 1e-10;
    const int subgrd = 200;
    int i, i0, j0, i1, j1, bidx0, bidx1, prog;
    double xcnt0, ycnt0, xcnt1, ycnt1, dstx, dsty;
    double xmin, ymin, xmax, ymax, x, y, r, r0, r1, dvvi;
    double dx0_2, dy0_2, dx1_2, dy1_2;
    RVector g0(2), g1(2);
    double xstep = 1.0/((bdim[0]-1.0)*subgrd);
    double ystep = 1.0/((bdim[1]-1.0)*subgrd);
    double dx = grid[0];
    double dy = grid[1];
    double radlimit2 = sup*sup;
    double radlimit22 = radlimit2*4.0;
    IVector npx(blen_pad);
    struct OVERLAP {
	double dst;
	double v;
    } *overlap;
    int noverlap = 0;
    
    std::cout << "Dvv: pass 1 [" << std::flush; prog = 0;
    for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	ycnt0 = bbmin[1]+j0*dy;
	for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
	    xcnt0 = bbmin[0]+i0*dx;
	    bidx0 = (i0+npad) + (j0+npad)*bdim_pad[0];
	    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
		ycnt1 = bbmin[1]+j1*dy;
		for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
		    xcnt1 = bbmin[0]+i1*dx;
		    dstx = xcnt0-xcnt1;
		    dsty = ycnt0-ycnt1;
		    if (dstx*dstx + dsty*dsty < radlimit22)
			npx[bidx0]++;
		}
	    }
	}
	while (((j0+1)*nprog)/bdim_pad[1] > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    int *rowptr = new int[blen_pad+1];
    rowptr[0] = 0;
    for (i = 0; i < blen_pad; i++)
	rowptr[i+1] = rowptr[i] + npx[i];
    int nz = rowptr[blen_pad];
    int *colidx = new int[nz];
    npx.Clear();

    std::cout << "Dvv: pass 2 [" << std::flush; prog = 0;
    for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	ycnt0 = bbmin[1]+j0*dy;
	for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
	    xcnt0 = bbmin[0]+i0*dx;
	    bidx0 = (i0+npad) + (j0+npad)*bdim_pad[0];
	    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
		ycnt1 = bbmin[1]+j1*dy;
		for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
		    xcnt1 = bbmin[0]+i1*dx;
		    bidx1 = (i1+npad) + (j1+npad)*bdim_pad[0];
		    dstx = xcnt0-xcnt1;
		    dsty = ycnt0-ycnt1;
		    if (dstx*dstx + dsty*dsty < radlimit22)
			colidx[rowptr[bidx0]+npx[bidx0]++] = bidx1;
		}
	    }
	}
	while (((j0+1)*nprog)/bdim_pad[1] > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    RCompRowMatrix *dvv = new RCompRowMatrix (blen_pad, blen_pad,
					      rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    std::cout << "Dvv: pass 3 [" << std::flush; prog = 0;
    for (j0 = -npad; j0 < bdim[1]+npad; j0++) {
	ycnt0 = bbmin[1]+j0*dy;
	for (i0 = -npad; i0 < bdim[0]+npad; i0++) {
	    xcnt0 = bbmin[0]+i0*dx;
	    bidx0 = (i0+npad) + (j0+npad)*bdim_pad[0];
	    for (j1 = -npad; j1 < bdim[1]+npad; j1++) {
		ycnt1 = bbmin[1]+j1*dy;
		for (i1 = -npad; i1 < bdim[0]+npad; i1++) {
		    xcnt1 = bbmin[0]+i1*dx;
		    bidx1 = (i1+npad) + (j1+npad)*bdim_pad[0];
		    dstx = xcnt0-xcnt1;
		    dsty = ycnt0-ycnt1;
		    if (dstx*dstx + dsty*dsty < radlimit22) {
			r = hypot(dstx,dsty);
			// check if we already have an entry
			for (i = 0; i < noverlap; i++) {
			    if (fabs(overlap[i].dst-r) < eps) {
				dvvi = overlap[i].v;
				break;
			    }
			}
			if (i == noverlap) { // not present - calculate
			    xmin = max(xcnt0,xcnt1)-sup;
			    xmax = min(xcnt0,xcnt1)+sup;
			    ymin = max(ycnt0,ycnt1)-sup;
			    ymax = min(ycnt0,ycnt1)+sup;
			    dvvi = 0.0;
			    for (x = xmin; x <= xmax; x += xstep) {
				dstx = x-xcnt0;
				dx0_2 = dstx*dstx;
				dstx = x-xcnt1;
				dx1_2 = dstx*dstx;
				for (y = ymin; y <= ymax; y += ystep) {
				    dsty = y-ycnt0;
				    dy0_2 = dsty*dsty;
				    if (dx0_2 + dy0_2 >= radlimit2) continue;
				    dsty = y-ycnt1;
				    dy1_2 = dsty*dsty;
				    if (dx1_2 + dy1_2 >= radlimit2) continue;
				    r0 = sqrt(dx0_2+dy0_2);
				    r1 = sqrt(dx1_2+dy1_2);
				    g0 = RadGradient(r0);
				    g1 = RadGradient(r1);
				    dvvi += g0[0]*g1[0] + g0[1]*g1[1];
				}
			    }			    
			    OVERLAP *tmp = new OVERLAP[noverlap+1];
			    if (noverlap) {
				memcpy(tmp,overlap,noverlap*sizeof(OVERLAP));
				delete []overlap;
			    }
			    overlap = tmp;
			    overlap[noverlap].dst = r;
			    overlap[noverlap].v = dvvi;
			    noverlap++;
			}
			(*dvv)(bidx0,bidx1) = dvvi;
		    }
		}
	    }
	}
	while (((j0+1)*nprog)/bdim_pad[1] > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;
    (*dvv) *= xstep*ystep;
    
    // diagonal conditioning
    RVector d = dvv->Diag();
    double dmean = mean(d);
    for (i = 0; i < blen; i++)
	(*dvv)(i,i) += dmean*dgscale;

    dvv->Shrink();
    return dvv;
}

// ===========================================================================

RCompRowMatrix *Raster_Blob2::CreateDuv_tri () const
{
    int i, j, k, ii, nv, m, el, nd, idx_i, idx_j, prog;
    int nel = meshptr->elen(), n = meshptr->nlen();
    bool intersect;
    double rx, ry, djac, v;
    RDenseMatrix der, IJ(2,2);
    RVector g;
    
    // grid range and spacing (assumes regular arrangement of blob
    // basis functions)
    double xrange = bbmax[0]-bbmin[0];
    double yrange = bbmax[1]-bbmin[1];
    double dx = xrange/(bdim[0]-1.0);
    double dy = yrange/(bdim[1]-1.0);
    double radlimit2 = sup*sup;
    
    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    int np = QRule_tri_4_6 (&wght, &absc);

    int *npx = new int[n];
    for (i = 0; i < n; i++) npx[i] = 0;

    // pass 1: determine matrix fill structure
    double px, py;
    std::cout << "Duv: pass 1 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j = -npad; j < bdim[1]+npad; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = -npad; i < bdim[0]+npad; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    for (k = 0; k < pel->nNode(); k++)
			npx[pel->Node[k]]++;
		}
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    int *rowptr = new int[n+1];
    rowptr[0] = 0;
    for (i = 0; i < n; i++)
	rowptr[i+1] = rowptr[i]+npx[i];
    int nz = rowptr[n];
    int *colidx = new int[nz];
    for (i = 0; i < n; i++)
	npx[i] = 0;

    // pass 2
    std::cout << "Duv: pass 2 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	for (j = -npad; j < bdim[1]+npad; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = -npad; i < bdim[0]+npad; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    for (k = 0; k < pel->nNode(); k++) {
			nd = pel->Node[k];
			colidx[rowptr[nd]+npx[nd]++] =
			    (i+npad) + (j+npad)*bdim_pad[0];
		    }
		}
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    RCompRowMatrix *duv = new RCompRowMatrix (n, blen_pad, rowptr, colidx);

    // pass 3: fill the matrix
    std::cout << "Duv: pass 3 [" << std::flush; prog = 0;
    for (el = 0; el < nel; el++) {
	Element *pel = meshptr->elist[el];
	const int npoly = 10;
	Point poly[npoly];
	Triangle3 t3;
	NodeList n3(3);
	for (i = 0; i < 3; i++) {
	    t3.Node[i] = i;
	    n3[i].New(2);
	}

	for (j = -npad; j < bdim[1]+npad; j++) {
	    py = bbmin[1] + j*dy;
	    for (i = -npad; i < bdim[0]+npad; i++) {
		px = bbmin[0] + i*dx;
		intersect = false;
		for (k = 0; k < pel->nNode(); k++) {
		    nd = pel->Node[k];
		    rx = px-meshptr->nlist[nd][0];
		    ry = py-meshptr->nlist[nd][1];
		    if (rx*rx + ry*ry < radlimit2) {
			intersect = true;
			break;
		    }
		}
		if (intersect) {
		    nv = SutherlandHodgman (el, px, py, poly, npoly);
		    idx_j = (i+npad) + (j+npad)*bdim_pad[0];

		    // split into nv-2 triangles
		    for (k = 0; k < nv-2; k++) {
			n3[0][0] = poly[0][0];
			n3[0][1] = poly[0][1];
			n3[1][0] = poly[k+1][0];
			n3[1][1] = poly[k+1][1];
			n3[2][0] = poly[k+2][0];
			n3[2][1] = poly[k+2][1];
			t3.Initialise(n3);

			for (m = 0; m < np; m++) {
			    djac = t3.IJacobian(absc[m], &n3, IJ);
			    Point glob = t3.Global(n3, absc[m]);
			    Point loc = pel->Local(meshptr->nlist, glob);
			    der = IJ*pel->LocalShapeD (loc);
			    g = Gradient_nomask(glob,idx_j,false);
			    g *= wght[m] * djac;
			    for (ii = 0; ii < 3; ii++) {
				idx_i = pel->Node[ii];
				v = g[0]*der(ii,0) + g[1]*der(ii,1);
				(*duv)(idx_i, idx_j) += v;
			    }
			}
		    }
		}    
	    }
	}
	if (((el+1)*nprog)/nel > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    delete []rowptr;
    delete []colidx;
    delete []npx;

    duv->Shrink();
    return duv;
}

// ==========================================================================
// Map directly between basis and a regular voxel image with piecewise
// constant basis functions

#ifdef UNDEF
RCompRowMatrix *Raster_Blob2::CreateBvw_tri (
    const IVector &pxdim) const
{
    int i, j, k, i0, j0, i1, j1, idx_i, idx_j, el, nd, nv, s, m, prog;
    int nel = meshptr->elen();
    int stride1 = bdim[0];
    int pstride1 = pxdim[0];
    int plen = pxdim[0]*pxdim[1];
    int *npx = new int[blen_pad];
    for (i = 0; i < blen_pad; i++) npx[i] = 0;
    double px0, py0, px1, py1, rx, ry, djac;
    double dx = grid[0];
    double dy = grid[1];
    double pdx = bbsize[0]/pxdim[0];
    double pdy = bbsize[1]/pxdim[1];
    double radlimit2 = sup*sup;
    bool intersect, intersect0;

    // quadrature rule for local triangle
    const double *wght;
    const Point *absc;
    int np = QRule_tri_9_19 (&wght, &absc);

    // pass 1: determine storage requirements
    std::cout << "Bvw: pass 1 [" << std::flush; prog = 0;
    for (idx_i = 0; idx_i < blen_pad; idx_i++) {
	j0 = idx_i / bdim_pad[0] - npad;
	i0 = idx_i % bdim_pad[0] - npad;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	intersect = false;
	for (el = 0; el < nel; el++) {
	    Element *pel = meshptr->elist[el];
	    for (s = 0; s < pel->nNode(); s++) {
		nd = pel->Node[s];
		rx = px0 - meshptr->nlist[nd][0];
		ry = py0 - meshptr->nlist[nd][1];
		if (rx*rx + ry*ry < radlimit2) {
		    intersect = true;
		    break;
		}
	    }
	    if (intersect)
		break;
	}
	if (intersect) {
	    for (idx_j = 0; idx_j < plen; idx_j++) {
		j1 = idx_j / pstride1;
		i1 = idx_j % pstride1;
		px1 = bbmin[0] + i1*pdx;
		py1 = bbmin[1] + j1*pdy;
		intersect = false;
		for (j = 0; j < 2; j++) {
		    for (i = 0; i < 2; i++) {
			rx = px0 - (px1+i*pdx);
			ry = py0 - (py1+j*pdy);
			if (rx*rx + ry*ry < radlimit2) {
			    intersect = true;
			    i = j = 2; // break;
			}
		    }
		}
		if (intersect) { // blob intersects pixel
		    npx[idx_i]++;
		    // note: we didn't check if the pixel intersects
		    // an element for performance reasons, so some
		    // unnecessary matrix entries may be allocated
		}
	    }
	}
	if (((idx_i+1)*nprog)/blen_pad > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    int *rowptr = new int[blen_pad+1];
    rowptr[0] = 0;
    for (i = 0; i < blen_pad; i++) {
        rowptr[i+1] = rowptr[i]+npx[i];
	npx[i] = 0;
    }
    int nz = rowptr[blen_pad];
    int *colidx = new int[nz];

    // pass 2: sparsity pattern
    std::cout << "Bvw: pass 2 [" << std::flush; prog = 0;
    for (idx_i = 0; idx_i < blen_pad; idx_i++) {
	j0 = idx_i / bdim_pad[0] - npad;
	i0 = idx_i % bdim_pad[0] - npad;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	intersect = false;
	for (el = 0; el < nel; el++) {
	    Element *pel = meshptr->elist[el];
	    for (s = 0; s < pel->nNode(); s++) {
		nd = pel->Node[s];
		rx = px0 - meshptr->nlist[nd][0];
		ry = py0 - meshptr->nlist[nd][1];
		if (rx*rx + ry*ry < radlimit2) {
		    intersect = true;
		    break;
		}
	    }
	    if (intersect)
		break;
	}
	if (intersect) {
	    for (idx_j = 0; idx_j < plen; idx_j++) {
		j1 = idx_j / pstride1;
		i1 = idx_j % pstride1;
		px1 = bbmin[0] + i1*pdx;
		py1 = bbmin[1] + j1*pdy;
		intersect = false;
		for (j = 0; j < 2; j++) {
		    for (i = 0; i < 2; i++) {
			rx = px0 - (px1+i*pdx);
			ry = py0 - (py1+j*pdy);
			if (rx*rx + ry*ry < radlimit2) {
			    intersect = true;
			    i = j = 2; // break;
			}
		    }
		}
		if (intersect) // blob intersects pixel
		    colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
	    }
	}
	if (((idx_i+1)*nprog)/blen_pad > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    delete []npx;
    RCompRowMatrix *bvw = new RCompRowMatrix (blen_pad, plen, rowptr, colidx);

    const int npoly = 10;
    Point poly[npoly];
    Triangle3 t3;
    NodeList n3(3);
    for (i = 0; i < 3; i++) {
	t3.Node[i] = i;
	n3[i].New(2);
    }

    // pass 3: fill the matrix
    std::cout << "Bvw: pass 3 [" << std::flush; prog = 0;
    for (idx_i = 0; idx_i < blen_pad; idx_i++) {
	j0 = idx_i / bdim_pad[0] - npad;
	i0 = idx_i % bdim_pad[0] - npad;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;

	for (el = 0; el < nel; el++) {
	    intersect = false;
	    Element *pel = meshptr->elist[el];
	    for (s = 0; s < pel->nNode(); s++) {
		nd = pel->Node[s];
		rx = px0 - meshptr->nlist[nd][0];
		ry = py0 - meshptr->nlist[nd][1];
		if (rx*rx + ry*ry < radlimit2) {
		    intersect = true;
		    break;
		}
	    }
	    if (!intersect) continue; // element doesn't intersect blob

	    nv = SutherlandHodgman (el, px0, py0, poly, npoly);
	    for (k = 0; k < nv-2; k++) {
		n3[0][0] = poly[0][0];
		n3[0][1] = poly[0][1];
		n3[1][0] = poly[k+1][0];
		n3[1][1] = poly[k+1][1];
		n3[2][0] = poly[k+2][0];
		n3[2][1] = poly[k+2][1];
		t3.Initialise(n3);

		for (m = 0; m < np; m++) {
		    djac = t3.DetJ(absc[m], &n3);
		    Point glob = t3.Global (n3, absc[m]);
		    px1 = glob[0];
		    py1 = glob[1];
		    i1 = max(0, min(pxdim[0]-1, (int)((px1-bbmin[0])/pdx)));
		    j1 = max(0, min(pxdim[1]-1, (int)((py1-bbmin[1])/pdy)));
		    idx_j = i1 + j1*pstride1;
		    (*bvw)(idx_i,idx_j) += Value_nomask(glob,idx_i,false)
			* djac * wght[m];
		}
	    }
	}
	if (((idx_i+1)*nprog)/blen_pad > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    delete []rowptr;
    delete []colidx;
    bvw->Shrink();
    return bvw;
}
#endif

//#ifdef UNDEF
RCompRowMatrix *Raster_Blob2::CreateBvw_tri (const IVector &pxdim) const
{
    int idx_i, idx_j, i, j, i0, j0, i1, j1, jj, prog;
    //int stride1 = bdim[0];
    int pstride1 = pxdim[0];
    int plen = pxdim[0]*pxdim[1];
    int subdiv = 10;
    double px0, py0, px1, py1, rx, ry, x, y, v;
    double radlimit2 = sup*sup;
    double dx = grid[0];
    double dy = grid[1];
    double pdx = bbsize[0]/pxdim[0];
    double pdy = bbsize[1]/pxdim[1];
    double subwght = (pdx*pdy)/(subdiv*subdiv);
    bool intersect;
    int *npx = new int[blen_pad];
    for (i = 0; i < blen_pad; i++) npx[i] = 0;
    Point p(2);

    // pass 1: determine storage requirements
    std::cout << "Bvw: pass 1 [" << std::flush; prog = 0;
    for (idx_i = 0; idx_i < blen_pad; idx_i++) {
	j0 = idx_i / bdim_pad[0] - npad;
	i0 = idx_i % bdim_pad[0] - npad;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	for (idx_j = 0; idx_j < plen; idx_j++) {
	    j1 = idx_j / pstride1;
	    i1 = idx_j % pstride1;
	    px1 = bbmin[0] + i1*pdx;
	    py1 = bbmin[1] + j1*pdy;
	    intersect = false;
	    for (j = 0; j < 2; j++) {
	        for (i = 0; i < 2; i++) {
		    rx = px0 - (px1+i*pdx);
		    ry = py0 - (py1+j*pdy);
		    if (rx*rx + ry*ry < radlimit2) {
		        intersect = true;
			i = j = 2; // break
		    }
		}
	    }
	    if (intersect)
	        npx[idx_i]++;
	}
	if (((idx_i+1)*nprog)/blen_pad > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    int *rowptr = new int[blen_pad+1];
    rowptr[0] = 0;
    for (i = 0; i < blen_pad; i++) {
        rowptr[i+1] = rowptr[i]+npx[i];
	npx[i] = 0;
    }
    int nz = rowptr[blen_pad];
    int *colidx = new int[nz];

    // pass 2: determine matrix sparsity pattern
    std::cout << "Bvw: pass 2 [" << std::flush; prog = 0;
    for (idx_i = 0; idx_i < blen_pad; idx_i++) {
	j0 = idx_i / bdim_pad[0] - npad;
	i0 = idx_i % bdim_pad[0] - npad;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	for (idx_j = 0; idx_j < plen; idx_j++) {
	    j1 = idx_j / pstride1;
	    i1 = idx_j % pstride1;
	    px1 = bbmin[0] + i1*pdx;
	    py1 = bbmin[1] + j1*pdy;
	    intersect = false;
	    for (j = 0; j < 2; j++) {
	        for (i = 0; i < 2; i++) {
		    rx = px0 - (px1+i*pdx);
		    ry = py0 - (py1+j*pdy);
		    if (rx*rx + ry*ry < radlimit2) {
		        intersect = true;
			i = j = 2; // break
		    }
		}
	    }
	    if (intersect)
		colidx[rowptr[idx_i]+npx[idx_i]++] = idx_j;
	}      
	if (((idx_i+1)*nprog)/blen_pad > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    delete []npx;
    RCompRowMatrix *bvw = new RCompRowMatrix (blen_pad, plen, rowptr, colidx);

    std::cout << "Bvw: pass 3 [" << std::flush; prog = 0;
    for (idx_i = 0; idx_i < blen_pad; idx_i++) {
	j0 = idx_i / bdim_pad[0] - npad;
	i0 = idx_i % bdim_pad[0] - npad;
	px0 = bbmin[0] + i0*dx;
	py0 = bbmin[1] + j0*dy;
	for (jj = rowptr[idx_i]; jj < rowptr[idx_i+1]; jj++) {
	    idx_j = colidx[jj];
	    j1 = idx_j / pstride1;
	    i1 = idx_j % pstride1;
	    px1 = bbmin[0] + i1*pdx;
	    py1 = bbmin[1] + j1*pdy;
	    v = 0.0;
	    for (j = 0; j < subdiv; j++) {
	        y = py1 + (j+0.5)/subdiv*pdy;
		for (i = 0; i < subdiv; i++) {
		    x = px1 + (i+0.5)/subdiv*pdx;

		    rx = x-px0;
		    ry = y-py0;
		    if (rx*rx + ry*ry < radlimit2) {
		        p[0] = x; p[1] = y;
			v += Value_nomask(p,idx_i,false);
		    }
		}
	    }
	    (*bvw)(idx_i,idx_j) = v*subwght;
	}
	if (((idx_i+1)*nprog)/blen_pad > prog) {
	    std::cout << "=" << std::flush;
	    prog++;
	}
    }
    std::cout << "]" << std::endl;

    delete []rowptr;
    delete []colidx;

    bvw->Shrink();
    return bvw;
}
//#endif

// ==========================================================================
// Mapping W->V for blob basis V and bilinear pixel basis W
RCompRowMatrix *Raster_Blob2::CreateBvw_linw_tri (const IVector &pxdim) const
{
    const int npad = 1; // generalise!
    const int subgrd = 10;
    
    int i, j, k, m, ip, jp, ib, jb, ii, jj, bidx, pidx;
    double x, y, x0, x1, y0, y1, xcnt, ycnt, dx, dy, r, phi;
    double v, w0, w1, w2, w3;
    bool overlap;
    int plen = prod(pxdim);
    RVector ssize(dim);
    for (i = 0; i < dim; i++) ssize[i] = bbsize[i] / (pxdim[i]*subgrd);
    double As = prod(ssize);
    RVector bvw_px(4), gbvw_px(4);
    
    double pdx = 1.0/(double)(pxdim[0]-1);
    double pdy = 1.0/(double)(pxdim[1]-1);
    double bdx = 1.0/(double)(bdim[0]-1);
    double bdy = 1.0/(double)(bdim[1]-1);
    double sup2 = sup*sup, sup22 = sup*sup*4.0;
    double isubgrd = 1.0/subgrd;
    IVector npx(blen_pad);
    
    // pass 1: matrix structure
    for (jp = 0; jp < pxdim[1]-1; jp++) {
	y0 = jp*pdy;
	y1 = y0+pdy;
	for (ip = 0; ip < pxdim[0]-1; ip++) {
	    x0 = ip*pdx;
	    x1 = x0+pdx;
	    for (jb = -npad; jb < bdim[1]+npad; jb++) {
		ycnt = jb*bdy;
		for (ib = -npad; ib < bdim[0]+npad; ib++) {
		    xcnt = ib*bdx;
		    bidx = (jb+npad)*bdim_pad[0] + (ib+npad);
		    overlap = false;
		    for (k = 0; k < 4; k++) {
			dx = xcnt - (k%2 ? x1:x0);
			dy = ycnt - (k/2 ? y1:y0);
			if (dx*dx + dy*dy < sup2) {
			    overlap = true;
			    break;
			}
		    }
		    if (overlap) npx[bidx] += 4;
		}
	    }
	}
    }
    int *rowptr0 = new int[blen_pad+1];
    rowptr0[0] = 0;
    for (i = 0; i < blen_pad; i++)
	rowptr0[i+1] = rowptr0[i] + npx[i];
    int nz = rowptr0[blen_pad];
    int *colidx0 = new int[nz];
    npx.Clear();
    
    // pass 2: matrix structure
    for (jp = 0; jp < pxdim[1]-1; jp++) {
	y0 = jp*pdy;
	y1 = y0+pdy;
	for (ip = 0; ip < pxdim[0]-1; ip++) {
	    x0 = ip*pdx;
	    x1 = x0+pdx;
	    for (jb = -npad; jb < bdim[1]+npad; jb++) {
		ycnt = jb*bdy;
		for (ib = -npad; ib < bdim[0]+npad; ib++) {
		    xcnt = ib*bdx;
		    bidx = (jb+npad)*bdim_pad[0] + (ib+npad);
		    overlap = false;
		    for (k = 0; k < 4; k++) {
			dx = xcnt - (k%2 ? x1:x0);
			dy = ycnt - (k/2 ? y1:y0);
			if (dx*dx + dy*dy < sup2) {
			    overlap = true;
			    break;
			}
		    }
		    if (overlap) {
			for (k = 0; k < 4; k++) {
			    pidx = (j+k/2)*pxdim[0] + i + k%2;
			    for (m = 0; m < npx[bidx]; m++)
				if (colidx0[rowptr0[bidx]+m] == pidx)
				    break;
			    if (m == npx[bidx])
				colidx0[rowptr0[bidx] + npx[bidx]++] = pidx;
			}
		    }
		}
	    }
	}
    }
    // cut off unused entries
    int *rowptr = new int[blen_pad+1];
    rowptr[0] = 0;
    for (i = 0; i < blen_pad; i++)
	rowptr[i+1] = rowptr[i] + npx[i];
    nz = rowptr[blen_pad];
    int *colidx = new int[nz];
    for (i = 0; i < blen_pad; i++)
	for (j = 0; j < npx[i]; j++)
	    colidx[rowptr[i]+j] = colidx0[rowptr0[i]+j];
    delete []rowptr0;
    delete []colidx0;
    RCompRowMatrix *bvw = new RCompRowMatrix(blen_pad,plen,rowptr,colidx);
    RVector bcount(blen_pad);

    // pass 3: fill matrix
    for (jp = 0; jp < pxdim[1]-1; jp++) {
	y0 = jp*pdy;
	y1 = y0+pdy;
	for (ip = 0; ip < pxdim[0]-1; ip++) {
	    x0 = ip*pdx;
	    x1 = x0+pdx;
	    for (jb = -npad; jb < bdim[1]+npad; jb++) {
		ycnt = jb*bdy;
		for (ib = -npad; ib < bdim[0]+npad; ib++) {
		    xcnt = ib*bdx;
		    bidx = (jb+npad)*bdim_pad[0] + (ib+npad);
		    overlap = false;
		    for (k = 0; k < 4; k++) {
			dx = xcnt - (k%2 ? x1:x0);
			dy = ycnt - (k/2 ? y1:y0);
			if (dx*dx + dy*dy < sup2) {
			    overlap = true;
			    break;
			}
		    }
		    bvw_px.Clear();
		    gbvw_px.Clear();
		    if (overlap) {
			for (jj = 0; jj < subgrd; jj++) {
			    y = y0 + (y1-y0)*(jj+0.5)*isubgrd;
			    for (ii = 0; ii < subgrd; ii++) {
				x = x0 + (x1-x0)*(ii+0.5)*isubgrd;
				dx = x-xcnt;
				dy = y-ycnt;
				if (dx*dx + dy*dy < sup2) {
				    phi = atan2(dy,dx);
				    r = hypot(dx,dy);
				    // generalise for any blob profile!
				    // currently: spline blob
				    if (r < 0.5*sup) {
					v = 2.0/3.0*pow(sup-r,3.0) -
					    8.0/3.0*pow(sup*0.5-r,3.0);
					//dv = 2.0*r*(3.0*r-2.0*sup);
					//dvx = cos(phi)*dv;
					//dvy = sin(phi)*dv;
				    } else {
					v = 2.0/3.0*pow(sup-r,3.0);
					//dv = -2.0*(r-sup)*(r-sup);
					//dvx = cos(phi)*dv;
					//dvy = sin(phi)*dv;
				    }
				    bcount[bidx] += v;
				    w0 = As * (x1-x)*(y1-y);
				    w1 = As * (x-x0)*(y1-y);
				    w2 = As * (x1-x)*(y-y0);
				    w3 = As * (x-x0)*(y-y0);
				    //dw0x = As * (y-y1);
				    //dw0y = As * (x-x1);
				    //dw1x = As * (y1-y);
				    //dw1y = As * (x0-x);
				    //dw2x = As * (y0-y);
				    //dw2y = As * (x1-x);
				    //dw3x = As * (y-y0);
				    //dw3y = As * (x-x0);
				    bvw_px[0] += w0*v;
				    bvw_px[1] += w1*v;
				    bvw_px[2] += w2*v;
				    bvw_px[3] += w3*v;
				    //gbvw_px[0] += dw0x*dvx + dw0y*dvy;
				    //gbvw_px[1] += dw1x*dvx + dw1y*dvy;
				    //gbvw_px[2] += dw2x*dvx + dw2y*dvy;
				    //gbvw_px[3] += dw3x*dvx + dw3y*dvy;
				}
			    }
			}
			for (k = 0; k < 4; k++) {
			    pidx = (j+k/2)*pxdim[0] + i + k%2;
			    (*bvw)(bidx,pidx) += bvw_px[k];
			    //gbvw(bidx,pidx) += gbvw_px[k];
			}
		    }	    
		}
	    }
	}
    }
    double *pbvw = bvw->ValPtr();
    for (i = 0; i < blen_pad; i++) {
	if (bcount[i]) {
	    double scale = 1.0/sqrt(bcount[i]);
	    for (j = rowptr[i]; j < rowptr[i+1]; j++)
		pbvw[j] *= scale;
	}
    }

    delete []rowptr;
    delete []colidx;

    return bvw;
}


int Raster_Blob2::SutherlandHodgman (const Point *p0, int np0,
				     const Point *p1, int np1,
				     Point *clip_poly, const int npoly) const
{
    const int maxpoly = 20;
    int i;

    static vec_t v0[maxpoly];
    for (i = 0; i < np0; i++) {
	v0[i].x = p0[i][0];
	v0[i].y = p0[i][1];
    }
    poly_t subject = {np0, 0, v0};

    static vec_t v1[maxpoly];
    for (i = 0; i < np1; i++) {
	v1[i].x = p1[i][0];
	v1[i].y = p1[i][1];
    }
    poly_t clipper = {np1, 0, v1};

    poly res = poly_clip (&subject, &clipper);
    int len = res->len;
    if (len > npoly) {
	len = -1;
    } else {
	for (i = 0; i < len; i++) {
	    clip_poly[i].New(2);
	    clip_poly[i][0] = res->v[i].x;
	    clip_poly[i][1] = res->v[i].y;
	}
    }
    poly_free(res);
    return len;
}



int Raster_Blob2::SutherlandHodgman (int el, double px, double py,
    Point *clip_poly, int npoly) const
{
    int i;

    Element *pel = meshptr->elist[el];
    int ni, nv = pel->nNode();
    if (nv > npoly) return -1;

    vec_t tri[3];
    for (i = 0; i < 3; i++) {
	tri[i].x = meshptr->nlist[pel->Node[i]][0];
	tri[i].y = meshptr->nlist[pel->Node[i]][1];
    }
    poly_t subject = {3, 0, tri};

    const int nseg = 32; // segments for approximating support circle
    vec_t clip[nseg];
    for (i = 0; i < nseg; i++) {
	double phi = (double)i/(double)nseg * Pi2;
	clip[i].x = cos(phi)*sup + px;
	clip[i].y = sin(phi)*sup + py;
    }
    poly_t clipper = {nseg, 0, clip};

    poly res = poly_clip (&subject, &clipper);
    int len = res->len;
    if (len > npoly) {
	len = -1;
    } else {
	for (i = 0; i < len; i++) {
	    clip_poly[i].New(2);
	    clip_poly[i][0] = res->v[i].x;
	    clip_poly[i][1] = res->v[i].y;
	}
    }
    poly_free(res);
    return len;
}
