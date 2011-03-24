#include <fstream>
#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

void usage()
{
  cerr << "Usage: pv_materials [out_pv_materials_file] [number_of_pvlabels (int)] [TC_bck TC_csf TC_gm TC_wm TC_hippo] <options>\n";
  cerr << "where <options> is one or more of the following:\n";
  cerr << "\t <-y Young's modulus BCK CSF GM WM hippos (double)> (default = 1.0 1.0 8.0 4.0 6.0)\n\n";
  cerr << "\t <-p Poisson's ratio (double)> (default = 0.45)\n\n";
  cerr << "\t <-d density (double)> (default = 1e-6)\n\n";
  cerr << "\t <-r n_inexistent_labels (int) list_inexistent_labels (int's)> \n\n";
  cerr << endl;
  cerr << "TC (double)= Thermal Coefficient (background, CSF, GM, WM, hippocampi).\n\n" << endl;
  cerr << "Construct a text file with the biomechanical properties (including thermal coefficient) \n" << endl;
  cerr << "corresponding to partial volume materials. \n" << endl;
  cerr << "The r option is used to remove inexistent labels. It needs the number of labels to remove . \n" << endl;
  cerr << "and a list of labels to remove. \n" << endl;
  cerr << endl;
  exit(1);
}

int main (int argc, char **argv) {
  
  char *matfile_out_name = NULL;
  
  //double TC_bck = 0.0, TC_csf = 0.0, TC_gm = 0.0, TC_wm = 0.0, TC_hippos = 0.0;
  int i = 0, j = 0, k = 0, l = 0, m = 0, n = 0;
  double poisson = 0.45, density = 1e-6;
  double TCs[31], y[31];
  int npvlabels = 31, nlabels_remove = 0, nfinal_labels = 31;
  int *labels_remove;
  
  bool ok = false, remove = false;
  
  for ( i = 0; i < 31; i++)
    {
      TCs[i] = 0.0;
      y[i]= 0.0;
    }
  y[0] = 1.0;
  y[1] = 1.0;
  y[2] = 8.0;
  y[3] = 4.0;
  y[4] = 6.0;

  srand(1234567);
  
  if (argc < 8){
    usage();
  }
  matfile_out_name = argv[1];
  argc--;
  argv++;
  npvlabels = atoi(argv[1]);
  argc--;
  argv++;
  //double *TCs = new double[npvlabels];
  TCs[0] = atof(argv[1]);
  argc--;
  argv++;
  TCs[1] = atof(argv[1]);
  argc--;
  argv++;
  TCs[2] = atof(argv[1]);
  argc--;
  argv++;
  TCs[3] = atof(argv[1]);
  argc--;
  argv++;
  TCs[4] = atof(argv[1]);
  argc--;
  argv++;
  
  while (argc > 1) {
    ok = false;
    
    if ((ok == false) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      //young = true; 
      ok = true;
      y[0] = atof(argv[1]);
      argc--;
      argv++;
      y[1] = atof(argv[1]);
      argc--;
      argv++;
      y[2] = atof(argv[1]);
      argc--;
      argv++;
      y[3] = atof(argv[1]);
      argc--;
      argv++;
      y[4] = atof(argv[1]);
      argc--;
      argv++;
      
    }    
    if ((ok == false) && (strcmp(argv[1], "-p") == 0)){
      argc--;
      argv++;
      //young = true; 
      ok = true;
      poisson = atof(argv[1]);
      argc--;
      argv++;
    }    
    if ((ok == false) && (strcmp(argv[1], "-d") == 0)){
      argc--;
      argv++;
      //young = true; 
      ok = true;
      density = atof(argv[1]);
      argc--;
      argv++;
    }    
    if ((ok == false) && (strcmp(argv[1], "-r") == 0)){
      argc--;
      argv++;
      ok = true;
      remove = true;
      nlabels_remove = atoi(argv[1]);
      labels_remove = new int [nlabels_remove];
      argc--;
      argv++;
      for (i = 0; i < nlabels_remove; i++)
	{
	  labels_remove[i] = atoi(argv[1]);
	  argc--;
	  argv++;
	}
    }    
    if ( !ok ) {
      usage();
    }
  }
  
  //computation of thermal coefficients of pv labels (here, only case of 5 original labels, then 29 pv labels
  
  //combination of 2 labels 
  i = 5;
  for ( j = 0; j < 4; j++ )
    {
      for ( k = j+1; k < 5; k++)
	{
	  TCs[i] = (TCs[j] + TCs[k]) / 2.0;
	  y[i] = (y[j] + y[k]) / 2.0;
	  i++;
	}
    }
  //combination of 3 labels 
  for ( j = 0; j < 4; j++ )
    {
      for ( k = j+1; k < 5; k++)
	{
	  for ( l = k+1; l < 5; l++)
	    {
	      TCs[i] = (TCs[j] + TCs[k] + TCs[l]) / 3.0;
	      y[i] = (y[j] + y[k] + y[l]) / 3.0;
	      //cout << "(j,k,l)=(" << j << "," << k << "," << l << ")" << endl;
	      i++;
	    }
	}
    }
  //combination of 4 labels 
  for ( j = 0; j < 4; j++ )
    {
      for ( k = j+1; k < 5; k++)
	{
	  for ( l = k+1; l < 5; l++)
	    {
	      for ( m = l+1; m < 5; m++)
		{
		  TCs[i] = (TCs[j] + TCs[k] + TCs[l] + TCs[m]) / 4.0;
		  y[i] = (y[j] + y[k] + y[l] + y[m]) / 4.0;
		  //cout << "(j,k,l,m)=(" << j << "," << k << "," << l << "," << m <<")"  <<  endl;
		  i++;
		}
	    }
	}
    }
  //combination of 5 labels 
  for ( j = 0; j < 4; j++ )
    {
      for ( k = j+1; k < 5; k++)
	{
	  for ( l = k+1; l < 5; l++)
	    {
	      for ( m = l+1; m < 5; m++)
		{
		  for ( n = m+1; n < 5; n++)
		    {
		      TCs[i] = (TCs[j] + TCs[k] + TCs[l] + TCs[m] + TCs[n]) / 5.0;
		      y[i] = (y[j] + y[k] + y[l] + y[m] + y[n]) / 5.0;
		      //cout << "(j,k,l,m)=(" << j << "," << k << "," << l << "," << m << "," << n <<")" <<  endl;
		      i++;
		    }
		}
	    }
	}
    }

  if (remove)
    {
      for ( i = 0; i < nlabels_remove; i++)
	TCs[labels_remove[i]] = 9999;
      nfinal_labels = nfinal_labels - nlabels_remove;
    }

  //write the pv_materials_file
    ofstream ofs (matfile_out_name);
    ofs << nfinal_labels << endl;
    for ( i = 0; i < npvlabels; i++)
      {
	if ( TCs[i] != 9999 )
	  {
	    ofs << y[i] << " " << poisson << " " << density << " " << TCs[i] << endl; 
	  } 
      }

  return 0;
}                                                                              

