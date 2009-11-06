TOAST- Image Reconstruction in Optical Tomography
2008-2009 (c) Martin Schweiger and Simon Arridge, University College London

========================================================================

TOAST is a software suite for image reconstruction in optical tomography.
It is being developed by Martin Schweiger and Simon Arridge at University
College London.

The TOAST suite is available for download online at

    http://web4.cs.ucl.ac.uk/research/vis/toast/

Please check the TOAST web site regularly for updates.

The TOAST License is described in the license.html document.


Downloading TOAST:
------------------
TOAST is available for different computer platforms:

- Windows
- Linux (32-bit)
- Linux (64-bit)

You need to download two archive files from the TOAST download page:

- The "Common" package containing platform-independent scripts and text
  files
- The platform-specific package for your computer, containing libraries
  and mex files.


Installing TOAST:
-----------------
Unzip or untar the two files into a common directory. This will create
a sub-directory "toast2009" containing the complete TOAST suite.

For Linux systems, a TOASTDIR environment variable must be defined pointing
to the TOAST root directory. The additional TOAST environment settings can
then be defined by executing the toastenv scripts.

Example (csh):

  setenv TOASTDIR $HOME/toast2009
  source $TOASTDIR/toastenv.csh

Example (bash):

  export TOASTDIR=$HOME/toast2009
  source $TOASTDIR/toastenv.sh

For all computer systems, the Matlab environment is set by running the
mtoast_install script:

- Start Matlab
- Change the working directory to your toast2009 root directory.
- Run mtoast_install.

This will add the TOAST mex and script directories to your Matlab search
path. You should save the new path definition to make it permanent.


Getting started:
----------------
After setting up the TOAST toolbox, you can start by running the included
demos. In Matlab, type

  demo toolbox toast

This will bring up a list of available TOAST demos, showing some of the
capabilities of the TOAST suite.

To view the help pages for the TOAST toolbox, type

  doc toast

This contains a brief introduction, a list of available functions, and links
to the TOAST web page.

For troubleshooting, please refer to the TOAST Frequently Asked Question
page and the message board, accessible from the TOAST home page. If you
cannot resolve a problem, you can contact the authers via the Contacts page.
