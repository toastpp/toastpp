Toast++- Image Reconstruction in Optical Tomography
2008-2016 (c) Martin Schweiger and Simon Arridge, University College London

========================================================================

Toast++ is an open-source software suite for image reconstruction in optical
tomography.
It is being developed by Martin Schweiger and Simon Arridge at University
College London.

The TOAST suite is available for download online at

    http://toastplusplus.org

Please check the TOAST web site regularly for updates.

Toast++ is distributed with a GPL v3 license. Please see file COPYING for
details.


Downloading Toast++:
--------------------

Source distribution: The Toast++ sources are available via an SVN repository
on sourceforge (https://sourceforge.net/projects/toastpp)

To check out the current version:
svn checkout http://svn.code.sf.net/p/toastpp/code/toast/trunk toastpp


Pre-compiled binaries: Binary packages for several computer platforms are
hosted on sourceforge:
https://sourceforge.net/projects/toastpp/files

You will need to download two zip files:
- The "Common" package containing device-independent scripts and examples
- One of the platform-dependent packages for Windows, Linux or Mac-OS.


Installing precompiled Toast++ packages:
----------------------------------------

If you are using the pre-compiled packages, just unzip both packages into THE
SAME directory. You should end up with a single "toast" root directory.

Note that some unzip utilities, in particular on Windows, store the zip contents
under a directory named after the zip file, e.g.

toast_v2.0.0_common.zip  -->  toast_v2.0.0_common/

In that case, you will end up with two separate directories, and must merge the
two "toast" subdirectories from both manually.


Compiling from sources:
-----------------------

Please follow the instructions under doc/install.


Preparing the Toast++ environment:
----------------------------------

(Required for both precompiled and compiled installations)

On Linux, type

cd <toast root dir>
export TOASTDIR=$PWD
source toastenv.sh

This step can be skipped for Windows installations.

For setting the Toast++ paths for the Matlab toolbox:
From the Matlab prompt:

cd <toast root dir>
mtoast2_install

Save the paths to make them permanent.


Getting started:
----------------
After setting up the TOAST toolbox, you can start by running the included
demos. In Matlab, type

  doc

Follow the link to 

This will bring up a list of available TOAST demos, showing some of the
capabilities of the TOAST suite.

To view the help pages for the TOAST toolbox, type

  doc toast

This contains a brief introduction, a list of available functions, and links
to the TOAST web page.

For troubleshooting, please refer to the TOAST Frequently Asked Question
page and the message board, accessible from the TOAST home page. If you
cannot resolve a problem, you can contact the authers via the Contacts page.
