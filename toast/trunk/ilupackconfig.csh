#should be called after setting $TOASTDIR
setenv CC /opt/gcc-4.3.2/bin/gcc
setenv CXX /opt/gcc-4.3.2/bin/g++
setenv F77 /opt/gcc-4.3.2/bin/gfortran
./configure --with-ilu='$(TOASTDIR)/numerics/ILUPACK'
setenv LD_LIBRARY_PATH /usr/opt/gcc-4.3.2/lib64:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

