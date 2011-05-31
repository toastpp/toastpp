#should be source'd after setting $TOASTDIR
setenv CC /opt/gcc-4.3.2/bin/gcc
setenv CXX /opt/gcc-4.3.2/bin/g++
setenv F77 /opt/gcc-4.3.2/bin/gfortran
setenv LD_RUN_PATH /opt/gcc-4.3.2/lib64/
./configure --with-ilu='$(TOASTDIR)/numerics/ILUPACK' --with-superlu='-lsuperlu'
setenv LD_LIBRARY_PATH /usr/opt/gcc-4.3.2/lib64:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

