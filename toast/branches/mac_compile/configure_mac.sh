#Set up for building TOAST for Snow Leopard with gcc4.3.5 from macports
#source configure_mac.sh
export TOASTDIR=$(PWD) 
#Change to your build tools
/opt/local/bin/autoconf
#export CC=/opt/local/bin/gcc-mp-4.3
export CC=gcc
#export LDFLAGS=-L/usr/lib/
#export LIBS=-lstdc++
#export CXX=/opt/local/bin/g++-mp-4.3
export CXX=g++
export F77=/opt/local/bin/gfortran-mp-4.3
#setenv LD_RUN_PATH /opt/gcc-4.3.2/lib64/
#./configure CC=$CC CXX=$CXX F77=$F77 --with-ilu='$(TOASTDIR)/numerics/ILUPACK'
./configure CC=$CC CXX=$CXX F77=$F77 
#LDFLAGS=$LDFLAGS LIBS=$LIBS
source $TOASTDIR/toastenv.sh
touch $TOASTDIR/include/malloc.h