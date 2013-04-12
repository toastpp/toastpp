# This pytoast example solves the forward problem
# for a homogeneous 2D circular problem

# Import various modules
import os
import numpy as np
from numpy import matrix
from scipy import sparse
from scipy.sparse import linalg
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#from matplotlib import ion

plt.ion()

# PyToast environment
execfile(os.getenv("TOASTDIR") + "/ptoast_install.py")

# Import the toast modules
from toast import mesh
from toast import raster

# Set the file paths
meshdir = os.path.expandvars("$TOASTDIR/test/2D/meshes/")
meshfile = meshdir + "circle25_32.msh"
qmfile = meshdir + "circle25_32x32.qm"

# Load the mesh and source/detector specs
hmesh = mesh.Read(meshfile)
mesh.ReadQM(hmesh,qmfile)

# Extract mesh geometry
nlist,elist,eltp = mesh.Data (hmesh)
nlen = nlist.shape[0]

grd = np.array([64,64])
hraster = raster.Make (hmesh,grd);

# Homogeneous parameter distributions
mua = np.ones ((1,nlen)) * 0.025
mus = np.ones ((1,nlen)) * 2.0
ref = np.ones ((1,nlen)) * 1.4
freq = 100

# Set up the linear system
qvec = mesh.Qvec (hmesh)
mvec = mesh.Mvec (hmesh)
nq = qvec.shape[1]

phi = mesh.Fields(hmesh,hraster,qvec,mvec,mua,mus,ref,freq)

fig1 = plt.figure()
ims = []
for q in range(nq):
    lphi = np.log (phi[:,q])
    bphi = raster.MapBasis (hraster, 'S->B', lphi)
    bphi = np.reshape (bphi, grd)
    ims.append ((plt.imshow(bphi.real),))

im_ani = animation.ArtistAnimation(fig1, ims, interval=50, repeat_delay=3000, blit=True)

plt.show()


# Display fields
plt.subplot(1,2,1)
im = plt.imshow(bphi.real)
plt.title('log amplitude')
plt.colorbar()
#plt.draw()
#plt.show()
plt.subplot(1,2,2)
im = plt.imshow(bphi.imag)
plt.title('phase')
plt.colorbar()
plt.show()

