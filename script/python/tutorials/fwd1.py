# This pytoast example solves the forward problem
# for a homogeneous 2D circular problem

import pdb

# Import various modules
import os
import sys
import numpy as np
from numpy import matrix
from scipy import sparse
from scipy.sparse import linalg
from numpy.random import rand
import matplotlib.pyplot as plt

# PyToast environment
execfile(os.getenv("TOASTDIR") + "/ptoast_install.py")

# Import the toast modules
from toast import mesh

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

# Homogeneous parameter distributions
mua = np.ones ((1,nlen)) * 0.025
mus = np.ones ((1,nlen)) * 2.0
ref = np.ones ((1,nlen)) * 1.4
freq = 100

# Set up the linear system
smat = mesh.Sysmat (hmesh, mua, mus, ref, freq)
qvec = mesh.Qvec (hmesh)
qvec = qvec.transpose()
mvec = mesh.Mvec (hmesh)

# Solve the linear system
nq = qvec.shape[1]
phi = np.zeros((nlen,nq),dtype=np.cdouble)
for q in range(nq):
    qq = qvec[:,q].todense()
    res = linalg.bicgstab(smat,qq,tol=1e-12)
    phi[:,q] = res[0]

# Project to boundary
y = mvec * phi
logy = np.log(y)

# Display as sinogram
plt.figure(1)
im = plt.imshow(logy.real,interpolation='none')
plt.title('log amplitude')
plt.xlabel('detector index')
plt.ylabel('source index')
plt.colorbar()
plt.draw()
#plt.show()

plt.figure(2)
im = plt.imshow(logy.imag,interpolation='none')
plt.title('phase')
plt.xlabel('detector index')
plt.ylabel('source index')
plt.colorbar()
plt.show()

