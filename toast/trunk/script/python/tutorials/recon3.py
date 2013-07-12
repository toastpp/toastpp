# This example builds on recon2.py by adding a
# regularisation term.
#
# Note: run this with
#
#     ipython -pylab recon3.py
#
# to avoid python blocking on opening the figure


import pdb

# Import various modules
import os
import math
import numpy as np
from numpy import matrix
from scipy import sparse
from scipy.sparse import linalg
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.ion()

itrmax = 100  # max number of nonlinear iterations
tolCG = 1e-7
resetCG = 10
grd = np.array([100,100])
noiselevel = 0.01
tau = 1e-3
beta = 0.01

# ---------------------------------------------------
# Objective function
def objective(proj,data,sd,logx):
    err_data = np.sum(np.power((data-proj)/sd,2))
    err_prior = regul.Value (hreg, logx)
    return err_data + err_prior


# ---------------------------------------------------
# Objective function for line search callback
def objective_ls(logx):
    x = np.exp(logx)
    slen = x.shape[0]/2
    scmua = x[0:slen]
    sckap = x[slen:2*slen]
    smua = scmua/cm
    skap = sckap/cm
    smus = 1/(3*skap) - smua
    mua = raster.MapBasis(hraster, 'S->M', smua)
    mus = raster.MapBasis(hraster, 'S->M', smus)
    phi = mesh.Fields(hmesh,-1,qvec,mvec,mua,mus,ref,freq,'d')
    p = projection(phi,mvec)
    return objective(p,data,sd,logx)


# ---------------------------------------------------
# Projections from fields
def projection(phi,mvec):
    gamma = mvec.transpose() * phi
    gamma = np.reshape(gamma,(-1,1),'F')
    lgamma = np.log(gamma)
    lnamp = lgamma.real
    phase = lgamma.imag
    return np.concatenate((lnamp,phase))


# ---------------------------------------------------
# Image error
def imerr(im1,im2):
    im1 = np.reshape(im1,-1,1)
    im2 = np.reshape(im2,-1,1)
    err = np.sum(np.power(im1-im2,2))/np.sum(np.power(im2,2))
    return err


# PyToast environment
execfile(os.getenv("TOASTDIR") + "/ptoast_install.py")

# Import the toast modules
from toast import mesh
from toast import raster
from toast import regul

# Set the file paths
meshdir = os.path.expandvars("$TOASTDIR/test/2D/meshes/")
meshfile1 = meshdir + "ellips_tri10.msh"  # mesh for target data generation
meshfile2 = meshdir + "circle25_32.msh"   # mesh for reconstruction
qmfile = meshdir + "circle25_32x32.qm"    # source-detector file
muafile = meshdir + "tgt_mua_ellips_tri10.nim" # nodal target absorption
musfile = meshdir + "tgt_mus_ellips_tri10.nim" # nodal target scattering

# A few general parameters
c0 = 0.3        # speed of light in vacuum [mm/ps]
refind = 1.4    # refractive index in medium (homogeneous)
cm = c0/refind; # speed of light in medium
freq = 100      # modulation frequency [MHz]

qtype  = 'Neumann'      # source type
qprof  = 'Gaussian'     # source profile
qwidth = 2              # source width
mprof  = 'Gaussian'     # detector profile
mwidth = 2              # detector width

# ---------------------------------------------------
# Generate target data

# Set up mesh geometry
hmesh_fwd = mesh.Read(meshfile1)
mesh.ReadQM(hmesh_fwd,qmfile)
qvec = mesh.Qvec (hmesh_fwd,type=qtype,shape=qprof,width=qwidth)
mvec = mesh.Mvec (hmesh_fwd,shape=mprof,width=mwidth)
nlen = mesh.NodeCount (hmesh_fwd)
nqm = qvec.shape[0] * mvec.shape[0]
ndat = nqm*2

# Target parameters
mua = mesh.ReadNim (muafile)
mus = mesh.ReadNim (musfile)
ref = np.ones((1,nlen)) * refind

# Parameter plotting ranges
mua_min = 0.015 # np.min(mua)
mua_max = 0.055 # np.max(mua)
mus_min = 1     # np.min(mus)
mus_max = 4.5   # np.max(mus)

# Solve forward problem
phi = mesh.Fields(hmesh_fwd,-1,qvec,mvec,mua,mus,ref,freq)
data = projection(phi,mvec)

# Add noise
data = data + data*noiselevel*np.random.normal(0,1,data.shape)

lnamp_tgt = data[0:nqm]
phase_tgt = data[nqm:nqm*2]

# Map target parameters to images for display
hraster_fwd = raster.Make (hmesh_fwd,grd)
bmua_tgt = np.reshape(raster.MapBasis (hraster_fwd, 'M->B', mua),grd)
bmus_tgt = np.reshape(raster.MapBasis (hraster_fwd, 'M->B', mus),grd)

# Clear objects
raster.Clear (hraster_fwd)
mesh.Clear (hmesh_fwd)


# ---------------------------------------------------
# Solve inverse problem

# Set up mesh geometry
hmesh = mesh.Read(meshfile2)
mesh.ReadQM(hmesh,qmfile)
qvec = mesh.Qvec (hmesh,type='Neumann',shape='Gaussian',width=2)
mvec = mesh.Mvec (hmesh,shape='Gaussian',width=2)
nlen = mesh.NodeCount (hmesh)

# Initial parameter estimates
mua = np.ones(nlen) * 0.025
mus = np.ones(nlen) * 2
kap = 1/(3*(mua+mus))
ref = np.ones(nlen) * refind

# Solution basis
hraster = raster.Make (hmesh,grd);

# Initial projections
phi = mesh.Fields(hmesh,-1,qvec,mvec,mua,mus,ref,freq)
proj = projection(phi,mvec)
lnamp = proj[0:nqm]
phase = proj[nqm:nqm*2]

# Data scaling
sd_lnamp = np.ones(lnamp.shape) * np.linalg.norm(lnamp_tgt-lnamp)
sd_phase = np.ones(phase.shape) * np.linalg.norm(phase_tgt-phase)
sd = np.concatenate((sd_lnamp,sd_phase))

# Map parameter estimates to solution basis
bmua = raster.MapBasis(hraster,'M->B',mua)
bmus = raster.MapBasis(hraster,'M->B',mus)
bkap = raster.MapBasis(hraster,'M->B',kap)
bcmua = bmua * cm
bckap = bkap * cm
scmua = raster.MapBasis(hraster,'B->S',bcmua)
sckap = raster.MapBasis(hraster,'B->S',bckap)

# Vector of unknowns
x = np.asmatrix(np.concatenate((scmua,sckap))).transpose()
logx = np.log(x)
slen = x.shape[0]/2

# Create regularisation object
#pdb.set_trace()
#hreg = regul.Make ("TK1", hraster, logx, tau);
hreg = regul.Make ("TV", hraster, logx, tau, beta=beta);

# Initial error
err0 = objective(proj,data,sd,logx)
err = err0
errp = 1e10
erri = np.array([err])
errmua = np.array([imerr(bmua,bmua_tgt)])
errmus = np.array([imerr(bmus,bmus_tgt)])

itr = 1
step = 1.0

hfig1=plt.figure(1)
plt.show()
plt.figure(2)
plt.show()

while itr <= itrmax and err > tolCG*err0 and errp-err > tolCG:
    errp = err
    r = -mesh.Gradient (hmesh,hraster,qvec,mvec,mua,mus,ref,freq,data,sd)
    for i in range(r.shape[0]):
        r[i] = r[i] * np.asscalar(x[i]) # parameter scaling

    rr = -regul.Gradient (hreg, logx)

    r = r + rr
    r = np.transpose(np.matrix(r))
    
    plt.figure(2)
    plt.clf()
    plt.subplot(2,2,1)
    im = plt.imshow (np.reshape (raster.MapBasis (hraster, 'S->B', rr[0:slen]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mua prior gradient")
    plt.colorbar()

    plt.subplot(2,2,2)
    im = plt.imshow (np.reshape (raster.MapBasis (hraster, 'S->B', rr[slen:slen*2]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("kap prior gradient")
    plt.colorbar()

    plt.subplot(2,2,3)
    im = plt.imshow (np.reshape (raster.MapBasis (hraster, 'S->B', r[0:slen]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mua tot gradient")
    plt.colorbar()

    plt.subplot(2,2,4)
    im = plt.imshow (np.reshape (raster.MapBasis (hraster, 'S->B', r[slen:slen*2]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("kap tot gradient")
    plt.colorbar()

    plt.draw()
    
    if itr > 1:
        delta_old = delta_new
        delta_mid = np.dot (r.transpose(),s)
        
    s = r # replace this with preconditioner

    if itr == 1:
        d = s
        delta_new = np.dot(r.transpose(),d)
        delta0 = delta_new
    else:
        delta_new = np.dot(r.transpose(),s)
        beta = (delta_new-delta_mid) / delta_old
        if itr % resetCG == 0 or beta <= 0:
            d = s
        else:
            d = s + d*beta

    delta_d = np.dot(d.transpose(),d)
    step,err = mesh.Linesearch (logx, d, step, err, objective_ls)

    logx = logx + d*step
    x = np.exp(logx)
    scmua = x[0:slen]
    sckap = x[slen:2*slen]
    smua = scmua/cm
    skap = sckap/cm
    smus = 1/(3*skap) - smua
    mua = raster.MapBasis(hraster, 'S->M', smua)
    mus = raster.MapBasis(hraster, 'S->M', smus)

    bmua = np.reshape(raster.MapBasis(hraster, 'S->B', smua),grd)
    bmus = np.reshape(raster.MapBasis(hraster, 'S->B', smus),grd)

    erri=np.concatenate((erri,[err]))
    errmua = np.concatenate((errmua,[imerr(bmua,bmua_tgt)]))
    errmus = np.concatenate((errmus,[imerr(bmus,bmus_tgt)]))
    print ("Iteration "+str(itr)+", objective "+str(err))

    plt.figure(1)
    plt.clf()
    hfig1.suptitle("Iteration "+str(itr))

    plt.subplot(2,3,1)
    im = plt.imshow(bmua_tgt, vmin=mua_min, vmax=mua_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mua target")
    plt.colorbar()

    plt.subplot(2,3,2)
    im = plt.imshow(bmus_tgt, vmin=mus_min, vmax=mus_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mus target")
    plt.colorbar()

    plt.subplot(2,3,4)
    im = plt.imshow(bmua, vmin=mua_min, vmax=mua_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mua recon")
    plt.colorbar()
    
    plt.subplot(2,3,5)
    im = plt.imshow(bmus, vmin=mus_min, vmax=mus_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mus recon")
    plt.colorbar()

    plt.subplot(2,3,3)
    im = plt.semilogy(erri)
    plt.title("objective function")
    plt.xlabel("iteration")
    
    plt.subplot(2,3,6)
    im = plt.semilogy(errmua)
    im = plt.semilogy(errmus)
    plt.title("rel. image error")
    plt.xlabel("iteration")
    
    plt.draw()
    
    itr = itr+1

plt.ioff()

