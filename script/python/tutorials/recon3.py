# This example builds on recon2.py by adding a
# regularisation term.
#
# Note: run this with
#
#     ipython -pylab recon3.py
#
# to avoid python blocking on opening the figure


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
grd = np.array([100, 100])
noiselevel = 0.01
tau = 1e-3
beta = 0.01

# ---------------------------------------------------
# Objective function
def objective(proj,data,sd,logx):
    err_data = np.sum(np.power((data-proj)/sd, 2))
    err_prior = regul.Value(logx)
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
    mua = basis_inv.Map('S->M', smua)
    mus = basis_inv.Map('S->M', smus)
    phi = mesh_inv.Fields(None, qvec, mua, mus, ref, freq)
    p = projection(phi, mvec)
    return objective(p, data, sd, logx)


# ---------------------------------------------------
# Projections from fields
def projection(phi, mvec):
    gamma = mvec.transpose() * phi
    gamma = np.reshape(gamma, (-1, 1), 'F')
    lgamma = np.log(gamma)
    lnamp = lgamma.real
    phase = lgamma.imag
    return np.concatenate((lnamp, phase))


# ---------------------------------------------------
# Image error
def imerr(im1, im2):
    im1 = np.reshape(im1, -1, 1)
    im2 = np.reshape(im2, -1, 1)
    err = np.sum(np.power(im1-im2, 2))/np.sum(np.power(im2, 2))
    return err


# PyToast environment
execfile(os.getenv("TOASTDIR") + "/ptoast_install.py")
import toast

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
cm = c0/refind  # speed of light in medium
freq = 100      # modulation frequency [MHz]

qtype  = 'Neumann'      # source type
qprof  = 'Gaussian'     # source profile
qwidth = 2              # source width
mprof  = 'Gaussian'     # detector profile
mwidth = 2              # detector width

# ---------------------------------------------------
# Generate target data

# Set up mesh geometry
mesh_fwd = toast.Mesh(meshfile1)
mesh_fwd.ReadQM(qmfile)
qvec = mesh_fwd.Qvec(type=qtype, shape=qprof, width=qwidth)
mvec = mesh_fwd.Mvec(shape=mprof, width=mwidth, ref=refind)
nlen = mesh_fwd.NodeCount()
nqm = qvec.shape[1] * mvec.shape[1]
ndat = nqm*2

# Target parameters
mua = mesh_fwd.ReadNim(muafile)
mus = mesh_fwd.ReadNim(musfile)
ref = np.ones((1, nlen)) * refind

# Parameter plotting ranges
mua_min = 0.015 # np.min(mua)
mua_max = 0.055 # np.max(mua)
mus_min = 1     # np.min(mus)
mus_max = 4.5   # np.max(mus)

# Solve forward problem
phi = mesh_fwd.Fields(None, qvec, mua, mus, ref, freq)
data = projection(phi, mvec)

# Add noise
data = data + data*noiselevel*np.random.normal(0, 1, data.shape)

lnamp_tgt = data[0:nqm]
phase_tgt = data[nqm:nqm*2]

# Map target parameters to images for display
basis_fwd = toast.Raster(mesh_fwd, grd)
bmua_tgt = np.reshape(basis_fwd.Map('M->B', mua), grd)
bmus_tgt = np.reshape(basis_fwd.Map('M->B', mus), grd)


# ---------------------------------------------------
# Solve inverse problem

# Set up mesh geometry
mesh_inv = toast.Mesh(meshfile2)
mesh_inv.ReadQM(qmfile)
qvec = mesh_inv.Qvec(type=qtype, shape=qprof, width=qwidth)
mvec = mesh_inv.Mvec(shape=mprof, width=mwidth, ref=refind)
nlen = mesh_inv.NodeCount()

# Initial parameter estimates
mua = np.ones(nlen) * 0.025
mus = np.ones(nlen) * 2
kap = 1/(3*(mua+mus))
ref = np.ones(nlen) * refind

# Solution basis
basis_inv = toast.Raster(mesh_inv, grd)

# Initial projections
phi = mesh_inv.Fields(None, qvec, mua, mus, ref, freq)
proj = projection(phi, mvec)
lnamp = proj[0:nqm]
phase = proj[nqm:nqm*2]

# Data scaling
sd_lnamp = np.ones(lnamp.shape) * np.linalg.norm(lnamp_tgt-lnamp)
sd_phase = np.ones(phase.shape) * np.linalg.norm(phase_tgt-phase)
sd = np.concatenate((sd_lnamp,sd_phase))

# Map parameter estimates to solution basis
bmua = basis_inv.Map('M->B', mua)
bmus = basis_inv.Map('M->B', mus)
bkap = basis_inv.Map('M->B', kap)
bcmua = bmua * cm
bckap = bkap * cm
scmua = basis_inv.Map('B->S', bcmua)
sckap = basis_inv.Map('B->S', bckap)

# Vector of unknowns
x = np.asmatrix(np.concatenate((scmua, sckap))).transpose()
logx = np.log(x)
slen = x.shape[0]/2

# Create regularisation object
#pdb.set_trace()
#hreg = regul.Make ("TK1", hraster, logx, tau);
regul = toast.Regul("TV", basis_inv, logx, tau, beta=beta);

# Initial error
err0 = objective(proj, data, sd, logx)
err = err0
errp = 1e10
erri = np.array([err])
errmua = np.array([imerr(bmua, bmua_tgt)])
errmus = np.array([imerr(bmus, bmus_tgt)])

itr = 1
step = 1.0

hfig1=plt.figure(1)
plt.show()
plt.figure(2)
plt.show()

while itr <= itrmax and err > tolCG*err0 and errp-err > tolCG:
    errp = err
    
    r = -toast.Gradient(mesh_inv.Handle(), basis_inv.Handle(),
                     qvec, mvec, mua, mus, ref, freq, data, sd)
    r = matrix(r).transpose()
    r = np.multiply(r, x)  # parameter scaling

    rr = -regul.Gradient(logx)
    rr = np.transpose(np.matrix(rr))

    r = r + rr
    
    plt.figure(2)
    plt.clf()
    plt.subplot(2,2,1)
    im = plt.imshow (np.reshape (basis_inv.Map('S->B', rr[0:slen]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mua prior gradient")
    plt.colorbar()

    plt.subplot(2,2,2)
    im = plt.imshow (np.reshape (basis_inv.Map('S->B', rr[slen:slen*2]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("kap prior gradient")
    plt.colorbar()

    plt.subplot(2,2,3)
    im = plt.imshow (np.reshape (basis_inv.Map('S->B', r[0:slen]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("mua tot gradient")
    plt.colorbar()

    plt.subplot(2,2,4)
    im = plt.imshow (np.reshape (basis_inv.Map('S->B', r[slen:slen*2]), grd))
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.title("kap tot gradient")
    plt.colorbar()

    plt.draw()
    
    if itr > 1:
        delta_old = delta_new
        delta_mid = np.dot(r.transpose(), s)
        
    s = r # replace this with preconditioner

    if itr == 1:
        d = s
        delta_new = np.dot(r.transpose(), d)
        delta0 = delta_new
    else:
        delta_new = np.dot(r.transpose(), s)
        beta = (delta_new-delta_mid) / delta_old
        if itr % resetCG == 0 or beta <= 0:
            d = s
        else:
            d = s + d*beta

    delta_d = np.dot(d.transpose(), d)
    step,err = toast.Linesearch(logx, d, step, err, objective_ls)

    logx = logx + d*step
    x = np.exp(logx)
    scmua = x[0:slen]
    sckap = x[slen:2*slen]
    smua = scmua/cm
    skap = sckap/cm
    smus = 1/(3*skap) - smua
    mua = basis_inv.Map('S->M', smua)
    mus = basis_inv.Map('S->M', smus)

    bmua = np.reshape(basis_inv.Map('S->B', smua), grd)
    bmus = np.reshape(basis_inv.Map('S->B', smus), grd)

    erri = np.concatenate((erri, [err]))
    errmua = np.concatenate((errmua, [imerr(bmua, bmua_tgt)]))
    errmus = np.concatenate((errmus, [imerr(bmus, bmus_tgt)]))
    print ("Iteration "+str(itr)+", objective "+str(err))

    plt.figure(1)
    plt.clf()
    hfig1.suptitle("Iteration "+str(itr))



    ax1 = hfig1.add_subplot(231)
    im = ax1.imshow(bmua_tgt, vmin=mua_min, vmax=mua_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    ax1.set_title("mua target")
    plt.colorbar(im)

    ax2 = hfig1.add_subplot(232)
    im = ax2.imshow(bmus_tgt, vmin=mus_min, vmax=mus_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    ax2.set_title("mus target")
    plt.colorbar(im)

    ax3 = hfig1.add_subplot(234)
    im = ax3.imshow(bmua, vmin=mua_min, vmax=mua_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    ax3.set_title("mua recon")
    plt.colorbar(im)

    ax4 = hfig1.add_subplot(235)
    im = ax4.imshow(bmus, vmin=mus_min, vmax=mus_max)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    ax4.set_title("mus recon")
    plt.colorbar(im)

    ax5 = hfig1.add_subplot(233)
    im = ax5.semilogy(erri)
    ax5.set_title("objective function")
    plt.xlabel("iteration")
    
    ax6 = hfig1.add_subplot(236)
    im = ax6.semilogy(errmua)
    im = ax6.semilogy(errmus)
    ax6.set_title("rel. image error")
    plt.xlabel("iteration")
    
#    plt.draw()
    plt.pause(0.05)
    
    itr = itr+1

#plt.ioff()

