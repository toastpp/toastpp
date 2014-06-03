from toast import toastmod

def Make (regtype,hraster,x0,tau,beta=1):
    return toastmod.Regul (regtype,hraster,x0,tau,beta)

def Value (hreg, x):
    return toastmod.RegValue (hreg, x)

def Gradient (hreg, x):
    return toastmod.RegGradient (hreg, x)

def HDiag (hreg, x):
    return toastmod.RegHDiag (hreg, x)

