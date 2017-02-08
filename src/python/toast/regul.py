from toast import toastmod

class Regul:

    def __init__(self, regtype, raster, x0, tau, beta=1):
        self.handle = toastmod.Regul(regtype, raster.handle, x0,
                                     tau, beta)

    def Value(self, x):
        return toastmod.RegValue(self.handle, x)

    def Gradient(self, x):
        return toastmod.RegGradient(self.handle, x)

    def HDiag(self, x):
        return toastmod.RegHDiag(self.handle, x)

