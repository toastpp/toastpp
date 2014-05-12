function p = toastRegul
%toastRegul           - Create a regulariser instance and return a handle.
%
% Synopsis: hReg = toastRegul (type,hBasis,x0,tau[,...])
%    type:   regularisation method (string)
%    hBasis: basis mapper (handle)
%    x0:     initial parameter estimate (vector)
%    tau:    regularisation hyperparameter (scalar)
%    hReg:   regulariser handle
%
% Creates a regulariser instance and returns a handle to it. The handle can
% be used subsequently to invoke the regulariser, e.g. to obtain the 
% regularisation value or gradient for a given parameter distribution.
%
% Depending on regularisaton method, additional parameters may be required
% (see below).
%
% The following regularisation methods are currently supported:
%
%   'DIAG': 0th order Tikhonov (obsolete)
%   -------------------------------------
%   Synopsis:
%       hReg = toastRegul ('DIAG', hBasis, x0, tau)
%
%   'LAPLACIAN': 1st order Tikhonov (obsolete)
%   ------------------------------------------
%   Synopsis:
%       hReg = toastRegul ('LAPLACIAN', hBasis, x0, tau)
%
%   'TV': total variation
%   ---------------------
%   Synopsis:
%       hReg = toastRegul ('TV', hBasis, x0, tau [, Property, value, ...])
%
%       Recognised properties:
%       'Beta':        Next argument is scalar regularisation parameter.
%                      Default: 1.0.
%       'KapRefImage'  Next argument is scalar array used for generating
%                      reference diffusivity field. Image dimensions must
%                      be a power of 2 due to FFT invokation.
%                      Default: none.
%       'KapRefScale'  Next argument is scalar value used for scaling the
%                      reference image. Default: 0.
%       'KapRefPMThreshold'
%                      Fraction of max to use for Perona-Malik threshold
%                      [0-1]. Default: 0
%       'KapRefTensor' Next argument is boolean flag that defines if
%                      reference diffusivity is set up as a tensor field.
%                      Only used in conjunction with KapRefImage.
%                      Default: false.
%       'KapRef'       Next argument is reference diffusivity array.
%                      Supplied either as scalar array (isotropic
%                      diffusivity) or cell array (anisotropic diffusivity),
%                      where each cell contains the diffusivity tensor for
%                      one pixel. This option cannot be combined with
%                      KapRefImage. Default: no diffusivity information.
%
%   'TK0': 0th order Tikhonov regularisation
%   ----------------------------------------
%   Synopsis:
%       hReg = toastRegul ('TK0', hBasis, x0, tau [, Property, value, ...])
%
%       Recognised properties:
%       'Xs'           Next argument is a scaling vector with the same
%                      length as x0.
%
%   Notes:
%       The regularisation value is calculated as
%       psi = tau * <dx,dx>, where dx = (x-x0)/Xs
%  
%   'TK1': 1st order Tikhonov (Laplacian) regularisation
%   ----------------------------------------------------
%   Synopsis:
%       hReg = toastRegul ('TK1', hBasis, x0, tau [, Property, value, ...])
%
%       Recognised properties:
%       'KapRefImage'  Next argument is scalar array used for generating
%                      reference diffusivity field. Image dimensions must
%                      be a power of 2 due to FFT invokation.
%                      Default: none.
%       'KapRefScale'  Next argument is scalar value used for scaling the
%                      reference image. Default: 0.
%       'KapRefPMThreshold'
%                      Fraction of max to use for Perona-Malik threshold
%                      [0-1]. Default: 0
%       'KapRefTensor' Next argument is boolean flag that defines if
%                      reference diffusivity is set up as a tensor field.
%                      Only used in conjunction with KapRefImage.
%                      Default: false.
%       'KapRef'       Next argument is reference diffusivity array.
%                      Supplied either as scalar array (isotropic
%                      diffusivity) or cell array (anisotropic diffusivity),
%                      where each cell contains the diffusivity tensor for
%                      one pixel. This option cannot be combined with
%                      KapRefImage. Default: no diffusivity information.
%
%   'HUBER':  Huber-function regularisation
%   ---------------------------------------
%   Synopsis:
%       hReg = toastRegul ('HUBER', hBasis, x0, tau [, Property, value, ...])
%
%       Recognised properties:
%       'Eps':         Next argument is scalar regularisation parameter.
%                      Default: 1.0.
%       'KapRefImage'  Next argument is scalar array used for generating
%                      reference diffusivity field. Image dimensions must
%                      be a power of 2 due to FFT invokation.
%                      Default: none.
%       'KapRefScale'  Next argument is scalar value used for scaling the
%                      reference image. Default: 0.
%       'KapRefPMThreshold'
%                      Fraction of max to use for Perona-Malik threshold
%                      [0-1]. Default: 0
%       'KapRefTensor' Next argument is boolean flag that defines if
%                      reference diffusivity is set up as a tensor field.
%                      Only used in conjunction with KapRefImage.
%                      Default: false.
%       'KapRef'       Next argument is reference diffusivity array.
%                      Supplied either as scalar array (isotropic
%                      diffusivity) or cell array (anisotropic diffusivity),
%                      where each cell contains the diffusivity tensor for
%                      one pixel. This option cannot be combined with
%                      KapRefImage. Default: no diffusivity information.
%
%   'PM': Perona-Malik regularisation
%   ---------------------------------
%   Synopsis:
%       hReg = toastRegul ('PM', hBasis, x0, tau [, Property, value, ...])
%
%       Recognised properties:
%       'T':           Next argument is scalar regularisation parameter.
%                      Default: 1.0.
%       'KapRefImage'  Next argument is scalar array used for generating
%                      reference diffusivity field. Image dimensions must
%                      be a power of 2 due to FFT invokation.
%                      Default: none.
%       'KapRefScale'  Next argument is scalar value used for scaling the
%                      reference image. Default: 0.
%       'KapRefPMThreshold'
%                      Fraction of max to use for Perona-Malik threshold
%                      [0-1]. Default: 0
%       'KapRefTensor' Next argument is boolean flag that defines if
%                      reference diffusivity is set up as a tensor field.
%                      Only used in conjunction with KapRefImage.
%                      Default: false.
%       'KapRef'       Next argument is reference diffusivity array.
%                      Supplied either as scalar array (isotropic
%                      diffusivity) or cell array (anisotropic diffusivity),
%                      where each cell contains the diffusivity tensor for
%                      one pixel. This option cannot be combined with
%                      KapRefImage. Default: no diffusivity information.
%
%   'QPM': Quadratic Perona-Malik regularisation
%   --------------------------------------------
%   Synopsis:
%       hReg = toastRegul ('QPM', hBasis, x0, tau [, Property, value, ...])
%
%       Recognised properties:
%       'T':           Next argument is scalar regularisation parameter.
%                      Default: 1.0.
%       'KapRefImage'  Next argument is scalar array used for generating
%                      reference diffusivity field. Image dimensions must
%                      be a power of 2 due to FFT invokation.
%                      Default: none.
%       'KapRefScale'  Next argument is scalar value used for scaling the
%                      reference image. Default: 0.
%       'KapRefPMThreshold'
%                      Fraction of max to use for Perona-Malik threshold
%                      [0-1]. Default: 0
%       'KapRefTensor' Next argument is boolean flag that defines if
%                      reference diffusivity is set up as a tensor field.
%                      Only used in conjunction with KapRefImage.
%                      Default: false.
%       'KapRef'       Next argument is reference diffusivity array.
%                      Supplied either as scalar array (isotropic
%                      diffusivity) or cell array (anisotropic diffusivity),
%                      where each cell contains the diffusivity tensor for
%                      one pixel. This option cannot be combined with
%                      KapRefImage. Default: no diffusivity information.
%
%   'Tukey': Tukey regularisation
%   -----------------------------
%   Synopsis:
%       hReg = toastRegul ('Tukey', hBasis, x0, tau [, Property, value, ...])
%
%       Recognised properties:
%       'T':           Next argument is scalar regularisation parameter.
%                      Default: 1.0.
%       'KapRefImage'  Next argument is scalar array used for generating
%                      reference diffusivity field. Image dimensions must
%                      be a power of 2 due to FFT invokation.
%                      Default: none.
%       'KapRefScale'  Next argument is scalar value used for scaling the
%                      reference image. Default: 0.
%       'KapRefPMThreshold'
%                      Fraction of max to use for Perona-Malik threshold
%                      [0-1]. Default: 0
%       'KapRefTensor' Next argument is boolean flag that defines if
%                      reference diffusivity is set up as a tensor field.
%                      Only used in conjunction with KapRefImage.
%                      Default: false.
%       'KapRef'       Next argument is reference diffusivity array.
%                      Supplied either as scalar array (isotropic
%                      diffusivity) or cell array (anisotropic diffusivity),
%                      where each cell contains the diffusivity tensor for
%                      one pixel. This option cannot be combined with
%                      KapRefImage. Default: no diffusivity information.
%
%
% Note: The 'kref' parameter, if provided, defines a diffusivity
% distribution for the regularisation. This can be either isotropic, in
% which case kref is a scalar field, passed as a real array, or anisotropic,
% in which case kref is a tensor field, passed as a cell array, where
% each cell i contains a 2-D matrix defining the diffusivity of pixel i.
  
