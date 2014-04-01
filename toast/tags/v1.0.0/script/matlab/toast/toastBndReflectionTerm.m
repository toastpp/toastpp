function p = toastBndReflectionTerm
%toastBndReflectionTerm - Returns refractive index mismatch term for DOT.
%
%   Syntax: zeta = toastBndReflectionTerm(refind,method)
%     refind (scalar): refractive index of medium
%     method (string): 'Keijzer' or 'Contini'
%
%   Returns a scalar value that represents the effect of the refractive index
%   mismatch between the scattering domain (refind) and the surrounding medium
%   (assumed to be 1) for the diffuse optical tomography problem.
%
%   The refractive index mismatch term is required for the assembly of the
%   boundary term in the global stiffness matrix of the DOT term.
%
%   Two methods are implemented. For 'Keijzer', see M. Keijzer et al. AO 27,
%   1820-1824 (1988). For 'Contini', see Contini et al. AO 36, 4587-4599 (1997)
