function p = toastGradient(hmesh,hbasis,qvec,mvec,mua,mus,ref,freq,data,sd,varargin)
% toastGradient        - Return gradient of objective function.
%
% Syntax: g = toastGradient(mesh,basis,qvec,mvec,mua,mus,ref,freq,data,sd, ...)
%
% Parameters:
%         mesh (toastMesh instance):
%             mesh object
%         basis (toastBasis instance):
%             basis mapper object
%         qvec (complex matrix n x nq):
%             array of nodal source vectors
%         mvec (complex matrix n x nm):
%             array of nodal measurement vectors
%         mua (real array n):
%             array of nodal absorption coefficients [1/mm]
%         mus (real array n):
%             array of nodal scattering coefficients [1/mm]
%         ref (real array n):
%             array of nodal refractive indices
%         freq (scalar):
%             modulation frequency [MHz]
%         data (real array nqm*2):
%             array of measurement data (lnamp, phase)
%         sd (real array nqm*2):
%             array of standard deviations (lnamp, phase)
%
% Optional parameters:
%         'Method', method (string):
%             forward solver method
%             (DIRECT|CG|BICG|BICGSTAB|GMRES|GAUSSSEIDEL)
%         'Tolerance', tol (scalar):
%             forward solver tolerance
%         'Fields', phi (complex matrix n x nq):
%             array of nodal fields for each source
%         'Projections', proj (real array nqm*2):
%             array of projection data (lnamp, phase)
%         'Unwrap', true|false (boolean):
%             flag for activating phase unwrapping
%
% Return values:
%         g (real array slen*2):
%             gradient of objective function (mua, mus) in inverse basis
%
% This function runs the forward solver with the nodal optical parameters
% provided, and calculates the derivative of the objective function with
% respect to the optical parameters.
% The returned gradient vector is with respect to the inverse solution basis
% defined by basis, rather than the nodal mesh basis.
%
% Note:  If the modulation frequency (freq) is set to 0, toastGradient
%        computes the gradient for a steady-state problem. In that case,
%        data, sd and proj (if provided) can be of length nqm and only
%        contain the log amplitude component, omitting the phase component.
%        (if a phase component is present, it is ignored).
%        Likewise, qvec, mvec and phi (if provided) can be real-valued, as
%        the imaginary component is ignored.

p = toastmex(uint32(52),hmesh.handle,hbasis.handle,qvec,mvec,mua,mus,ref,freq,data,sd,varargin{:});
