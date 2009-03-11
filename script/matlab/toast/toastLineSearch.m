function [s, pmin] = toastLineSearch (x0, d, s0, p0, func, varargin)
%toastLineSearch      - 1-D minimisation along a given search direction.
%
%   Synopsis: [s,pmin] = toastLineSearch (x0, d, s0, p0, func, ...)
%     x0:    initial parameter vector
%     d:     search direction
%     s0:    initial step size
%     p0:    initial objective value
%     func:  callback function for objective evaluation. Format:
%            function p = func(x,...)
%     ...:   any additional parameters to be passed to func
%     s:     final step size
%     pmin:  final objective value
%
%   Performs an approximate minimisation of a function Q(x0 + s*d),
%   where x0 is the initial parameter set, d is the search direction, and
%   s is the step size being optimised.
%  
%   x0 and d must be provided as vectors of equal size. The actual size n
%   is of no importance to toastLineSearch, because it simply passes the
%   trial step x0+s*d on to the callback function func.
%
%   toastLineSearch is independent of the actual realisation of Q. The
%   calling function must provide a pointer to a function which evaluates
%   Q for given x.  
%
%   toastLineSearch works by first bracketing the minimum, and then
%   performing a quadratic interpolation step. d is required to be a
%   downhill direction in the vicinity of x0.
  
sl = 0;
pl = p0;
sh = s0;
x = x0 + d*sh;
ph = func (x, varargin{:});
fprintf (1, '--> Step: %f, objective: %f\n', sh, ph);

% bracket the minimum

if ph < pl  % increase step size
    sm = sh;
    pm = ph;
    sh = sh*2;
    x = x0 + d*sh;
    ph = func (x, varargin{:});
    fprintf (1, '--> Step: %f, objective: %f\n', sh, ph);

    while ph < pm
        sl = sm; pl = pm;
        sm = sh; pm = ph;
        sh = sh*2;
        x = x0 + d*sh;
        ph = func (x, varargin{:});
        fprintf (1, '--> Step: %f, objective: %f\n', sh, ph);
    end
else        % decrease step size
    sm = sh/2;
    x = x0 + d*sm;
    pm = func (x, varargin{:});
    fprintf (1, '--> Step: %f, objective: %f\n', sm, pm);

    while pm > pl && sh > 1e-8*s0
        sh = sm; ph = pm;
        sm = sm/2;
        x = x0 + d*sm;
        pm = func (x, varargin{:});
        fprintf (1, '--> Step: %f, objective: %f\n', sm, pm);
    end
end

% quadratic interpolation

a = ((pl-ph)/(sl-sh) - (pl-pm)/(sl-sm)) / (sh-sm);
b = (pl-ph)/(sl-sh) - a*(sl+sh);
s = -b/(2*a);
x = x0 + d*s;
pmin = func (x, varargin{:});
if pmin > pm   % no improvement
    s = sm;
    pmin = pm;
end
fprintf (1, '==> Step: %f, objective: %f [final]\n', s, pmin);
