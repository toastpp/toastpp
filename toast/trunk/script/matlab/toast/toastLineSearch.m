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
%   Function 'func' is called by toastLineSearch to evaluate the objective
%   function at a given parameter vector x. It has the format
%
%               p = func(x,...)
%
%   where first argument x is a parameter vector, and any subsequent
%   arguments are the optional arguments passed to toastLineSearch. The
%   return value p can have one of two formats:
%
%   - a real number: this is the value of the objective function at x.
%   - a structure: in this case, the value of the objective function must
%     be returned in field 'of', and an optional boolean field 'valid'
%     can be used to indicate if the parameter distribution x is within
%     the valid range. If field 'valid' is not present, or if p is a
%     number, then valid=true is assumed.
%
%   If p0 is empty ([]), it is evaluated by a call to p0 = func(x0,..)
%
%   toastLineSearch is independent of the actual realisation of Q. The
%   calling function must provide a pointer to a function which evaluates
%   Q for given x.  
%
%   toastLineSearch works by first bracketing the minimum, and then
%   performing a quadratic interpolation step. d is required to be a
%   downhill direction in the vicinity of x0.
  
verbose = false;

nvar = length(varargin);
i = 1;
while i <= nvar
    if ischar(varargin{i}) && i < nvar
        if strcmpi(varargin{i},'verbose')
            verbose = varargin{i+1};
            varargin = varargin([1:i-1,i+2:end]);
            nvar = nvar-2;
            i = i-1;
        end
    end
    i = i+1;
end

if isempty(p0)
    [p0,valid] =  get_of (x0, varargin{:});
end

sl = 0;
pl = p0;
sh = s0;
x = x0 + d*sh;
[ph,valid] = get_of (x, varargin{:});
if verbose
    fprintf (1, '--> Step: %f, objective: %f\n', sh, ph);
end

% bracket the minimum

if ph < pl && valid  % increase step size
    sm = sh;
    pm = ph;
    sh = sh*2;
    x = x0 + d*sh;
    [ph,valid] = get_of (x, varargin{:});
    if verbose
        fprintf (1, '--> Step: %f, objective: %f\n', sh, ph);
    end

    while ph < pm && valid
        sl = sm; pl = pm;
        sm = sh; pm = ph;
        sh = sh*2;
        x = x0 + d*sh;
        [ph,valid] = get_of (x, varargin{:});
        if verbose
            fprintf (1, '--> Step: %f, objective: %f\n', sh, ph);
        end
    end
else        % decrease step size
    sm = sh/2;
    x = x0 + d*sm;
    [pm,valid] = get_of (x, varargin{:});
    if verbose
        fprintf (1, '--> Step: %f, objective: %f\n', sm, pm);
    end
    
    while (pm > pl || ~valid) && sh > 1e-8*s0
        sh = sm; ph = pm;
        sm = sm/2;
        x = x0 + d*sm;
        [pm,valid] = get_of (x, varargin{:});
        if verbose
            fprintf (1, '--> Step: %f, objective: %f\n', sm, pm);
        end
    end
end

if ph < pm
    
    % minimum is beyond furthest sample (i.e outside valid range)
    % so we take the last known valid sample instead
    pmin = pm;
    s = sm;
    
else
    % quadratic interpolation

    a = ((pl-ph)/(sl-sh) - (pl-pm)/(sl-sm)) / (sh-sm);
    b = (pl-ph)/(sl-sh) - a*(sl+sh);
    s = -b/(2*a);
    if verbose
        fprintf (1, '==> Quadratic interpolation in LS, step: %f\n', s);
    end
    x = x0 + d*s;
    [pmin,valid] = get_of (x, varargin{:});
    if pmin > pm || ~valid   % no improvement
        s = sm;
        pmin = pm;
    end
end

if verbose
    fprintf (1, '==> Step: %f, objective: %f [final]\n', s, pmin);
end

    % ===================================================
    function [p,valid] = get_of(x,varargin)
        p = func(x,varargin{:});
        if isstruct(p)
            if isfield(p,'valid')
                valid = p.valid;
            else
                valid = true;
            end
            p = p.of;
        elseif isnan(p)
            valid = false;
        else
            valid = true;
        end
    end
    % ===================================================

end
