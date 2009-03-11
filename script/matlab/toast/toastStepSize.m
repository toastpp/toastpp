function [step,err] = toastStepSize (x0, d, s0, p0, func, varargin)
%toastStepSize        - Find downhill step in a given search direction.
%
%   Synopsis: [s,pmin] = toastStepSize (x0, d, s0, p0, func, ...)
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
%   Returns a step size <= s0 which reduces the initial objective function
%   p0.
%
%   Direction d is assumed to be a downhill direction in the vicinity of
%   x0.
%
%   This function only decreases, but never increases the step size.
%   It can be used as a simple alternative to toastLineSearch.

maxstep = 5;
step = s0;

for trystep = 1:maxstep
    err = func (x0 + d*step, varargin{:});
    fprintf (1, '--> Step: %f, objective: %f\n', step, err);
    if err < p0
        return;
    end
    step = step/2; % no improvement
end
