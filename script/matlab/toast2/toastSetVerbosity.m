%toastSetVerbosity   - set the level of diagnostic output
%
% Synopsis: toastSetVerbosity(level)
%
%   level: verbosity level (integer >= 0), where 0 is no output, larger
%          numbers provide more output.

function toastSetVerbosity(level)
toast(uint32(1000),uint32(level));
end
