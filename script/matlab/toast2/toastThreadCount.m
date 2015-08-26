function [th,ncore] = toastThreadCount (th_set)
% toastThreadCount - Set and return number of threads used by Toast toolbox
%
% Syntax: [th,ncore] = toastThreadCount
%         toastThreadCount(th_set)
%
% Parameters:
%         th_set: new number of threads (see notes)
%
% Return values:
%         th: current number of worker threads
%         ncore: number of cores available on the machine
%
% Notes:  If th_set is an integer >= 1, create this number of worker threads
%         If th_set==0, create as many workers as cores available on the machine
%         If 0 < th_set < 1, it defines the number of workers as a fraction of
%         the total number of cores available.

[th0,ncore0] = toastmex(uint32(1001));
if nargin > 0
    if th_set == 0
        th_set = ncore0;
    elseif th_set < 1
        th_set = max(1, round(th_set*ncore0));
    end
    toastmex(uint32(1001),th_set);
else
    th = th0;
    if nargout > 1
        ncore = ncore0;
    end
end
