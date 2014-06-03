function p = toastDataLinkList
%toastDataLinkList    - Returns permutation index array for data vector.
%
% Synopsis: perm = toastDataLinkList (hMesh)
%
%    hMesh:  mesh handle
%    perm:   permutation index vector
%
%
% Given a measurement setup with n sources and m detectors, a maximum of
% M = n x m measurements could be obtained if all detectors were recorded
% for every source. In practice, the measurement set is often sparse, i.e.
% only a subset of detectors is recorded for any given source and the total
% number of measurements is M' <= M.
%
% toastDataLinkList allows to map between the full data array of size M and
% the sparse data array of size M'. The permutation array perm is of size
% M' and contains the indices in the range 1...M into the full data array.
%
% Therefore
%
%   perm = toastDataLinkList(hMesh);
%   D1 = D(perm);
%
% extracts the actual measurement data from a full data array D into a
% condensed data array D1.
%
% Conversely,
%
%   D = zeros(n*m,1);
%   D(perm) = D1;
%
% fills the actually utilised measurements of the full data vector D with
% the values of the condensed vector D1, leaving the unused entries at
% zero.
