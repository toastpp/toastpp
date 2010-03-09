function p = toastWriteQM
%toastWriteQM         - Write source-detector definitions to a file.
%
% Synposis: toastWriteQM (qpos, mpos, linklist, fname)
%     qpos:     vector of source positions (nq x dim)
%     mpos:     vector of measurement positions (nm x dim)
%     linklist: source-detector link information (nm x nq)
%     fname:    QM output file name
%
% qpos and mpos contain the boundary positions of the source and detector
% locations.
% linklist contains the connectivity information between sources and
% detectors. It should be a sparse binary nm x nq matrix, where the nonzero
% entries for column j contain the indices (1-based) of all detectors
% containing data for source j.
%
% To create a fully populated link list (all sources connect to all
% detectors), use
%
%   linklist = sparse(ones(nm,nq));
  
  
