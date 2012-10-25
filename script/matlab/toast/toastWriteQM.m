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
% detectors. It should be a sparse binary nm x nq matrix, where any nonzero
% entries linklist(i,j) indicate that a measurement from detector i for
% source j is present.
%
% To create a fully populated link list (all sources connect to all
% detectors), use
%
%   linklist = sparse(ones(nm,nq));
  
  
