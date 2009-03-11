function p = toastWriteQM
%toastWriteQM         - Write source-detector definitions to a file.
%
% Synposis: toastWriteQM (qpos, mpos, linklist, fname)
%     qpos:     vector of source positions (nq x dim)
%     mpos:     vector of measurement positions (nm x dim)
%     linklist: source-detector link information (nq x nm)
%     fname:    QM output file name
%
% qpos and mpos contain the boundary positions of the source and detector
% locations.
% linklist contains the connectivity information between sources and
% detectors. It should be a sparse nq x nm matrix, where the nonzero
% entries for column i contain the indices (1-based) of all detectors
% containing data for source i.
  
  
