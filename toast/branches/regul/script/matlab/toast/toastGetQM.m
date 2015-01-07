function p = toastGetQM
%toastGetQM           - Extract QM link list from a QMMesh as a sparse matrix.
%
% Synopsis: LL = toastGetQM(hMesh)
%    hMesh:  mesh handle
%    LL:     link list (sparse matrix)
%
% Each column of the returned matrix corresponds to a source, and each row to
% a detector of the measurement setup defined in the QM file loaded for the
% mesh.
% Nonzero entries exist only for the source-detector combinations for which
% measurements were taken and which are registered in the link list.
% The nonzero entries are numbered consecutively (starting from 1) and denote
% the measurement position in the linear measurement vector.
%
% Usage examples:
%
% find(LL(:,i))
% Returns the detector indices used by source i. This is equivalent to row i
% of the QM file link list, except that indices are 1-based.
%
% nonzeros(LL(:,i))
% Contains the indices in the linear data vector of the measurements for
% source i.
%
% nonzeros(LL(j,:))
% Contains the indices in the linear data vector of the measurements for
% detector j.
%
% data_q1 = data(nonzeros(LL(:,1)))
% Extracts the measurements for source 1 from linear data vector 'data'.
