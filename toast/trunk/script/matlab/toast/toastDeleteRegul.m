function p = toastDeleteRegul
%toastDeleteRegul     - De-allocate a regularisation object from dynamic heap.
%
% Synopsis: toastDeleteRegul (hReg)
%    hReg: regularisation instance handle
%
% Deletes a regularisation object from dynamic heap. hReg is a handle to the
% regularisation instance obtained by a previous call to toastRegul.
% The handle is no longer valid after this call.
