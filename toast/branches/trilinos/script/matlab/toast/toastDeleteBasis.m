function p = toastDeleteBasis
%toastDeleteBasis     - De-allocate a basis mapper object from dynamic heap.
%
% Synopsis: toastDeleteBasis (hBasis)
%    hBasis: basis mapper handle
%
% Deletes a basis mapper from dynamic heap. hBasis is a handle to the mapper
% obtained by a previous call to toastSetBasis. The handle is no longer valid
% after this call.
