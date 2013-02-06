function vec = toastReadVector (fname,idx)
    if nargin < 2
        idx = -1;
    end
    vec = toast(uint32(42), fname,idx);
end