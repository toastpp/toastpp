function bndterm = toastDotBndterm(refind,method)
if nargin > 1
    bndterm = toastmex(uint32(62),refind,method);
else
    bndterm = toastmex(uint32(62),refind);
end

end
