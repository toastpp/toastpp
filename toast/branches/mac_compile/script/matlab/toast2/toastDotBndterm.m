function bndterm = toastDotBndterm(refind,method)
if nargin > 1
    bndterm = toast(uint32(62),refind,method);
else
    bndterm = toast(uint32(62),refind);
end

end
