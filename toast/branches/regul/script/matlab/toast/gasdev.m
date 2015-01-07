function p = gasdev (sigma)

persistent iset
persistent gset
if isempty(iset)
    iset = 0;
end

r = 10;

if iset == 0
    while r >= 1
        v1 = 2*rand - 1;
        v2 = 2*rand - 1;
        r = v1^2 + v2^2;
    end
    fac = sigma * sqrt(-2*log(r)/r);
    gset = v1*fac;
    iset = 1;
    p = v2*fac;
else
    iset = 0;
    p = gset;
end