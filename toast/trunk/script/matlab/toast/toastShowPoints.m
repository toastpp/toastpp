function p = toastShowPoints (P, sz, col)

n   = size(P,1);
dim = size(P,2);

if dim == 3

    [spx spy spz] = sphere(10);
    spx = spx.*sz;
    spy = spy.*sz;
    spz = spz.*sz;
    %spc = ones(size(spx))*col;

    for i=1:n
        fvc = surf2patch(spx+P(i,1),spy+P(i,2),spz+P(i,3));
        p = patch(fvc,'FaceLighting', 'Phong');
        set(p,'FaceColor',col,'EdgeColor','none');
    end
end
