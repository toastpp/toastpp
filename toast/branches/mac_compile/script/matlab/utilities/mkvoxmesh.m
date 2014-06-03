function [vtx idx eltp] = mkvoxmesh (img)

grd = size(img);
bb = [[0 0 0];grd];

% create the full bounding box slab mesh
[vtx idx eltp] = mkslab(bb,grd);
ne = size(idx,1);
nn = size(vtx,1);

% remove unused elements
for i=1:grd(1)
    for j=1:grd(2)
        for k=1:grd(3)
            if img(i,j,k) == 0
                ii = (i-1)*grd(2)*grd(3) + (j-1)*grd(3) + k;
                eltp(ii) = -1;
            end
        end
    end
end
eperm = find(eltp ~= -1);
eltp = eltp(eperm);
idx = idx(eperm,:);

% remove unused vertices
idxused = unique(idx);
vtx = vtx(idxused,:);

nperm = zeros(nn,1);
for i=1:length(idxused)
    nperm(idxused(i)) = i;
end
for i=1:size(idx,1)
    for j=1:size(idx,2)
        if nperm(idx(i,j)) == 0
            error('Consistency problem');
        end
        idx(i,j) = nperm(idx(i,j));
    end
end