function unim = unwrap_phase (mesh, nim, seed)
% Given a mesh and a phase field 'nim', unwrap the phase data starting
% from a seed point 'seed'.

[vtx,idx] = mesh.Data;
n = mesh.NodeCount;
nn = mesh.NodeNeighbour;
nnn = zeros(n,1);
for i=1:n
    nnn(i) = length(find(nn(i,:)));
end

for i=1:n
    dst = norm(vtx(i,:)-seed);
    if i==1 || dst < mindst
        mindst = dst;
        seednd = i;
    end
end

unim = nim;
processed = zeros(n,1);
processed(seednd) = 1;

frontline = seednd;

while length(frontline) > 0
    frontline = process_neighbours(frontline);
end


    function newfront = process_neighbours (front)
        newfront = [];
        for i_=1:length(front)
            s = front(i_);
            for j_=1:nnn(s)
                nb = nn(s,j_);
                if processed(nb) == 0
                    while unim(nb) >= unim(s)+pi
                        unim(nb) = unim(nb) - 2*pi;
                    end
                    while unim(nb) < unim(s)-pi
                        unim(nb) = unim(nb) + 2*pi;
                    end
                    processed(nb) = 1;
                    newfront = [newfront, nb];
                end
            end
        end
    end

end