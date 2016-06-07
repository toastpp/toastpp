function toastWriteMeshGmshPos(mesh,fname,prm)

[vtx,idx,eltp] = mesh.Data;
nel = size(idx,1);
fid = fopen(fname,'wt');
fprintf (fid, 'View "%s" {\n', fname);

for i=1:nel
    switch eltp(i)
        case 15 % 3-noded triangle
            fprintf(fid, 'ST(%f,%f,%f,%f,%f,%f,%f,%f,%f)', ...
                vtx(idx(i,1),1), vtx(idx(i,1),2), 0, ...
                vtx(idx(i,2),1), vtx(idx(i,2),2), 0, ...
                vtx(idx(i,3),1), vtx(idx(i,3),2), 0);
            fprintf(fid, '{%f,%f,%f};\n', ...
                prm(idx(i,1)), prm(idx(i,2)), prm(idx(i,3)));
        case 3 % 4-noded tetrahedron
            fprintf(fid, 'SS(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f)', ...
                vtx(idx(i,1),1), vtx(idx(i,1),2), vtx(idx(i,1),3), ...
                vtx(idx(i,2),1), vtx(idx(i,2),2), vtx(idx(i,2),3), ...
                vtx(idx(i,3),1), vtx(idx(i,3),2), vtx(idx(i,3),3), ...
                vtx(idx(i,4),1), vtx(idx(i,4),2), vtx(idx(i,4),3));
            fprintf(fid, '{%f,%f,%f,%f};\n', ...
                prm(idx(i,1)), prm(idx(i,2)), prm(idx(i,3)), prm(idx(i,4)));            
    end
end
fprintf (fid, '};\n');
fclose(fid);

end