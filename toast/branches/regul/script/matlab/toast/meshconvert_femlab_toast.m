nlen = size(mesh.p,2);
elen = size(mesh.t,2);
dof = 1;

fid = fopen ('convert.msh','w');
fprintf (fid,'MeshData 5.0\n\n');
fprintf (fid,'NodeList %d %d\n', nlen, dof);

for i=1:nlen
    fprintf (fid, 'N[');
    fprintf (fid, '%f %f %f', mesh.p(:,i)');
    fprintf (fid, ']\n');
end

fprintf (fid, '\nElementList %d\n', elen);

for i=1:elen
    fprintf (fid, 'c');
    fprintf (fid, ' %d %d %d %d', mesh.t(1:4,i)');
    fprintf (fid, '\n');
end

fprintf (fid, '\n[ParameterList]\n');
fprintf (fid, 'Size %d\n', nlen);
fprintf (fid, 'Param1 MUA\n');
fprintf (fid, 'Param2 KAPPA\n');
fprintf (fid, 'Param3 N\n');
fprintf (fid, 'Data\n');
for i=1:nlen
    fprintf (fid, '%f %f %f\n', 0.01, 0.330033, 1.4);
end