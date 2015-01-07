clear all
close all

generate_ref = false;

% Make sure the toast paths are available
cwd = pwd;
cd (getenv('TOASTDIR'));
mtoast_install(true);
cd (cwd);

if ~generate_ref
    load unit_test.mat
end

meshdir = '../meshes/';

fprintf('toast toolbox unit tests\n\n')

%% test: toastReadRealVector+toastWriteRealVector
a = [1:3:300]';
toastWriteRealVector('tmp.dat', a);
b = toastReadRealVector('tmp.dat');
err = norm(a-b)/norm(a);
if err < 1e-8
    fprintf('PASSED ReadRealVector/WriteRealVector: rel. error %f\n', err);
else
    error('ReadRealVector/WriteRealVector: rel. error %f', err);
end
clear a
clear b
delete tmp.dat

%% test: toastReadComplexVector+toastWriteComplexVector
a = ([1:3:300] + 1i*[1:3:300]).';
toastWriteComplexVector('tmp.dat',a);
b = toastReadComplexVector('tmp.dat');
err = norm(a-b)/norm(a);
if err < 1e-8
    fprintf('PASSED ReadComplexVector/WriteComplexVector: rel. error %f\n', err);
else
    error('ReadComplexVector/WriteComplexVector: rel. error %f', err);
end
clear a
clear b
delete tmp.dat

%% test: toastReadMesh+toastMeshData
hmesh = toastReadMesh([meshdir 'circle25_32.msh']);
[vtx idx eltp] = toastMeshData(hmesh);

if generate_ref
    test.ReadMeshMeshData1.vtx = vtx;
    test.ReadMeshMeshData1.idx = idx;
    test.ReadMeshMeshData1.eltp = eltp;
else
    err = max (norm (test.ReadMeshMeshData1.vtx-vtx));
    if err < 1e-8
        fprintf('PASSED ReadMesh/MeshData: vertex coordinate max displacement %f\n', err)
    else
        error('ReadMesh/MeshData: vertex coordinates displaced (%f)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.idx-idx)));
    if err==0
        fprintf('PASSED ReadMesh/MeshData: element node indices match\n')
    else
        error('ReadMesh/MeshData: element node index mismatch (%d)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.eltp-eltp)));
    if err==0
        fprintf('PASSED ReadMesh/MeshData: element type indices match\n')
    else
        error('ReadMesh/MeshData: element type index mismatch (%d)', err)
    end
end

%% test: toastSurfData
[svtx sidx perm] = toastSurfData(hmesh);

if generate_ref
    test.SurfData1.svtx = svtx;
    test.SurfData1.sidx = sidx;
    test.SurfData1.perm = perm;
else
    err = max (norm (test.SurfData1.svtx-svtx));
    if err < 1e-8
        fprintf ('PASSED SurfData: vertex coordinate max displacement: %f\n', err)
    else
        error('SurfData: vertex coordinates displaced (%f)', err)
    end
    err = max(max(abs(test.SurfData1.sidx-sidx)));
    if err==0
        fprintf('PASSED SurfData: element node indices match\n')
    else
        error('SurfData: element node index mismatch (%d)', err)
    end
    err = max(max(abs(test.SurfData1.perm-perm)));
    if err==0
        fprintf('PASSED SurfData: permutation indices match\n')
    else
        error('SurfData: permutation index mismatch (%d)', err)
    end
end

%% test: toastMakeMesh
hmesh2 = toastMakeMesh (vtx, idx, eltp);
[vtx2 idx2 eltp2] = toastMeshData (hmesh2);

if ~generate_ref
    err = max (norm (test.ReadMeshMeshData1.vtx-vtx2));
    if err < 1e-8
        fprintf('PASSED MakeMesh: vertex coordinate max displacement %f\n', err)
    else
        error('MakeMesh: vertex coordinates displaced (%f)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.idx-idx2)));
    if err==0
        fprintf('PASSED MakeMesh: element node indices match\n')
    else
        error('MakeMesh: element node index mismatch (%d)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.eltp-eltp2)));
    if err==0
        fprintf('PASSED MakeMesh: element type indices match\n')
    else
        error('MakeMesh: element type index mismatch (%d)', err)
    end
end
toastDeleteMesh(hmesh2);

%% test: toastWriteMesh
toastWriteMesh(hmesh,'tmp.msh');
hmesh2 = toastReadMesh('tmp.msh');
[vtx2 idx2 eltp2] = toastMeshData (hmesh2);

if ~generate_ref
    err = max (norm (test.ReadMeshMeshData1.vtx-vtx2));
    if err < 1e-8
        fprintf('PASSED WriteMesh: vertex coordinate max displacement %f\n', err)
    else
        error('WriteMesh: vertex coordinates displaced (%f)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.idx-idx2)));
    if err==0
        fprintf('PASSED WriteMesh: element node indices match\n')
    else
        error('WriteMesh: element node index mismatch (%d)', err)
    end
    err = max(max(abs(test.ReadMeshMeshData1.eltp-eltp2)));
    if err==0
        fprintf('PASSED WriteMesh: element type indices match\n')
    else
        error('WriteMesh: element type index mismatch (%d)', err)
    end
end
toastDeleteMesh(hmesh2);
delete tmp.msh

%% test: toastReadQM/toastQPos/toastMPos
toastReadQM(hmesh,[meshdir 'circle25_32x32.qm']);
qpos = toastQPos(hmesh);
mpos = toastMPos(hmesh);

if generate_ref
    test.ReadQMQMPos1.qpos = qpos;
    test.ReadQMQMPos1.mpos = mpos;
else
    err = max (norm (test.ReadQMQMPos1.qpos-qpos));
    if err < 1e-8
        fprintf('PASSED ReadQM/QPos/MPos: source coordinate max displacement %f\n', err)
    else
        error('ReadQM/QPos/MPos: source coordinates displaced (%f)', err)
    end
    err = max (norm (test.ReadQMQMPos1.mpos-mpos));
    if err < 1e-8
        fprintf('PASSED ReadQM/QPos/MPos: detector coordinate max displacement %f\n', err)
    else
        error('ReadQM/QPos/MPos: detector coordinates displaced (%f)', err)
    end
end

%% Output reference data
if generate_ref
    save unit_test.mat test
    fprintf('\ntest reference data saved\n')
else
    fprintf('\nall tests passed\n')
end