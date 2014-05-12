% This scripts demonstrates the use of the toastMesh class.

%% Constructing a mesh object

hmesh1 = toastMesh;
% create a mesh object without an associated mesh

hmesh2 = toastMesh('../../test/2D/meshes/circle25_32.msh');
% create a mesh object from a mesh file

[vtx idx eltp] = mkslab([-10, -10, -10; 10, 10, 10],[10 10 10]);
hmesh3 = toastMesh(vtx,idx,eltp);
% create a mesh object from geometry information


%% Replacing an existing mesh

hmesh2.Read('../../test/2D/meshes/circle25_64.msh');
% replace mesh with data read from file

[vtx idx eltp] = mkcircle(25,6,32,2);
hmesh3.Make(vtx,idx,eltp);
% replace mesh with new geometry


%% Retrieving mesh information

fprintf ('Mesh dimension is %d\n', hmesh3.Dimension);
% mesh dimension

fprintf ('Node count is %d\n', hmesh3.NodeCount);
% number of nodes

fprintf ('Element count is %d\n', hmesh3.ElementCount);
% number of elements

fprintf ('Mesh bounding box is\n');
hmesh3.BoundingBox
% mesh bounding box

fprintf ('Mesh size (volume or area) is %f\n', hmesh3.Size)
% mesh volume (3D meshes) or area (2D meshes)

elsize = hmesh3.ElementSize;
fprintf ('Element sizes %f to %f\n', min(elsize), max(elsize));
fprintf ('Sum of element sizes is %f\n', sum(elsize));
% element sizes

[mvtx midx mtp] = hmesh3.Data;
% retrieve mesh geometry

[svtx sidx perm] = hmesh3.SurfData;
% retrieve surface mesh geometry

[evtx eidx etp] = hmesh3.ElementData(50);
% retrieve geometry for an individual element

elid = hmesh3.FindElement([5,5;10,-3]);
% returns the indices of elements containing a set of points


%% Computing element matrices
intF   = hmesh3.Elmat(1,'F')   % integral of shape functions over element
intFF  = hmesh3.Elmat(1,'FF')  % integral of product of two shape functions
intFFF = hmesh3.Elmat(1,'FFF') % integral of product of two shape functions
intDD  = hmesh3.Elmat(1,'DD')  % integral of product of two shape derivatives
intFD  = hmesh3.Elmat(1,'FD')  % integral of product of shape function and derivative
intFDD = hmesh3.Elmat(1,'FDD') % integral of product of shape function and two derivatives
intdd  = hmesh3.Elmat(1,'dd')  % integral of product of partial derivatives
bintFF = hmesh3.Elmat(1,'BndFF',1) % integral of product of two shape functions over an edge


%% Displaying meshes

hmesh3.Display
% Show the mesh geometry


%% Copying mesh handles

hmesh1 = hmesh2;
% Copying a mesh handle creates a shallow copy (by reference)

hmesh1.NodeCount
hmesh2.Make(vtx,idx,eltp);
hmesh1.NodeCount
% Changes made to one of the meshes affect also the copy

clear hmesh2
hmesh1.NodeCount
% If one of the meshes is cleared, the other remains valid
% Only when all references to a mesh object are cleared is the
% mesh deallocated.


%% Deleting mesh objects

clear hmesh1 hmesh2 hmesh3
% delete mesh objects from memory