function [] = generateInputFile(hMesh, fname, output_file_prefix, source_node_num, g, MHz_freq, is_isotropic, directionVector, mua_3layers, mus_3layers)
% This MATLAB script generates the input file required to run variable order P_{N} code on the 3 layered sphere geometry
% Surya sprerapa@cs.ucl.ac.uk 27/05/2010
% 
% Inputs
%  	hMesh -> toast Mesh handle
%
%	fname -> output file name
%
%	source_node_num -> boundary node number around which the source is placed (Currently deals with only a single source)
%			   e.g. node number 1669 is chosen as source node for our test cases.
%
%	g -> value of 'g' in Henvey-Greenstein phase function.
%
%	MHz_freq -> Modulation frequency in MHz.
%
%	is_isotropic -> Flag to indicate if the source is isotropic (1 if isotropic, 0 otherwise).
%
%	directionVector -> if 'is_isotropic' is set to 0, this variable specifies the direction vector along which the source 
%			   is directed. NOTE: THIS VECTOR IS NORMALIZED TO BE A UNIT VECTOR INTERNALLY.
%			   e.g. Inward directed source is specified as [0 0 0] - vtx(source_node_num, :).
%
%	mua_3layers -> absorption coefficient (mm^{-1}) for each of the 3 layers
%		       e.g. [0.01, 0.025, 0.01].
%
%	mus_3layers -> scattering coefficient (mm^{-1}) for each of the 3 layers
%			e.g. [5, 0.1, 5].
[vtx, idx] = toastMeshData(hMesh); 
[mua, mus, ref] = assignMatValues(vtx, idx, mua_3layers, mus_3layers, [1, 1, 1]);


fid =  fopen(fname, 'w');
fprintf(fid, '%s\n%d\n%f\n%d\n%d\n%f\n%f\n%f\n', output_file_prefix, source_node_num, g, MHz_freq, is_isotropic, directionVector(1), directionVector(2), directionVector(3));
for i = 1  : size(idx, 1)
fprintf(fid, '%f\n%f\n%f\n', mua(i), mus(i), ref(i));
end;

fclose(fid);
