function [] = generateInputFile(hMesh, fname, output_file_prefix, specify_QM, ns, source_nodes, g, MHz_freq, is_cosine, directionVector, mua, mus, ref)
%function [] = generateInputFile(hMesh, fname, output_file_prefix, specify_QM, ns, source_nodes, g, MHz_freq, is_cosine, directionVector, mua, mus, ref)
% This MATLAB script generates the input file required to run variable order P_{N} code
% Surya sprerapa@cs.ucl.ac.uk 11/01/2011
% 
% Inputs
%  	hMesh -> toast Mesh handle
%
%	fname -> output file name
%    
%       output_file_prefix -> prefix for output files
%
%	specify_QM -> Is a QM file being specified (1. Yes, 0. No)
%
%	ns -> number of sources ( If a QM file is not specified)
%
%	source_nodes -> A vector of source node numbers (If a QM file is not specified. Zero based indexing)
%
%	g -> value of 'g' in Henvey-Greenstein phase function.
%
%	MHz_freq -> Modulation frequency in MHz.
%
%	is_cosine -> Flag to indicate if the source is cosine or directed (1 if cosine, 0 directed).
%
%	directionVector -> if 'is_cosine' is set to 0, this variable specifies the direction vector along which the source 
%			   is directed. NOTE: THIS VECTOR IS NORMALIZED TO BE A UNIT VECTOR INTERNALLY.
%			   e.g. Inward directed source is specified as [0 0 0] - vtx(source_node_num, :).
%
%	mua -> absorption coefficient (mm^{-1}) for each element
%
%	mus -> scattering coefficient (mm^{-1}) for each element
%
%	ref -> Refractive index for each element
%
[vtx, idx] = toastMeshData(hMesh); 
fid =  fopen(fname, 'w');
fprintf(fid, '%s\n', output_file_prefix);
if(specify_QM == 0)
	fprintf(fid, '%d\n', ns);
	for i = 1 : ns
		fprintf(fid, '%d\n', source_nodes(i));
	end;
end;
fprintf(fid, '%f\n%f\n%d\n%f\n%f\n%f\n', g, MHz_freq, is_cosine, directionVector(1), directionVector(2), directionVector(3));
for i = 1  : size(idx, 1)
fprintf(fid, '%f\n%f\n%f\n', mua(i), mus(i), ref(i));
end;

fclose(fid);
