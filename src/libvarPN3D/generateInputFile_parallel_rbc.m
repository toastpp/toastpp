function [] = generateInputFile_parallel_rbc(hMesh, fname, output_file_prefix, nQ, nM, specify_QM, source_types, source_profiles, source_widRadii,source_nodes, meas_profiles, meas_widRadii, meas_nodes, g, MHz_freq, is_cosine, directionVector, mua, mus, ref, refout, sphOrder)
%function [] = generateInputFile(hMesh, fname, output_file_prefix,  ns, specify_QM, source_types, source_profiles, source_widRadii,source_nodes, g, MHz_freq, is_cosine, directionVector, mua, mus, ref, sphOrder)
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
%	nQ -> number of sources 
%
%	nM -> number of measurements
%
%	specify_QM -> Is a QM file being specified (1. Yes, 0. No)
%
%	source_types -> nQ x 1 vector containing the source types where 'nQ' is the number of sources.
%		      (1) Neumann boundary source (2) Isotropic point source
%                     NOTE: Ignored when specify_QM = 0. Use [].
%
%	source_profiles -> nQ x 1 vector containing the source profiles where 'nQ' is the number of sources.
%			  (1) Point (2) Gaussian (3) Cosine
%                     NOTE: Ignored when specify_QM = 0. Use [].
%
%	source_widRadii -> nQ x 1 vector containing the source radii (if a 'Point' source) or width (if a 'Gaussian' source) or ignored (if 'Cosine') 
%			   where 'nQ' is the number of sources.
%		      	   NOTE:Ignore when specify_QM = 0. Use [].
%
%	source_nodes -> A vector of source node numbers (If a QM file is not specified)
%
%	meas_profiles -> nM x 1 vector containing the source profiles where 'nM' is the number of measurements.
%			 (1) Point (2) Gaussian (3) Cosine
%                        NOTE: Ignored when specify_QM = 0. Use [].
%
%	meas_widRadii -> nM x 1 vector containing the source radii (if a 'Point' source) or width (if a 'Gaussian' source) or ignored (if 'Cosine') 
%			 where 'nM' is the number of measurements.
%		         NOTE:Ignore when specify_QM = 0 and when meas_profiles = 1. Use [].
%
%	meas_nodes -> A vector of measurement node numbers (If a QM file is not specified)
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
%	refout -> Refractive index of the media outside the domain
%
%	sphOrder -> spherical harmonics order for each node number
[vtx, idx] = toastMeshData(hMesh); 
fid =  fopen(fname, 'w');
fprintf(fid, '%s\n', output_file_prefix);
if(specify_QM == 0)
	fprintf(fid, '%d\n', nQ);
	for i = 1 : nQ
		fprintf(fid, '%d\n', source_nodes(i));
	end;
	fprintf(fid, '%d\n', nM);
	for i = 1 : nM
		fprintf(fid, '%d\n', meas_nodes(i));
	end;
else
	for i = 1 : nQ
		fprintf(fid, '%d\n%d\n', source_types(i), source_profiles(i));
		if(source_profiles(i) == 1 | source_profiles(i) == 2)
			fprintf(fid, '%f\n', source_widRadii(i));
		end;
	end;
	for i = 1 : nM
		fprintf(fid, '%d\n', meas_profiles(i));
		if((meas_profiles(i)-1) > 0)
			disp('Entered here')
			meas_profiles(i) - 1
			fprintf(fid, '%f\n', meas_widRadii(i));
		end;
	end;

end;
fprintf(fid, '%f\n%f\n%d\n%f\n%f\n%f\n', g, MHz_freq, is_cosine, directionVector(1), directionVector(2), directionVector(3));
for i = 1  : size(idx, 1)
fprintf(fid, '%f\n%f\n%f\n', mua(i), mus(i), ref(i));
end;

fprintf(fid, '%f\n', refout);

for i = 1 : size(vtx, 1)
fprintf(fid, '%d\n', sphOrder(i));
end; 
fclose(fid);
