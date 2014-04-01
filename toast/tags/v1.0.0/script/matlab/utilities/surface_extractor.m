function [srf] = surface_extractor(tri, vtx, compartment)
% Usage: [srf] = surface_extractor(tri, vtx, compartment);
%
% e.g.
%
%[srf_out_1] = surface_extractor(tri, vtx, [m{4} m{3} m{2} m{1}]);
%[srf_out_2] = surface_extractor(tri, vtx, [m{4} m{3} m{2}]);
%[srf_out_3] = surface_extractor(tri, vtx, [m{4} m{3}]); 
%[srf_out_4] = surface_extractor(tri, vtx, m{4});
%
% General:
% Auxiliary function that extract the boundary faces of a given 3D volume
%
% Input: 
% tri - simplices matrix {k x 4}
% vtx - vertices matrix {n x 3}
% compartment - list of simplces for extraction {no. of elements in comparment x 1}
%
% Output:
% srf - outer boundary surfaces {no. of external surfaces x 3}
%------------------------------------------------------------------------------------------------------------------------

% create a list of all element facetts & sort
S = [tri(compartment,[1 2 3]); tri(compartment,[1 2 4]); tri(compartment,[2 3 4]); tri(compartment,[1 3 4])];
clear tri
N = sort(S,2);
clear S;
M = sortrows(N);
clear N;

i = 1;
inc = 1;

while i < size(M,1);
    ithrow = M(i,:);
    jthrow = M(i+1,:); 
    if ithrow == jthrow 
        i = i + 2;
    else
        % put in boundary node list
        srf(inc,:) = M(i,:);
        inc = inc + 1;
        i = i + 1;
    end
end

fig_title='extracted surface';
f1=figure('Name',fig_title);

trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'FaceAlpha',0.5);
colormap([0 0 0]);
daspect([1 1 1]);

hold on;
axis image;
col = rand(1,3)/2 + [ 0.2 0.2 0.2];
set(gcf,'Colormap',col);
grid off
view(3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This might become a part of the EIDORS suite
% Copyright (c) Lior Horesh 2005, EIT group, Medical Physics and Bioengineering, UCL, UK
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details
% EIDORS 3D version 2
% MATLAB Version 6.1.0.450 (R12.1) on PCWIN
% MATLAB License Number: 111672
% Operating System: Microsoft Windows Server 2003 Standard Edition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
