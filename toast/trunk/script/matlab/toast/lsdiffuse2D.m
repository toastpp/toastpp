function [outmat] = lsdiffuse2D (inmat,difflength,diffsteps) 
%lsdiffuse2D          - diffuse smoothing in 2D
%
% input: inmat (typically ls_function)
% output: smoothed form

     size_x=size(inmat,1);
     size_y=size(inmat,2);

     outmat = inmat;

     for diff_count = 1: diffsteps

     mask = zeros(size_x,size_y);
     mask(2:size_x-1,2:size_y-1) =  ...
           outmat(1:size_x-2,2:size_y-1) + ...
           outmat(3:size_x,2:size_y-1) +  ...   
           outmat(2:size_x-1,1:size_y-2) +  ...  
           outmat(2:size_x-1,3:size_y) - ... 
           4.* outmat(2:size_x-1,2:size_y-1) ;


     outmat = outmat + difflength .* mask;

     end





