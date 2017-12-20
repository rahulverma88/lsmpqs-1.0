%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:    Masa Prodanovic
% Copyright: (c) 2009 The University of Texas at Austin. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function writeGeneralSphereBinary(x,y,z,r,x_lo,x_hi,fname)
% 
% Writes information about a number of spheres into a binary file.
%
% Arguments
% x_lo  -(array with 3 values) - lower ends of geometry in x,y,z directions
% x_hi  -(array with 3 values) - upper ends of geometry in x,y,z directions
% x     - array with x coords of sphere centers
% y     - array with y coords of sphere centers
% z     - array with z coords of sphere centers
% r     - array with sphere radii
% fname - file name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeGeneralSphereBinary(x,y,z,r,x_lo,x_hi,fname)

fid = fopen(fname,'w');
n_sphere = numel(x)

fwrite(fid,n_sphere,'int');
fwrite(fid,x_lo,'float');
fwrite(fid,x_hi,'float');
fwrite(fid,x,'float');
fwrite(fid,y,'float');
fwrite(fid,z,'float');
fwrite(fid,r,'float');

fclose(fid);
