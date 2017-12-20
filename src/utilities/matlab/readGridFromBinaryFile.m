%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:    Masa Prodanovic
% Copyright: (c) 2009 The University of Texas at Austin. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function grid = readGridFromBinaryFile(filename)
%
% Arguments
% filename - file name

% Returns
% grid      - Grid structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function grid = readGridFromBinaryFile(filename)

[zip,filename1] = checkUnzipFile(filename);

fid = fopen(filename1,'r');

grid.num_dims = fread(fid,1,'int');

grid.x_lo = fread(fid,3,'float');
grid.x_hi = fread(fid,3,'float');

grid.x_lo_ghostbox = fread(fid,3,'float');
grid.x_hi_ghostbox = fread(fid,3,'float');

grid.grid_dims = fread(fid,3,'int');
grid.grid_dims_ghostbox = fread(fid,3,'int');

grid.dx = fread(fid,3,'float');
grid.num_gridpts = fread(fid,1,'int');

% note that the rest of the LSMLIB grid structure is not read, 
% as is typically not relevant to basic visualization (fill boxes etc)

fclose(fid);

if(zip) 
    gzip(filename1); 
    delete(filename1); 
end

