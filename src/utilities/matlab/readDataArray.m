%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:    Masa Prodanovic
% Copyright: (c) 2009 The University of Texas at Austin. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [data,nx,ny,nz] = readDataArray(filename)
%
% Reads a binary data file which contains 3 binary integers nx, ny, nz 
% indicating sizes of the array in each dimension, followed by an array
% of nx*ny*nz single precision real numbers.
% 
% Arguments:
% filename - file name
%
% Returns
% data      - data array read from the file
% nx,ny,nz  - array dimensions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,nx,ny,nz] = readDataArray(filename)

data = [];
nx = 0; ny = 1; nz = 1;

[zip,filename1] = checkUnzipFile(filename);

if( exist(filename1,'file') )
    fid = fopen(filename1,'r');
    nx = fread(fid,1,'int');
    ny = fread(fid,1,'int');
    nz = fread(fid,1,'int');

    data = fread(fid,nx*ny*nz,'float');
    data = reshape(data,[nx ny nz]);
    fclose(fid);

    if(zip) 
	gzip(filename1); delete(filename1);
    end;
else
    data = [];
    nx = 0; ny = 0; nz = 0;
end

return
