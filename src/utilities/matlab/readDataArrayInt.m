
function [data,nx,ny,nz] = readDataArrayInt(filename)

nx = 500; ny = 500; nz = 500;

[zip,filename1] = checkUnzipFile(filename);

if( exist(filename1,'file') )
    fid = fopen(filename1,'r');
    nx = fread(fid,1,'int');
    ny = fread(fid,1,'int');
    nz = fread(fid,1,'int');

    data = fread(fid,nx*ny*nz,'int');
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
