function [data] = readDataArrayUchar(filename)

nx = 799; ny = 799; nz = 799;

[zip,filename1] = checkUnzipFile(filename);

if( exist(filename1,'file') )
    fid = fopen(filename1,'r');
    %nx = fread(fid,1,'int32');
    %ny = fread(fid,1,'int32');
    %nz = fread(fid,1,'int32');

    data = fread(fid,nx*ny*nz,'unsigned char');
    data = reshape(data,[nx ny nz]);
    fclose(fid);

    if(zip) 
	gzip(filename1); delete(filename1);
    end;
else
    data = [];
end

return