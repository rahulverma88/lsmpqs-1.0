function writeDataArrayUchar(data,filename)

fid = fopen(filename,'w');
fwrite(fid,data);
fclose(fid);


return