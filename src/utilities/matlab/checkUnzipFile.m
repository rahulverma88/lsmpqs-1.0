%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:    Masa Prodanovic
% Copyright: (c) 2009 The University of Texas at Austin. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [zip,filename1] = checkUnzipFile(filename)
%
%  Check whether the filename has .gz extension and
%  unzip the file as necessary.
%
%  Arguments
%  filename - file name to be checked
%  
%  Returns
%  zip       - 1 if the extension was found, 0 otherwise
%  filename1 - name of the unzipped file (same as filename if
%              the file was not zipped to begin with)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [zip,filename1] = checkUnzipFile(filename)
    if( findstr('.gz',filename) )
        zip = 1;
        filename = gunzip(filename);
        filename1 = filename{1};
    else
        %checks if a zipped version of the filename exists
        filename_zip = sprintf('%s.gz',filename);
    
        if(exist(filename_zip,'file'))
	   tmp_fname = gunzip(filename_zip);
	   filename1 = tmp_fname{1};
	   zip = 1;
        else   
	   filename1 = filename;
	   zip = 0;
        end  
    end
    
return
