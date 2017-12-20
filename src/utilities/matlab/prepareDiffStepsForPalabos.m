%{
pairs = [ [7,18];
          [11,13];
          [11,14];
          [12,15];
          [13,15];
          [14,15];
          [12,16];
          [12,17];
          [13,17];
          [14,17];
          [14,18];
          [15,18];
          [16,18]];
%  
[8,9];
          [9,10];
          [10,11];
          [11,12];
          [13,14];
          [17,18]];
%}
npairs = size(pairs,1);
interp = 1;


for i = 1:npairs
    fname_1 = ['data_step_' num2str(pairs(i,1)) '.gz'];
    data_1 = readDataArray(fname_1);
    data_1 = prune_data(data_1);
    
    [nx, ny, nz] =  size(data_1);
    data_rot_1 = zeros(nz,ny,nx);
    for j=1:nz
        data_rot_1(:,:,j) = reshape(data_1(j,:,:),ny,nz);
    end
    
    data_1 = data_rot_1;    
    
    data_1 = data_1 + 0.5*0.04;
    
    fname_2 = ['data_step_' num2str(pairs(i,2)) '.gz'];
    data_2 = readDataArray(fname_2);
    data_2 = prune_data(data_2);
    
    %data_2 = data_2 - 0.5*0.04;
    
    [nx, ny, nz] =  size(data_2);
    data_rot_2 = zeros(nz,ny,nx);
    for j=1:nz
        data_rot_2(:,:,j) = reshape(data_2(j,:,:),ny,nz);
    end
    
    data_2 = data_rot_2;
    
    % make data_1 a mask
    data_1 = -data_1;
    
    data_diff = max(data_1, data_2);
    
    if (interp)
        data_interp = interp3(data_diff,'spline');
        data_interp_uchar = data_interp;
        data_interp_uchar(data_interp < 0) = 1;
        data_interp_uchar(data_interp >= 0) = 0;
        data_interp_comp = bwlabeln(data_interp_uchar);
        slice_start = data_interp_comp(:,:,1);
        slice_end = data_interp_comp(:,:,end);
        if(size(intersect(slice_start,slice_end),1) > 1)
            disp('Percolation in');
            disp(pairs(i,1)); disp(pairs(i,2));
            percolate = 1;
        else
            percolate = 0;
        end
        fname = ['diff_steps_interp2_' num2str(pairs(i,1)) '_' num2str(pairs(i,2)) '.dat'];
        if (percolate)
            createDATImage(data_interp_uchar, size(data_interp_uchar,3), fname);
        end
    else
        data_uchar = data_diff;
        data_uchar(data_diff < 0) = 1;
        data_uchar(data_diff >= 0) = 0;
        data_comp = bwlabeln(data_uchar);
        slice_start = data_comp(:,:,1);
        slice_end = data_comp(:,:,end);
        if(size(intersect(slice_start,slice_end),1) > 1)
            disp('Percolation in');
            disp(pairs(i,1)); disp(pairs(i,2));
            percolate = 1;
        else
            percolate = 0;
        end
        fname = ['diff_steps_' num2str(pairs(i,1)) '_' num2str(pairs(i,2)) '.dat'];
        if (percolate)
            createDATImage(data_uchar, size(data_uchar,3), fname);
        end
    end
end
