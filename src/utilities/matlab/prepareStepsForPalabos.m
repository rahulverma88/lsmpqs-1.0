function prepareStepsForPalabos(steps, phase, interp)

mask = readDataArray('mask.gz');
mask(1:3,:,:) = [];
mask(end-2:end,:,:) = [];
mask(:,1:3,:) = [];
mask(:,end-2:end,:) = [];
mask(:,:,1:3) = [];
mask(:,:,end-2:end) = [];
    
[nx, ny, nz] =  size(mask);
mask_rot = zeros(nz,ny,nx);

for j=1:nz
    mask_rot(:,:,j) = reshape(mask(j,:,:),ny,nz);
end

mask = mask_rot;

for i=steps
    fname = ['data_step_' num2str(i) '.gz'];
    data = readDataArray(fname);
    data(1:3,:,:) = [];
    data(end-2:end,:,:) = [];
    data(:,1:3,:) = [];
    data(:,end-2:end,:) = [];
    data(:,:,1:3) = [];
    data(:,:,end-2:end) = [];
    
    [nx, ny, nz] =  size(data);
    data_rot = zeros(nz,ny,nx);
    for j=1:nz
        data_rot(:,:,j) = reshape(data(j,:,:),ny,nz);
    end
    
    data = data_rot;    
    if (strcmp(phase,'nw'))
        if (interp)
            data_interp = interp3(data,'spline');
            data_interp_uchar = data_interp;
            data_interp_uchar(data_interp < 0) = 1;
            data_interp_uchar(data_interp >= 0) = 0;
            data_interp_comp = bwlabeln(data_interp_uchar);
            slice_start = data_interp_comp(:,:,1);
            slice_end = data_interp_comp(:,:,end);
            if(size(intersect(slice_start,slice_end),1) > 1)
                disp('Percolation in step');
                disp(i);
                percolate = 1;
            else
                percolate = 0;
            end
            fname = ['data_step_nw_' num2str(i) '_interp.dat'];
            if(percolate)
                createDATImage(data_interp_uchar, size(data_interp,3), fname);
            end
        else
            data_uchar = data;
            data_uchar(data < 0) = 1;
            data_uchar(data >= 0) = 0;
            fname = ['data_step_nw_' num2str(i) '.dat'];
            createDATImage(data_uchar, size(data_uchar,3), fname);
        end
    else
        data = -data;
        data(mask >= 0) = mask(mask >= 0);
        if (interp)
            data_interp = interp3(data,'spline');
            data_interp_uchar = data_interp;
            data_interp_uchar(data_interp < 0) = 1;
            data_interp_uchar(data_interp >= 0) = 0;
            data_interp_comp = bwlabeln(data_interp_uchar);
            slice_start = data_interp_comp(:,:,1);
            slice_end = data_interp_comp(:,:,end);
            if(size(intersect(slice_start,slice_end),1) > 1)
                disp('Percolation in step');
                disp(i);
                percolate = 1;
            else
                percolate = 0;
            end
            fname = ['data_step_w_' num2str(i) '_interp.dat'];
            if(percolate)
                createDATImage(data_interp_uchar, size(data_interp,3), fname);
            end
        else
            data_uchar = data;
            data_uchar(data < 0) = 1;
            data_uchar(data >= 0) = 0;
            fname = ['data_step_w_' num2str(i) '.dat'];
            createDATImage(data_uchar, size(data_uchar,3), fname);
        end
    end
end
