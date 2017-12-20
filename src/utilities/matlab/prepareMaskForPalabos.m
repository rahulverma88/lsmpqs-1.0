interp = 1;
[mask, ~, ~, ~] = readDataArray('mask.gz');
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

dx = 0.04;
eps = 0*dx;

if (interp)
    mask_interp = interp3(mask,'spline');%,'spline');
    mask_interp_uchar = mask_interp;
    mask_interp_uchar(mask_interp < -eps) = 1;
    mask_interp_uchar(mask_interp >= -eps) = 0;
    createDATImage(mask_interp_uchar, size(mask_interp,3), 'mask_interp.dat');
else
    mask_uchar = mask;
    mask_uchar(mask < -eps) = 1;
    mask_uchar(mask >= -eps) = 0;
    createDATImage(mask_uchar, size(mask,3), 'mask.dat');
end
