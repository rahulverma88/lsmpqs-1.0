clear; close all; clc;
[data_step, ~, ~, ~] = readDataArray('data_step_2.gz');
[mask,nx,ny,nz] = readDataArray('mask.gz');
im = data_step;
im(data_step < 0) = 2;
im(data_step >= 0) = 0;
im_mask = mask;
im_mask(mask >= 0) = 4;
im_mask(mask < 0) = 0;

%cmap = gray;
cmap(1,:) = [1 1 1];
cmap(2,:) = [1 0 0];
%cmap(3,:) = [

im_rgb = ind2rgb(im, cmap);
%im_mask_rgb = ind2rgb(im_mask, cmap);

hfig = figure;
hIm = imshow(im_rgb,[]);
hSp = imscrollpanel(hfig, hIm);
%axis equal;
%imshow(data_step(2000:3000,2000:3000))