clc; clear; close all;

grid = readGridFromBinaryFile('grid.gz');
[mask,nx,ny,nz] = readDataArray('mask.gz');
saveFile_pre = 'drain2d_';
readFile = 'data_init.gz';
saveFile_post = strrep(readFile,'data_','');
saveFile_post = strrep(saveFile_post,'.gz','');
saveFile = [saveFile_pre saveFile_post];

[data_step,~,~,~] = readDataArray(readFile);

nw_color = [1    0.270    0]; % orange
mask_color = [0.411 0.411 0.411]; % gray

data_plot_1 = data_step; data_plot_2 = data_step; data_plot_3 = data_step;

% White (blank) for w phase everything
data_plot_1(data_step >= 0) = 1;
data_plot_2(data_step >= 0) = 1;
data_plot_3(data_step >= 0) = 1;

% Color in mask 
data_plot_1(mask >= 0) = mask_color(1); 
data_plot_2(mask >= 0) = mask_color(2); 
data_plot_3(mask >= 0) = mask_color(3);

% Color in nw phase
data_plot_1(data_step <= 0) = nw_color(1); 
data_plot_2(data_step <= 0) = nw_color(2); 
data_plot_3(data_step <= 0) = nw_color(3);

%concatenate everything into a single colored image and image it
data_plot = cat(3, data_plot_1, data_plot_2, data_plot_3);

h = image(rot90(data_plot));

axis equal;
axis off;
%{
dstep_bin = data_step;
dstep_bin(data_step < 0) = 1;
dstep_bin(data_step > 0) = 0;
imshow(dstep_bin);
%}





