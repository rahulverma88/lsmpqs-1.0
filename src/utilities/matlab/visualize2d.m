clc; clear; close all;

grid = readGridFromBinaryFile('grid.gz');
[mask,nx,ny,nz] = readDataArray('mask.gz');
%[data_step,~,~,~] = readDataArray('data_init.gz');
%[normal_vel,~,~,~] = readDataArray('normal_vel.gz');
%[curv_coeff,~,~,~] = readDataArray('curvature_coeff.gz');
%[vel_x,~,~,~] = readDataArray('external_velocity_x.gz');

% Determine file name to be saved
saveFile_pre = 'drain2d_';
readFile = 'data_step_1.gz';
saveFile_post = strrep(readFile,'data_','');
saveFile_post = strrep(saveFile_post,'.gz','');
saveFile = [saveFile_pre saveFile_post];

[data_step,~,~,~] = readDataArray(readFile);

%[mask_w,~,~,~] = readDataArray('mask_w.gz');
%[mask_nw,~,~,~] = readDataArray('mask_nw.gz');

% For trapped blob masks, uncomment this line
%data_step = -1*data_step;

data_step = max(mask,data_step);
g.dim = grid.num_dims;
g.dx = grid.dx(1);
g.min = [grid.x_lo_ghostbox(1); grid.x_lo_ghostbox(2)];


% There is a mismatch between lsmpqs and Matlab toolbox
% If there is an error in grid xs, reduce g.max by g.dx for both dims 
g.max = [grid.x_hi_ghostbox(1)-g.dx; grid.x_hi_ghostbox(2)-g.dx];4
g = processGrid(g);

% Increasing x0 moves window right
x0 = 1;
% Increasing y0 moves window up
y0 = 1;
x0_width = 55;
y0_width = 30;

bluemap = [0, 0, 0.3
    0, 0, 0.4
    0, 0, 0.5
    0, 0, 0.6];
 bluemap=   [0, 0, 0.9
    0, 0, 1.0];

if (size(g.xs{1}) == size(mask))
    %subplot(1,2,1)
    [C, hf] = contourf(g.xs{1}(x0:x0+x0_width,y0:y0+y0_width), g.xs{2}(x0:x0+x0_width,y0:y0+y0_width),...
        mask(x0:x0+x0_width,y0:y0+y0_width), [0 0], 'k-'); colormap gray
    set(hf, 'LineColor','none');
    %hold on;
    axis equal
    axis off;
    %saveas(gcf, 'testMask.tif');
    %fig = gcf;
    %fig.PaperPositionMode = 'auto';
    print('testMask','-dpng','-r300')

    %export_fig 'testMask.png' -m4
    
    figure;
    [~, hf_1] = contourf(g.xs{1}(x0:x0+x0_width,y0:y0+y0_width), g.xs{2}(x0:x0+x0_width,y0:y0+y0_width),...
        -data_step(x0:x0+x0_width,y0:y0+y0_width), [0 0], '-'); colormap hot;
    
    set(hf_1, 'LineColor','none');
    
    %set(hf, 'LineWidth', 1.0);
    %set(hf_1, 'LineWidth', 1.0);

    axis equal
    axis off
    %axis(g.axis);
    
    %saveas(gcf, 'testFig.tif');
    %fig = gcf;
    %fig.PaperPositionMode = 'auto';
    print('testFig','-dpng','-r300')
    %export_fig 'testFig.png' -m4
    %title('new')
 %{
    subplot(1,2,2)
     [C, hf] = contour(g.xs{1}(x0:x0+width,y0:y0+width), g.xs{2}(x0:x0+width,y0:y0+width),...
        mask(x0:x0+width,y0:y0+width), [0 0], 'k-'); %colormap gray
    hold on;
    [~, hf_1] = contourf(g.xs{1}(x0:x0+width,y0:y0+width), g.xs{2}(x0:x0+width,y0:y0+width),...
        -data_step_old(x0:x0+width,y0:y0+width), [0 0], '-'); colormap hot;
    set(hf, 'LineWidth', 1.0);
    set(hf_1, 'LineWidth', 1.0);

    axis equal
    %axis(g.axis);
    title('old')
     %}
else
    disp('adjust grid for processGrid');
end

close all;
bg = imread('testMask.png');
im = imread('testFig.png');
slice = im(:,:,3);
imAlphaData = ones(size(slice));
imAlphaData(im(:,:,3) >= 100) = 0;
imAlphaData(im(:,:,3) == 0) = 1;

figure;
ibg2 = image(bg);
axis off
axis equal
hold on
%   Overlay the image, and set the transparency previously calculated
iim2 = image(im);
set(iim2,'AlphaData',imAlphaData);
axis equal
axis off

%saveas(gcf, 'overlayImage.png');
print(saveFile,'-dpng','-r300');
%}