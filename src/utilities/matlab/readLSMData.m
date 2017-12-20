clc; clear; close all;

grid = readGridFromBinaryFile('grid.gz');
[mask,nx,ny,nz] = readDataArray('mask.gz');
g.dim = grid.num_dims;
g.dx = grid.dx(1);
g.min = [grid.x_lo_ghostbox(1); grid.x_lo_ghostbox(2); grid.x_lo_ghostbox(3)];

% There is a mismatch between lsmpqs and Matlab toolbox
% If there is an error in grid xs, reduce g.max by g.dx for both dims 
g.max = [grid.x_hi_ghostbox(1)-g.dx; grid.x_hi_ghostbox(2); grid.x_hi_ghostbox(3)-g.dx];
g = processGrid(g);

nsteps = 27;
k_avg = zeros(nsteps,1);
Pc = k_avg;

sigma = 0.072;
dx_sim = g.dx(1);
dx_real = 5.345*1e-6;
pa_to_psig = 0.000145;

for i = 1:nsteps
    fname = ['data_step_' num2str(i) '.gz'];
    [data_step,~,~,~] = readDataArray(fname);
    [k_avg(i), ~, ~] = findCurvatureNearIntfc(data_step,g,mask);
    Pc(i) = k_avg(i)*sigma*dx_sim/dx_real*pa_to_psig;
end




%}