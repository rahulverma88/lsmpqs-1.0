function k_avg = getLSMCurvs(step_start, step_end)

grid = readGridFromBinaryFile('grid.gz');
[mask,~,~,~] = readDataArray('mask.gz');

% remove boundaries from mask
mask = prune_data(mask);

g.dim = grid.num_dims;
g.dx = grid.dx(1);
g.min = [grid.x_lo(1); grid.x_lo(2); grid.x_lo(3)];

% There is a mismatch between lsmpqs and Matlab toolbox
% If there is an error in grid xs, reduce g.max by g.dx for both dims 
g.max = [grid.x_hi(1); grid.x_hi(2); grid.x_hi(3)];
g = processGrid(g);

nsteps = step_end - step_start;
k_avg = zeros(nsteps,1);
Pc = k_avg;

sigma = 0.072;
dx_sim = g.dx(1);
dx_real = 5.345*1e-6;
pa_to_psig = 0.000145;

for i = step_start:step_end
    fname = ['data_step' num2str(i) '.gz'];
    [data_step,~,~,~] = readDataArray(fname);
    data_step = prune_data(data_step);
    [k_avg(i - step_start + 1), ~, ~] = findCurvatureNearIntfc(data_step,g,mask);
    Pc(i) = k_avg(i - step_start + 1)*sigma*dx_sim/dx_real*pa_to_psig;
end




%}