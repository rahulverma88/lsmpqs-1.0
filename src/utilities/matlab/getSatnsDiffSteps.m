grid = readGridFromBinaryFile('grid.gz');

mask = readDataArray('mask.gz');
mask = prune_data(mask);

g.dim = 3;
g.dx = grid.dx(1);
g.min = [grid.x_lo(1); grid.x_lo(2); grid.x_lo(3)];
g.max = [grid.x_hi(1)-grid.dx(1); grid.x_hi(2)-grid.dx(1); grid.x_hi(3)-grid.dx(1)];
g = processGrid(g);


npairs = size(pairs,1);

layer = 0.5*grid.dx(1);

for i = 1:npairs
    fname_1 = ['data_step_' num2str(pairs(i,1)) '.gz'];
    data_1 = readDataArray(fname_1);
    data_1 = prune_data(data_1);
        
    fname_2 = ['data_step_' num2str(pairs(i,2)) '.gz'];
    data_2 = readDataArray(fname_2);
    data_2 = prune_data(data_2);
    
    data_2 = data_2 - layer;
    
    % make data_1 a mask
    data_1 = -data_1;
    
    data_diff = max(data_1, data_2);
    
    satn = size(find(data_diff < 0),1)/size(find(mask < 0),1);
    disp(satn);
    
end