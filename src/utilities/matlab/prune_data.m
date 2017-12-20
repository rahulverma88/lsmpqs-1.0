function data_return = prune_data(data)

[r,c,p] = size(data);

data(1:3,:,:) = [];
data(:,1:3,:) = [];
data(:,:,1:3) = [];

data((r-5):(r-3),:,:) = [];
data(:,(c-5):(c-3),:) = [];
data(:,:,(p-5):(p-3)) = [];

data_return = data;