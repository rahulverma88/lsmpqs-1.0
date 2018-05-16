clc; clear; close all;
data = readDataArray('data_step_10.gz');
data_bin = data;
data_bin(data < 0) = 1;
data_bin(data > 0) = 0;
imshow(data_bin);