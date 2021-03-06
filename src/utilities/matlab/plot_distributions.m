clc;clear; %close all;
data = readDataArray('data_step_4.gz');

data(:,:,[1:3,end-2:end])=[];
data(:,[1:3,end-2:end],:)=[];
data([1:3,end-2:end],:,:)=[];

nw_data = data(data < 0);

figure;
subplot(2,1,1);
histfit(nw_data);

subplot(2,1,2);
pd = fitdist(nw_data,'Normal');
x_values = min(nw_data):0.01:max(nw_data);
y = cdf(pd,x_values);
plot(x_values,y,'LineWidth',2)
hold on;
ecdf(nw_data);