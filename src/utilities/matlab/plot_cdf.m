function plot_cdf(filename)

data = readDataArray(filename);

data(:,:,[1:3,end-2:end])=[];
data(:,[1:3,end-2:end],:)=[];
data([1:3,end-2:end],:,:)=[];

nw_data = data(data < 0);

pd = fitdist(nw_data,'Normal');
x_values = min(nw_data):0.01:max(nw_data);
y = cdf(pd,x_values);
plot(x_values,y,'LineWidth',2)
hold on;
ecdf(nw_data);

end