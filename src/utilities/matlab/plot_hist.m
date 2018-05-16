function plot_hist(filename)

data = readDataArray(filename);

data(:,:,[1:3,end-2:end])=[];
data(:,[1:3,end-2:end],:)=[];
data([1:3,end-2:end],:,:)=[];

nw_data = data(data < 0);

histfit(nw_data);
xlabel('LS values');
ylabel('Frequency');

end