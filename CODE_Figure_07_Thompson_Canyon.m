clear all
close all
clc

rng('default')

[filename pathname] = uigetfile({'*'},'File Selector'); %load the supplemental file with zircon age eHfT data

if ispc == 1
	fullpathname = char(strcat(pathname, '\', filename));
end
if ismac == 1
	fullpathname = char(strcat(pathname, '/', filename));
end

% x range 
xmin = 0;
xmax = 2500;
xint = 1;

% kernel bandwidth
kernel = 15;

x = [xmin:xint:xmax]';

% Read in data, format is name header and two columns of info, for our example we use age + Hf, but any 2D data will work
[numbers text1, data] = xlsread(fullpathname);
numbers = num2cell(numbers);

% Filter out any data that are not pairs of numbers
for i = 1:size(numbers,1)
	for j = 1:size(numbers,2)
		if cellfun('isempty', numbers(i,j)) == 0
			if cellfun(@isnan, numbers(i,j)) == 1
				numbers(i,j) = {[]};
			end
		end
	end
end

% pull the names from the headers
for i = 1:(size(data,2)+1)/2
	Name(i,1) = data(1,i*2-1);
end

data_tmp = numbers(1:end,:); %use temporary variable
N = size(data_tmp,2)/2; % figure out how many samples

for i = 1:N
	data_tmp2 = numbers(:,i*2-1:i*2);
	data_tmp2 = data_tmp2(~any(cellfun(@isempty, data_tmp2),2), :);
	n(i,1) = length(data_tmp2(:,1));
	
	All(i).data = data_tmp2;
	All(i).data_R100 = datasample(All(i).data,100);
	All(i).data_R300 = datasample(All(i).data,300);	
	
	KDE_All(:,i) = kde1(cell2num(All(i).data(:,1)), kernel*ones(length(All(i).data)), xmin, xmax, xint);
	KDE_R300(:,i) = kde1(cell2num(All(i).data_R300(:,1)), kernel*ones(length(All(i).data_R300)), xmin, xmax, xint);
	KDE_R100(:,i) = kde1(cell2num(All(i).data_R100(:,1)), kernel*ones(length(All(i).data_R100)), xmin, xmax, xint);
end

figure
hold on
patch([x;xmax;0], [KDE_All(:,1);0;0],[166/255 124/255 82/255])
plot(x, KDE_All(:,1),'k','linewidth',2)
title(Name(1))

figure
hold on
patch([x;xmax;0], [KDE_All(:,2);0;0],[0 1 0])
plot(x, KDE_All(:,2),'k','linewidth',2)
title(Name(2))

figure
hold on
patch([x;xmax;0], [KDE_All(:,3);0;0],[1 1 0])
plot(x, KDE_All(:,3),'k','linewidth',2)
title(Name(3))
	
figure
hold on
patch([x;xmax;0], [KDE_All(:,4);0;0],[1 1 0])
plot(x, KDE_All(:,4),'k','linewidth',2)
title(Name(4))

figure
hold on
patch([x;xmax;0], [KDE_All(:,5);0;0],[166/255 124/255 82/255])
plot(x, KDE_All(:,5),'k','linewidth',2)
title(Name(5))

count = 1;
for j = 1:length(KDE_All(1,:))
	for k = 1:length(KDE_All(1,:))
		name_comp(count,1) = strcat(Name(j,1), {' vs '}, Name(k,1));
		Rsquared(j,k) = ((sum((KDE_All(:,j) - mean(KDE_All(:,j))).*(KDE_All(:,k) - mean(KDE_All(:,k)))))/(sqrt((sum((KDE_All(:,j) - mean(KDE_All(:,j)))...
			.*(KDE_All(:,j) - mean(KDE_All(:,j)))))*(sum((KDE_All(:,k) - mean(KDE_All(:,k))).*(KDE_All(:,k) - mean(KDE_All(:,k))))))))^2;
		count = count+1;
	end
end

% Multidimensional Scaling (MDS) Figure 6
symbols = [{'s'},{'s'},{'s'},{'o'},{'o'},    {'s'},{'s'},{'s'},{'o'},{'o'},    {'s'},{'s'},{'s'},{'o'},{'o'}];
size = 900;
colors = [[166/255 124/255 82/255];[0 1 0];[1 1 0];[1 1 0];[166/255 124/255 82/255];     [166/255 124/255 82/255];[0 1 0];[1 1 0];[1 1 0];[166/255 124/255 82/255];...
	[166/255 124/255 82/255];[0 1 0];[1 1 0];[1 1 0];[166/255 124/255 82/255]];  
	
KDE_ALL = [KDE_All,KDE_R300,KDE_R100];

Name_All = [Name; strcat(Name,{' '},'n=300'); strcat(Name,{' '},'n=100')];

count = 1;
for j = 1:length(Name_All)
	for k = 1:length(Name_All)
		name_comp(count,1) = strcat(Name_All(j,1), {' vs '}, Name_All(k,1));
		Rsquared(j,k) = ((sum((KDE_ALL(:,j) - mean(KDE_ALL(:,j))).*(KDE_ALL(:,k) - mean(KDE_ALL(:,k)))))/(sqrt((sum((KDE_ALL(:,j) - mean(KDE_ALL(:,j)))...
			.*(KDE_ALL(:,j) - mean(KDE_ALL(:,j)))))*(sum((KDE_ALL(:,k) - mean(KDE_ALL(:,k))).*(KDE_ALL(:,k) - mean(KDE_ALL(:,k))))))))^2;
		count = count+1;
	end
end
[XY,stress,disparities] = mdscale_new(1-Rsquared,2,'Criterion','stress'); % 2D MDS plot and metric squared stress only
X=XY(:,1);
Y=XY(:,2);
dx = .01;
figure
hold on
for m = 1:length(Name_All)
	scatter(X(m),Y(m),size, 's', 'MarkerFaceColor',colors(m,:),'MarkerEdgeColor','black', 'linewidth', 2);
end
[rubbish,md] = sort(Rsquared,1,'descend');
text(X+dx,Y+dx,Name_All, 'FontSize', 16, 'Interpreter', 'none', 'FontWeight', 'Bold');
set(gca,'FontSize',16)
xlabel('Dimension I')
ylabel('Dimension II')
%axis([-0.2 0.6 -0.15 0.35])
%xlim([-.2 0.6])
%ylim([-.2 0.4])




for i = 1:N
	figure
	
	subplot(3,1,3);
	plot(x, KDE_All(:,i),'k','linewidth',2)
	
	subplot(3,1,2);
	plot(x, KDE_R300(:,i),'k','linewidth',2)
	
	subplot(3,1,1);
	plot(x, KDE_R100(:,i),'k','linewidth',2)
	title(Name(i))
end




% MDAs










