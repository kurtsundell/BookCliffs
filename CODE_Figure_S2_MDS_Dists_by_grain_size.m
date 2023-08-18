% Supporting code for Sundell and Macdonald (submitted, EPSL) for bivariate KDE and two-dimensional
% similarity comparisons of global eHfT data from Puetz et al. (2021)

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

% Range of 2D data
xmin = 0;
xmax = 2500;
xint = 1;

% 1D kernel bandwidth
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
	
	nQ2(i,1) = floor(n(i,1)/2);
	nQ1(i,1) = floor(nQ2(i,1)/2);
	nQ3(i,1) = nQ1(i,1)+nQ2(i,1);
	
	data_tmp2 = sortrows(data_tmp2,2);
	
	All(i).data = data_tmp2;
	
	All(i).dataQ01 = data_tmp2(1:nQ1(i,1),:);
	All(i).dataQ12 = data_tmp2(nQ1(i,1)+1:nQ2(i,1),:);
	All(i).dataQ23 = data_tmp2(nQ2(i,1)+1:nQ3(i,1),:);
	All(i).dataQ34 = data_tmp2(nQ3(i,1)+1:end,:);

	All(i).data_R100 = datasample(All(i).data,100);
	All(i).data_R300 = datasample(All(i).data,300);
	
	All(i).data_S100 = data_tmp2(1:100,:);
	All(i).data_S300 = data_tmp2(1:300,:);

	All(i).data_L100 = data_tmp2(end-99:end,:);
	All(i).data_L300 = data_tmp2(end-299:end,:);
		
	KDE_All(:,i) = kde1(cell2num(All(i).data(:,1)), kernel*ones(length(All(i).data)), xmin, xmax, xint);
	KDE_Q01(:,i) = kde1(cell2num(All(i).dataQ01(:,1)), kernel*ones(length(All(i).dataQ01)), xmin, xmax, xint);	
	KDE_Q12(:,i) = kde1(cell2num(All(i).dataQ12(:,1)), kernel*ones(length(All(i).dataQ12)), xmin, xmax, xint);	
	KDE_Q23(:,i) = kde1(cell2num(All(i).dataQ23(:,1)), kernel*ones(length(All(i).dataQ23)), xmin, xmax, xint);	
	KDE_Q34(:,i) = kde1(cell2num(All(i).dataQ34(:,1)), kernel*ones(length(All(i).dataQ34)), xmin, xmax, xint);
	
	KDE_R300(:,i) = kde1(cell2num(All(i).data_R300(:,1)), kernel*ones(length(All(i).data_R300)), xmin, xmax, xint);
	KDE_R100(:,i) = kde1(cell2num(All(i).data_R100(:,1)), kernel*ones(length(All(i).data_R100)), xmin, xmax, xint);
		
	KDE_S300(:,i) = kde1(cell2num(All(i).data_S300(:,1)), kernel*ones(length(All(i).data_S300)), xmin, xmax, xint);
	KDE_S100(:,i) = kde1(cell2num(All(i).data_S100(:,1)), kernel*ones(length(All(i).data_S100)), xmin, xmax, xint);
	
	KDE_L300(:,i) = kde1(cell2num(All(i).data_L300(:,1)), kernel*ones(length(All(i).data_L300)), xmin, xmax, xint);
	KDE_L100(:,i) = kde1(cell2num(All(i).data_L100(:,1)), kernel*ones(length(All(i).data_L100)), xmin, xmax, xint);
end

%MDS 

colorst = jet(7);
colorst(2,:) = [];
colors(1,1:3) = [0 0 0];
colors(2:7,:) = colorst;

symbols = [{'p'},{'<'},{'v'},{'o'},{'s'},{'^'},{'>'}];
size = [3500,500,500,500,500,500,500];

for i = 1:N
	KDE_ALL = [KDE_All(:,i),KDE_Q01(:,i),KDE_Q12(:,i),KDE_Q23(:,i),KDE_Q34(:,i),KDE_S100(:,i),KDE_L100(:,i)];
	Name_All = [{'All'}, {'S100'}, {'<Q1'}, {'Q1-Q2'}, {'Q2-Q3'}, {'Q3-Q4'}, {'L100'}]';
	count = 1;
	for j = 1:length(Name_All)
		for k = 1:length(Name_All)
			name_comp(count,1) = strcat(Name_All(j,1), {' vs '}, Name_All(k,1));
			Rsquared(j,k) = ((sum((KDE_ALL(:,j) - mean(KDE_ALL(:,j))).*(KDE_ALL(:,k) - mean(KDE_ALL(:,k)))))/(sqrt((sum((KDE_ALL(:,j) - mean(KDE_ALL(:,j)))...
				.*(KDE_ALL(:,j) - mean(KDE_ALL(:,j)))))*(sum((KDE_ALL(:,k) - mean(KDE_ALL(:,k))).*(KDE_ALL(:,k) - mean(KDE_ALL(:,k))))))))^2;
			count = count+1;
		end
	end
	R2_out(:,:,i) = Rsquared;
	[XY,stress,disparities] = mdscale_new(1-Rsquared,2,'Criterion','stress'); % 2D MDS plot and metric squared stress only
	X=XY(:,1);
	Y=XY(:,2);
	dx = .02;
	figure
	hold on
	for m = 1:length(Name_All)
		scatter(X(m),Y(m),size(1,m), char(symbols(1,m)), 'MarkerFaceColor',colors(m,:),'MarkerEdgeColor','black', 'linewidth', 2);
	end
	[rubbish,md] = sort(Rsquared,1,'descend');
	text(X+dx,Y+dx,Name_All, 'FontSize', 16, 'Interpreter', 'none', 'FontWeight', 'Bold');
	set(gca,'FontSize',16)
	xlabel('Dimension I')
	ylabel('Dimension II')
	xlim([-.45 0.52])
	ylim([-.35 0.4])
	title(Name(i))

end


