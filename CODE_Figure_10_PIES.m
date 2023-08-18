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
extra = 100;

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

	data_tmp2 = sortrows(data_tmp2,2);
	All(i).data = data_tmp2;
	All(i).data_S100 = data_tmp2(1:100,:);
	All(i).data_L100 = data_tmp2(end-99:end,:);	
	KDE_S100(:,i) = kde1(cell2num(All(i).data_S100(:,1)), kernel*ones(length(All(i).data_S100)), xmin, xmax, xint);
	KDE_L100(:,i) = kde1(cell2num(All(i).data_L100(:,1)), kernel*ones(length(All(i).data_L100)), xmin, xmax, xint);

	
% 	Archean ages as 4000 – 2500 Ma, 
% 	Trans-Hudson ages as 2000 – 1800 Ma, 
% 	Yavapai – Mazatzal ages as 1800 – 1600 Ma, 
% 	Midcontinent ages as 1550 – 1300 Ma, 
% 	Grenville ages as 1250 – 900 Ma, 
% 	Appalachian - Peri-Gondwana ages as 750 – 300 Ma, 
% 	Western Cordillera ages as < 300 Ma

	
	
	Pies_L100(i).WC = All(i).data_L100(cell2num(All(i).data_L100(:,1)) <= 300);
	Pies_L100(i).AppPG = All(i).data_L100(300 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 750);
	Pies_L100(i).Gren = All(i).data_L100(900 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 1200);
	Pies_L100(i).MC = All(i).data_L100(1300 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 1550);
	Pies_L100(i).YM = All(i).data_L100(1600 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 1800);
	Pies_L100(i).TH = All(i).data_L100(1800 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 2000);
	Pies_L100(i).Arch = All(i).data_L100(2400 < cell2num(All(i).data_L100(:,1)));
	
	Pies_S100(i).WC = All(i).data_S100(cell2num(All(i).data_S100(:,1)) <= 300);
	Pies_S100(i).AppPG = All(i).data_S100(300 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 750);
	Pies_S100(i).Gren = All(i).data_S100(900 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 1200);
	Pies_S100(i).MC = All(i).data_S100(1300 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 1550);
	Pies_S100(i).YM = All(i).data_S100(1600 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 1800);
	Pies_S100(i).TH = All(i).data_S100(1800 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 2000);
	Pies_S100(i).Arch = All(i).data_S100(2400 < cell2num(All(i).data_S100(:,1)));	
	
		

	Pie_L100(i,:) = [length(Pies_L100(i).WC), length(Pies_L100(i).AppPG), length(Pies_L100(i).Gren), length(Pies_L100(i).MC), length(Pies_L100(i).YM), length(Pies_L100(i).TH), ...
	length(Pies_L100(i).Arch)];
	
	Pie_S100(i,:) = [length(Pies_S100(i).WC), length(Pies_S100(i).AppPG), length(Pies_S100(i).Gren), length(Pies_S100(i).MC), length(Pies_S100(i).YM), length(Pies_S100(i).TH), ...
	length(Pies_S100(i).Arch)];




end



for i = 1:N
	base = 0;
	figure
	hold on
	
	patch([x;xmax;0], [KDE_S100(:,i)+base;base;base],'k')
	plot(x, KDE_S100(:,i)+base,'color','k','linewidth',1)	
	text(xmax,base+max(KDE_S100(:,i))/2,'Smallest 100','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_S100(:,i));
	
	patch([x;xmax;0], [KDE_L100(:,i)+base;base;base],'k')
	plot(x, KDE_L100(:,i)+base,'color','k','linewidth',2)	
	text(xmax,base+max(KDE_L100(:,i))/2,'Largest 100','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_L100(:,i));
	
	title(Name(i))
	
	ylim([0 base])
end






labels = {'WC', 'AppPG', 'Gren', 'MC', 'YM', 'TH', 'Arch'};





for i = 1:N
	figure

	subplot(1,2,1);
	pie(Pie_S100(i,:),'%.0f%%')
	subplot(1,2,2);
	pie(Pie_L100(i,:),'%.0f%%')
		title(Name(i))
end

