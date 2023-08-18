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

	Pies(i).WC = All(i).data(cell2num(All(i).data(:,1)) <= 275);
	Pies(i).App = All(i).data(275 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 500);
	Pies(i).PG = All(i).data(500 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 700);
	Pies(i).Gren = All(i).data(900 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 1250);
	Pies(i).MC = All(i).data(1300 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 1550);
	Pies(i).YM = All(i).data(1600 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 1800);
	Pies(i).TH = All(i).data(1800 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 2000);
	Pies(i).Arch = All(i).data(2500 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 4000);
	
	
	Pies(i).Other = All(i).data(700 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 900 |...
		1250 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 1300 |...
		1550 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 1600 |...
		2000 < cell2num(All(i).data(:,1)) & cell2num(All(i).data(:,1)) <= 2500);
	
	
% 	Archean ages as 4000 – 2500 Ma, 
% 	Trans-Hudson ages as 2000 – 1800 Ma, 
% 	Yavapai – Mazatzal ages as 1800 – 1600 Ma, 
% 	Midcontinent ages as 1550 – 1300 Ma, 
% 	Grenville ages as 1250 – 900 Ma, 
% 	Peri-Gondwana ages as 700 – 500 Ma, 
% 	Appalachian ages as 500 – 275 Ma, and 
% 	Western Cordillera ages as < 275 Ma
	
	
	
	
	PiesQ01(i).WC = All(i).dataQ01(cell2num(All(i).dataQ01(:,1)) <= 275);
	PiesQ01(i).App = All(i).dataQ01(275 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 500);
	PiesQ01(i).PG = All(i).dataQ01(500 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 730);
	PiesQ01(i).Gren = All(i).dataQ01(930 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 1250);
	PiesQ01(i).MC = All(i).dataQ01(1300 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 1550);
	PiesQ01(i).YM = All(i).dataQ01(1600 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 1800);
	PiesQ01(i).TH = All(i).dataQ01(1800 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 2000);
	PiesQ01(i).Arch = All(i).dataQ01(2350 < cell2num(All(i).dataQ01(:,1)));
	PiesQ01(i).Other = All(i).dataQ01(730 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 930 | ... 
		2000 < cell2num(All(i).dataQ01(:,1)) & cell2num(All(i).dataQ01(:,1)) <= 2350);
	
	PiesQ34(i).WC = All(i).dataQ34(cell2num(All(i).dataQ34(:,1)) <= 275);
	PiesQ34(i).App = All(i).dataQ34(275 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 500);
	PiesQ34(i).PG = All(i).dataQ34(500 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 730);
	PiesQ34(i).Gren = All(i).dataQ34(930 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 1250);
	PiesQ34(i).MC = All(i).dataQ34(1300 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 1550);
	PiesQ34(i).YM = All(i).dataQ34(1600 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 1800);
	PiesQ34(i).TH = All(i).dataQ34(1800 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 2000);
	PiesQ34(i).Arch = All(i).dataQ34(2350 < cell2num(All(i).dataQ34(:,1)));
	PiesQ34(i).Other = All(i).dataQ34(730 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 930 | ... 
		2000 < cell2num(All(i).dataQ34(:,1)) & cell2num(All(i).dataQ34(:,1)) <= 2350);	
	
	
	
	
	
	
	Pies_L100(i).WC = All(i).data_L100(cell2num(All(i).data_L100(:,1)) <= 275);
	Pies_L100(i).App = All(i).data_L100(275 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 500);
	Pies_L100(i).PG = All(i).data_L100(500 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 730);
	Pies_L100(i).Gren = All(i).data_L100(930 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 1250);
	Pies_L100(i).MC = All(i).data_L100(1300 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 1550);
	Pies_L100(i).YM = All(i).data_L100(1600 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 1800);
	Pies_L100(i).TH = All(i).data_L100(1800 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 2000);
	Pies_L100(i).Arch = All(i).data_L100(2350 < cell2num(All(i).data_L100(:,1)));
	Pies_L100(i).Other = All(i).data_L100(730 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 930 | ... 
		2000 < cell2num(All(i).data_L100(:,1)) & cell2num(All(i).data_L100(:,1)) <= 2350);
	
	
	
	
	
	
	
	Pies_S100(i).WC = All(i).data_S100(cell2num(All(i).data_S100(:,1)) <= 275);
	Pies_S100(i).App = All(i).data_S100(275 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 500);
	Pies_S100(i).PG = All(i).data_S100(500 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 730);
	Pies_S100(i).Gren = All(i).data_S100(930 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 1250);
	Pies_S100(i).MC = All(i).data_S100(1300 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 1550);
	Pies_S100(i).YM = All(i).data_S100(1600 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 1800);
	Pies_S100(i).TH = All(i).data_S100(1800 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 2000);
	Pies_S100(i).Arch = All(i).data_S100(2350 < cell2num(All(i).data_S100(:,1)));
	Pies_S100(i).Other = All(i).data_S100(730 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 930 | ... 
		2000 < cell2num(All(i).data_S100(:,1)) & cell2num(All(i).data_S100(:,1)) <= 2350);
	
	Pie_n(i,:) = [length(Pies(i).WC), length(Pies(i).App), length(Pies(i).PG), length(Pies(i).Gren), length(Pies(i).MC), length(Pies(i).YM), length(Pies(i).TH), ...
	length(Pies(i).Arch), length(Pies(i).Other)];
	
	Pie_nQ01(i,:) = [length(PiesQ01(i).WC), length(PiesQ01(i).App), length(PiesQ01(i).PG), length(PiesQ01(i).Gren), length(PiesQ01(i).MC), length(PiesQ01(i).YM), length(PiesQ01(i).TH), ...
	length(PiesQ01(i).Arch), length(PiesQ01(i).Other)];

	Pie_nQ34(i,:) = [length(PiesQ34(i).WC), length(PiesQ34(i).App), length(PiesQ34(i).PG), length(PiesQ34(i).Gren), length(PiesQ34(i).MC), length(PiesQ34(i).YM), length(PiesQ34(i).TH), ...
	length(PiesQ34(i).Arch), length(PiesQ34(i).Other)];
	
	Pie_L100(i,:) = [length(Pies_L100(i).WC), length(Pies_L100(i).App), length(Pies_L100(i).PG), length(Pies_L100(i).Gren), length(Pies_L100(i).MC), length(Pies_L100(i).YM), length(Pies_L100(i).TH), ...
	length(Pies_L100(i).Arch), length(Pies_L100(i).Other)];
	
	Pie_S100(i,:) = [length(Pies_S100(i).WC), length(Pies_S100(i).App), length(Pies_S100(i).PG), length(Pies_S100(i).Gren), length(Pies_S100(i).MC), length(Pies_S100(i).YM), length(Pies_S100(i).TH), ...
	length(Pies_S100(i).Arch), length(Pies_S100(i).Other)];




end


colorst = jet(7);
colorst(2,:) = [];
colors(1,1:3) = [0 0 0];
colors(2:7,:) = colorst;


for i = 1:N
	base = 0;
	figure
	hold on
	
	%colors = [0 0 1;0 0.5 1;0 1 1;0 0.5 0;1 0.5 0;1 0 0];
	patch([x;xmax;0], [KDE_All(:,i);0;0],colors(1,:))
	plot(x, KDE_All(:,i),'k','linewidth',2)
	base = max(KDE_All(:,i));
	text(xmax,base/2,'All','fontsize',16, 'horizontalAlignment', 'right')
	
	patch([x;xmax;0], [KDE_S100(:,i)+base;base;base],colors(2,:))
	plot(x, KDE_S100(:,i)+base,'color','k','linewidth',1)	
	text(xmax,base+max(KDE_S100(:,i))/2,'Smallest 100','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_S100(:,i));
	
	patch([x;xmax;0], [KDE_Q01(:,i)+base;base;base],colors(3,:))
	plot(x, KDE_Q01(:,i)+base,'color','k','linewidth',2)
	text(xmax,base+max(KDE_Q01(:,i))/2,'< Q1','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_Q01(:,i));	
	
	patch([x;xmax;0], [KDE_Q12(:,i)+base;base;base],colors(4,:))
	plot(x, KDE_Q12(:,i)+base,'color','k','linewidth',2)
	text(xmax,base+max(KDE_Q12(:,i))/2,'Q1-Q2','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_Q12(:,i));	
	
	patch([x;xmax;0], [KDE_Q23(:,i)+base;base;base],colors(5,:))
	plot(x, KDE_Q23(:,i)+base,'color','k','linewidth',2)
	text(xmax,base+max(KDE_Q23(:,i))/2,'Q2-Q3','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_Q23(:,i));	
	
	patch([x;xmax;0], [KDE_Q34(:,i)+base;base;base],colors(6,:))
	plot(x, KDE_Q34(:,i)+base,'color','k','linewidth',2)
	text(xmax,base+max(KDE_Q34(:,i))/2,'> Q3','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_Q34(:,i));	
	
	patch([x;xmax;0], [KDE_L100(:,i)+base;base;base],colors(7,:))
	plot(x, KDE_L100(:,i)+base,'color','k','linewidth',2)	
	text(xmax,base+max(KDE_L100(:,i))/2,'Largest 100','fontsize',16, 'horizontalAlignment', 'right')
	base = base + max(KDE_L100(:,i));
	
	title(Name(i))
	
	ylim([0 base])
end






labels = {'WC', 'App', 'PG', 'Gren', 'MC', 'YM', 'TH', 'Arch', 'Other'};





for i = 1:N
	figure

	subplot(1,2,1);
	pie(Pie_S100(i,:),labels)
	subplot(1,2,2);
	pie(Pie_L100(i,:),labels)
		title(Name(i))
end





% 
% for i = 1:N
% 	figure
% 	title(Name(i))
% 	subplot(1,3,1);
% 	pie(Pie_n(i,:),labels)
% 	subplot(1,3,2);
% 	pie(Pie_S100(i,:),labels)
% 	subplot(1,3,3);
% 	pie(Pie_L100(i,:),labels)
% 	title(Name(i))
% end


% 
% figure 
% for i = 1:N
% 	subplot(N,1,i);
% 	pie(Pie_S100(i,:),labels)
% end
% 
% figure 
% for i = 1:N
% 	subplot(N,1,i);
% 	pie(Pie_L100(i,:),labels)
% end
% 
% 
% 
% figure 
% for i = 1:N
% 	subplot(N,1,i);
% 	pie(Pie_S100(i,:))
% end
% 
% figure 
% for i = 1:N
% 	subplot(N,1,i);
% 	pie(Pie_L100(i,:))
% end