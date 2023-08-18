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
ymin = 500;
ymax = 10000;
extra = 100;

sub = 100;

% 1D kernel bandwidth
kernel = 15;

x = [xmin:xint:xmax]';

% 2D kernel bandwidths
bandwidth_x = 25;
bandwidth_y = 75;

% how many pixels for the images, has to be in powers of 2, no need to go over go over 2^12, results look the same
gridspc = 2^9;

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
		
end


%%%%%% Bivariate KDEs

% Filter out any data not in the range set above
for k = 1:N
	for i = 1:length(data_tmp(:,1))
		if cellfun('isempty', data_tmp(i,k*2-1)) == 0 && cellfun('isempty', data_tmp(i,k*2)) == 0
			if cell2num(data_tmp(i,k*2-1)) >= xmin && cell2num(data_tmp(i,k*2-1)) <= xmax && ...
					cell2num(data_tmp(i,k*2)) >= ymin && cell2num(data_tmp(i,k*2)) <= ymax
				data1(i,k*2-1:k*2) = cell2num(data_tmp(i,k*2-1:k*2));
			end
		end
	end
end

% set min/max ranges for kde2d function
MIN_XY=[xmin-extra,ymin];
MAX_XY=[xmax+extra,ymax];

cmap = cmocean('balance',300);
cmap(1:3,:) = 1; %clip zero vals at 99% and set as white space

% Make and plot bivariate kdes for samples, save as 3D matrix block
for k = 1:N
	data2 = data1(:,k*2-1:k*2);
	data2 = data2(any(data2 ~= 0,2),:);
	[bandwidth1,density1(:,:,k),X1,Y1] = kde2d_set_kernel(data2, gridspc, MIN_XY, MAX_XY, bandwidth_x, bandwidth_y);
	density1(:,:,k) = density1(:,:,k)./sum(sum(density1(:,:,k)));
	
	figure
	hold on
	
	%plot density
	surf(X1,Y1,density1(:,:,k));
	
	%make contours
	max_density1 = max(max(density1(:,:,k)));
	perc = 0.99; % percent of bivariate KDE from peak
	max_density_conf(k,1) = max_density1*(1-perc); % contour from peak density 
	
	
	% format plots
	colormap(cmap)
	shading interp
	view(2)
	title(Name(k,1),'FontSize',40)
	xlabel('Age (Ma)','FontSize',20)
	ylabel('Grainsize (µm^2)','FontSize',20)
	axis([xmin xmax ymin ymax])
	set(gca,'FontSize',20)
	
	%plot contours
	F1 = figure;
	hold on
	contour3(X1,Y1,density1(:,:,k),[max_density_conf(k,1) max_density_conf(k,1)],'k', 'LineWidth', 3);
	grid off
	view(2)
	axis([xmin xmax ymin ymax])
	%[file,path] = uiputfile('*.eps','Save file'); print(F1,'-depsc','-painters',[path file]); epsclean([path file]); % save simplified contours
	

	
% 	plot3([xmin-extra,xmax+extra],[max(cell2num(All(k).data_S100(:,2))),max(cell2num(All(k).data_S100(:,2)))],[1000,1000],'k','linewidth',2)
% 	plot3([xmin-extra,xmax+extra],[max(cell2num(All(k).dataQ01(:,2))),max(cell2num(All(k).dataQ01(:,2)))],[1000,1000],'k','linewidth',2)
% 	plot3([xmin-extra,xmax+extra],[max(cell2num(All(k).dataQ12(:,2))),max(cell2num(All(k).dataQ12(:,2)))],[1000,1000],'k','linewidth',2)
% 	plot3([xmin-extra,xmax+extra],[max(cell2num(All(k).dataQ23(:,2))),max(cell2num(All(k).dataQ23(:,2)))],[1000,1000],'k','linewidth',2)
% 	plot3([xmin-extra,xmax+extra],[min(cell2num(All(k).data_L100(:,2))),min(cell2num(All(k).data_L100(:,2)))],[1000,1000],'k','linewidth',2)
% 	axis([xmin xmax ymin ymax])
% 	view(2)
	
	%text(xmax+extra,max(cell2num(All(k).data_S100(:,2))),1000,'Smallest 100','fontsize',16, 'horizontalAlignment', 'right')
	
end




