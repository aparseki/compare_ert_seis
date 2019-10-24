clear; close all;clc;

%% terrain
topo =  load('cp_elev.txt'); % need to have a file defined for topography
doi  =  load('cp_ert_doi.txt'); % need to have an ERT DOI defined. Assumed SSR DOI is always greater than ERT
xaxis = 0:0.35:95; %the grid. Need to manually define based on site and desired resolution
yaxis = -10:0.1:0; %same as above
%% load seismic
data = load('cp_seis_result.xzv'); % col1 = x, col2 = z, col3 = velocity, m/s
data(:,1) = data(:,1)+0; % Geophone 1 = elecrode #1. If the ERT and Sesismic are not co-registered so GP1=Elx1, this can be fixed here
data(:,3) = data(:,3)/1000; % convert to km/s
%% interpolate seismic
[Xgrid, Ygrid, vpField]=griddata(data(:,1),data(:,2),data(:,3),xaxis,yaxis'); %interpolate data onto rectilinear grid
vpField(isnan(vpField)) = -999; % just changes NANs to numbers
vpField = fliplr(vpField); %if the lines were measured in opposite direction, use this, otherwise comment out.
%% ert (measured at the same time or similar as seismic)
%data = load('2016_05_23/f001_res.dat'); % col1 = x, col2 = z, col3 = resistivity, ohm m, col4 = resistivity log10(ohm m)
data = load('f001_res.dat'); %  col1 = x, col2 = z, col3 = resistivity, ohm m, col4 = resistivity log10(ohm m)
%data(:,2) = data(:,2)+2565.4757;

[Xgrid, Ygrid, logresField1]=griddata(data(:,1),data(:,2),data(:,4),xaxis,yaxis'); %regridding the ERT data
logresField1(isnan(logresField1)) = -999; %just getting rid of NANS

%% removing data not in the zone of interest
cnt = 1;
for i = 1:length(xaxis) % loop though and check to be sure the values are below toporaphy
    for j = 1:length(yaxis)
        a1 = topo(:,1)-xaxis(i);
        a2 = find(a1==min(abs(a1)));
        if yaxis(j)<topo(a2,2)
            datOut(cnt,:) = [xaxis(i) yaxis(j) logresField1(j,i) vpField(j,i) yaxis(j)-topo(a2,2)]; %tabulates X, Y, rho, VP, and distance below the surface
            cnt = cnt+1;
        end
    end
end

cnt = 1;
for i = 1:length(datOut); %loops through to make sure points only within DOI are kept
    b1 = doi(:,1)-datOut(i,1);
    b2 = find(b1==min(abs(b1)));
    if datOut(i,2)>doi(b2,2)
        datOut2(cnt,:) = datOut(i,:);
        cnt = cnt+1;
    end
end

%% Makes the first figure just showing both datasets
close all

figure
imagesc(xaxis,yaxis,logresField1); caxis([1 3]); set(gca,'ydir','normal'); hold on
plot(topo(:,1),topo(:,2),'-k'); plot(doi(:,1),doi(:,2),'--k'); hold on
[c,h] =contour(xaxis,yaxis,vpField,0.5:0.05:0.8,'color','w');
clabel(c,h,'color','w')
patch([topo(1,1)-1; topo(:,1); topo(1,1)],[topo(1,2); topo(:,2); topo(end,2)],'w')
ylim([min(doi(:,2))-2 max(topo(:,2))]) %automatically set the y-axis limits based on topo and ERT doi
xlabel('distance, m')
ylabel('elevation, m')
title('resistivity, ohm m; Vp, km s^-^1')
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',11 ) 
set(findall(gcf,'-property','FontName'),'FontName','Avenir' ) 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 8]) 
outname = 'rho_vp.jpg';
print(outname,'-djpeg','-r600')
close all
%% Makes the second figure showing the first ratio result

% figure
% imagesc(xaxis,yaxis,vpField); caxis([.5 1.2]);set(gca,'ydir','normal'); hold on
% plot(topo(:,1),topo(:,2),'-k')

figure
imagesc(xaxis,yaxis,vpField./logresField1);caxis([.1 .7]); set(gca,'ydir','normal'); hold on
plot(topo(:,1),topo(:,2),'-k'); plot(doi(:,1),doi(:,2),'--k')
patch([topo(1,1)-1; topo(:,1); topo(1,1)],[topo(1,2); topo(:,2); topo(end,2)],'w')
ylim([min(doi(:,2))-2 max(topo(:,2))]) %automatically set the y-axis limits based on topo and ERT doi
xlabel('distance, m')
ylabel('elevation, m')
title('ratio Vp/log_1_0(\rho), km s^-^1/log_1_0(ohm m)')
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',11 ) 
set(findall(gcf,'-property','FontName'),'FontName','Avenir' ) 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 8]) 
outname = 'ratio1.jpg';
print(outname,'-djpeg','-r600')
close all
%% Makes the third figure showing the alternative ratio result
figure
imagesc(xaxis,yaxis,real(log10(1000.*vpField))./logresField1);caxis([.5 2]); set(gca,'ydir','normal'); hold on
plot(topo(:,1),topo(:,2),'-k'); plot(doi(:,1),doi(:,2),'--k')
patch([topo(1,1)-1; topo(:,1); topo(1,1)],[topo(1,2); topo(:,2); topo(end,2)],'w')
ylim([min(doi(:,2))-2 max(topo(:,2))]) %automatically set the y-axis limits based on topo and ERT doi
xlabel('distance, m')
ylabel('elevation, m')
title('ratio log_1_0(Vp)/log_1_0(\rho), log_1_0(m s^-^1)/log_1_0(ohm m)')
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',11 ) 
set(findall(gcf,'-property','FontName'),'FontName','Avenir' ) 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 8]) 
outname = 'ratio2.jpg';
print(outname,'-djpeg','-r600')
close all

%% Makes the Rock Physics Scatterplot, color coded by depth, following Meju 2003
%
% Meju, M. A., L. A. Gallardo, and A. K. Mohamed, Evidence for correlation of 
% electrical resistivity and seismic velocity in heterogeneous near-surface
% materials, Geophys. Res. Lett., 30(7), 1373, doi:10.1029/2002GL016048, 2003.

figure
subplot(1,2,1)
scatter(datOut2(:,3),datOut2(:,4),5,datOut2(:,5),'filled','o')
xlim([1 3.5])
ylim([0 .9])
xlabel('log10(\rho)')
ylabel('Vp')
caxis([min(datOut2(:,5)) max(datOut2(:,5))])

subplot(1,2,2)
caxis([min(datOut2(:,5)) max(datOut2(:,5))])
colorbar west
axis off

set(findall(gcf,'-property','FontSize'),'FontSize',11 ) 
set(findall(gcf,'-property','FontName'),'FontName','Avenir' ) 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 6]) 
outname = 'rpplot.jpg';
print(outname,'-djpeg','-r600')
close all