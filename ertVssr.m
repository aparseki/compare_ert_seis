function ertVssr(Topo, Doi,Data, zero_offset, Er_data, xaxis, yaxis, VpContRng, VpThresh,rpRange,corr_flg,site)

%% ERT vs. SSR: compares colocated ERT and Vp images
% topo         = elevation [x | z] meters
% doi          = ERT depth of investigation [x | z] meters depth
% data         = ssr inverted result, xzv format from Geogigga, = m/s
% zero_offfset = difference between start of ssr and ert lines.
% er_data      = name of inverted ERT result file, output from R2
% xaxis        = set the Yaxis for the mesh with interval, eg min#:intvl#:max#
% yaxis        = set the Yaxis for the mesh with interval, eg min#:intvl#:max#
% VpThresh     = threshold minimum Vp values to ignore, suggest ~0.3km/s
% rpRange      = plotting ranges [minRho maxRho minVp maxVp]
% corr_flg     = flag indicating (1 = yes) if topo correction needed [DOI ERT_mesh]
% site         = name of site for figure title
%% terrain
topo =  load(Topo);
%topo = [topo(:,1)+31 topo(:,2)];
doi  =  load(Doi);
if corr_flg(1) == 1;
    doi(:,2)  =  doi(:,2)+min(topo(:,2)); %correct DOI for elevation
end
%% load seismic
data = load(Data); % col1 = x, col2 = z, col3 = velocity, m/s
data(:,1) = data(:,1)+zero_offset; % Geophone 1 = elecrode #1.
vVals = data(:,3)/1000; % convert to km/s
vVals(vVals<0) = max(vVals);
data(:,3) = vVals;
%% interpolate seismic
[~, ~, vpField]=griddata(data(:,1),data(:,2),data(:,3),xaxis,yaxis','natural');
%vpField(isnan(vpField)) = max(max(vpField));
vpField(isnan(vpField)) = 0;
%% ert (measured at the same time or similar as seismic)
er_data = load(Er_data); %  col1 = x, col2 = z, col3 = resistivity, ohm m, col4 = resistivity log10(ohm m)
if corr_flg(2) == 1;
    er_data(:,2) = er_data(:,2)+min(topo(:,2)); %if needed, correct ERT mesh to topography
end
[~, ~, logresField1]=griddata(er_data(:,1),er_data(:,2),er_data(:,4),xaxis,yaxis');
logresField1(isnan(logresField1)) = -999;
%% removing data not in the zone of interest
cnt = 1;
% first loop removes anything above the topography
for i = 1:length(xaxis)
    for j = 1:length(yaxis)
        a1 = topo(:,1)-xaxis(i);
        a2 = find(a1==min(abs(a1)));
        if yaxis(j)<topo(a2,2)
            datOut(cnt,:) = [xaxis(i) yaxis(j) logresField1(j,i) vpField(j,i) yaxis(j)-topo(a2,2)];
            cnt = cnt+1;
        end
    end
end

%second loop removes anything below ERT DOI
cnt = 1;
for i = 1:length(datOut);
    b1 = doi(:,1)-datOut(i,1);
    b2 = find(b1==min(abs(b1)));
    if datOut(i,2)>doi(b2,2)
        datOutA(cnt,:) = datOut(i,:);
        cnt = cnt+1;
    end
end

%third loop removes anything outside SSR DOI
cnt = 1;
for i = 1:length(datOutA);
    if datOutA(i,4)>VpThresh
        datOut2(cnt,:) = datOutA(i,:);
        cnt = cnt+1;
    end
end

%% rho and Vp data co-plotted
close all

figure
imagesc(xaxis,yaxis,logresField1); caxis([1 3]); set(gca,'ydir','normal'); hold on
plot(topo(:,1),topo(:,2),'-k'); plot(doi(:,1),doi(:,2),'--k'); hold on
[c,h] =contour(xaxis,yaxis,vpField,VpContRng,'color','w');
clabel(c,h,'color','w')
patch([topo(1,1); topo(:,1); topo(1,1)],[topo(1,2); topo(:,2); topo(end,2)],'w')
xlabel('distance, m')
ylabel('elevation, m')
title('resistivity, ohm m; Vp, km s^-^1')
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',11 )
set(findall(gcf,'-property','FontName'),'FontName','Avenir' )
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 8])
outname = 'rho_vp.jpg';
print([site outname],'-djpeg','-r600')
close all

%% Vp/er ratio #2
figure
imagesc(xaxis,yaxis,real(log10(1000.*vpField))./logresField1);caxis([.5 2]); set(gca,'ydir','normal'); hold on
plot(topo(:,1),topo(:,2),'-k'); plot(doi(:,1),doi(:,2),'--k')
patch([topo(1,1); topo(:,1); topo(1,1)],[topo(1,2); topo(:,2); topo(end,2)],'w')
xlabel('distance, m')
ylabel('elevation, m')
title('ratio log_1_0(Vp)/log_1_0(\rho), log_1_0(m s^-^1)/log_1_0(ohm m)')
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',11 )
set(findall(gcf,'-property','FontName'),'FontName','Avenir' )
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 8])
outname = 'ratio2.jpg';
print([site outname],'-djpeg','-r600')
close all

%% RockPhysics Plot (rrplot)
figure
scatter(datOut2(:,4),datOut2(:,3),20,datOut2(:,5),'filled','o')
ylim(rpRange(1:2))
xlim(rpRange(3:4))
ylabel('log10(\rho) [Ohm m]')
xlabel('Vp [km s^-^1]')
colorbar('location','north'); caxis([-5 0])
title('depth below surface [m]')

set(findall(gcf,'-property','FontSize'),'FontSize',11 )
set(findall(gcf,'-property','FontName'),'FontName','Avenir' )
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8 8])
outname = 'rpplot.jpg';
print([site outname],'-djpeg','-r600')
close all