clc
clear all
close all
%%
load('RocketMeteorCage.mat')

hold on
col=['g.' 'b.' 'c.'];
for t=1:300
    for ii=1:3
        skyplot(mean(files(ii).azimuth((round(files(ii).sec(:,1))==t-1))),files(ii).elevation((round(files(ii).sec(:,1))==t-1)),col(ii));
    end
     pause(.1)
end

% n = 100;
% colors = jet(100);
% x = linspace(0,10,n);
% y = sin(x);
% cla
% xlim([min(x) max(x)])
% ylim([min(y) max(y)])
% hold on
% for i = 1:n-1
%     xx=90-mean(files(ii).elevation((round(files(ii).sec(:,1))==t-1))).*sin(mean(files(ii).azimuth((round(files(ii).sec(:,1))==t-1)))/180*pi;
%     yy=90-mean(files(ii).elevation((round(files(ii).sec(:,1))==t-1))).*cos(mean(files(ii).azimuth((round(files(ii).sec(:,1))==t-1)))/180*pi;
%     plot(x(i:i+1),y(i:i+1),'color',colors(i,:))
%     pause(0.1)
% end
% hold off


% % Make the first frame: 
% f = figure;
% skyplot(0,0,'k');
% hold on
% for ii=1:2
% % transform data to Cartesian coordinates.
% yy = (90-files(ii).elevation(1)).*cos(files(ii).azimuth(1)/180*pi);
% xx = (90-files(ii).elevation(1)).*sin(files(ii).azimuth(1)/180*pi);
% 
% files(ii).mag(isnana(files(ii).mag),2) = 20;
% 
% h(ii) = scatter(xx,yy,60,files(ii).mag(1,2),'filled'); 
% end 
% 
% % set color axis limits: 
% caxis([-10 10]) 
% 
% cb = colorbar; 
% ylabel(cb,'Magnitude') 
% 
% colormap(f, flipud(colormap(f)))
% f.WindowState = 'maximized';
% 
% % write the first frame: 
% gif('temperaturedata.gif') % <-uncomment to make a gif
% 
% % Loop through each subsequent time step: 
% for  t = 2:300
%     for ii=1:3
%         if isnan(files(ii).mag(round(files(ii).sec(:,1))==t-1))
%          set(h(ii),'xdata',(90-mean(files(ii).elevation((round(files(ii).sec(:,1))==t-1))).*sin(mean(files(ii).azimuth((round(files(ii).sec(:,1))==t-1)))/180*pi),...
%             'ydata',(90-mean(files(ii).elevation((round(files(ii).sec(:,1))==t-1))).*cos(mean(files(ii).azimuth((round(files(ii).sec(:,1))==t-1)))/180*pi),...
%             'cdata',mean(files(ii).mag((round(files(ii).sec(:,1))==t-1))),2))
%         end
%     end
%      %pause(0.1) % <-not necessary for making a gif
%      gif  <-uncomment to make a gif
%  end