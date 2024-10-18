clc
clear all
close all
%% CHANGE FOR SIMULATION

sec_asteroid=300;%130
% sec_asteroid=180;%110 or 130 5
initial_altitude=130;
% asteroidNum=10;

load('meteorRocketMass.mat')
load('meteorcageMass.mat')
load('meteorBallMass.mat')
load('meteorMissionMass.mat')

for asteroidNum=1:length(files_mission) %10%[4 10 13 19] %1:length(files_mission) 
files(1)=rocket;
if asteroidNum==10
    load('meteorBallMass.mat')
    files(2)=ball;
else
    files(2)=files_mission(asteroidNum);
end
files(3)=cage;
fontSize = 15;
switch files(2).material
    case "iron"
        titles=append('Iron Meteor ',num2str(files(2).init_size),' cm');
    case "basaltmax"
        titles=append('Basalt Acidic Meteor ',num2str(files(2).init_size),' cm');
    case "minbasalt"
        titles=append('Basalt Basic Meteor ',num2str(files(2).init_size),' cm');
    case "ordchon"
        titles=append('Ordinary Chondrite Meteor ',num2str(files(2).init_size),' cm');
    case "carchon"
        titles=append('Carbonaceous Chondrite Meteor ',num2str(files(2).init_size),' cm');
    case "granite"
        titles=append('Granite Meteor ',num2str(files(2).init_size),' cm');
    case "sandstone"
        titles=append('Sandstone Meteor ',num2str(files(2).init_size),' cm');
end

%% DATA

lum_eff = [0.1 0.001];

newcolors = [0 1 0
             0 1 0
             0 0 1
             0 0 1
             0 1 1
             0 1 1];
         
newcolors2= [0 1 0
             0 0 1
             0 1 1];

r_planet = 6371000;
                                                                               
lat_OBS=23; %best 23 or closer to the second stage 22
lon_OBS=-5; %best -5

z_init= 10000;  %OBS_height(ii);%up 
x_init= (z_init+r_planet)*deg2rad(0-lon_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
y_init= (z_init+r_planet)*deg2rad(0-lat_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north

%% Rocket
ii=1;

tnew=files(ii).time;
tnew.Format='mm:ss';
files(ii).time=tnew;
load('RocketSARA.mat')

lonRocket=-abs(files(ii).lon(files(ii).alt<initial_altitude)-mean(files(ii).lon(round(files(ii).alt)==initial_altitude)));
lonRocket(lonRocket<-180)=-(360+lonRocket(lonRocket<-180));
latRocket=-(files(ii).lat(files(ii).alt<initial_altitude)-mean(files(ii).lat(round(files(ii).alt)==initial_altitude)));
Up_pos = (files(ii).alt(files(ii).alt<initial_altitude).*1000 - z_init);
Nort_pos=(files(ii).alt(files(ii).alt<initial_altitude).*1000+r_planet).*deg2rad(latRocket-lat_OBS);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
West_pos=(files(ii).alt(files(ii).alt<initial_altitude).*1000+r_planet).*deg2rad(lonRocket-lon_OBS);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
Distance = sqrt((West_pos).^2 + (Nort_pos).^2 + (Up_pos).^2 );  % distance respec to the observer

files(ii).REAL_lum_ene = 0.5 .* files(ii).vel(files(ii).alt<initial_altitude).^2 .* files(ii).mass_rate(files(ii).alt<initial_altitude).* lum_eff .* 1e10 ./ (Distance.^2);
files(ii).REAL_lum_ene(files(ii).REAL_lum_ene==0)=nan;
files(ii).REAL_mag = 6.8 - 1.086 * log(files(ii).REAL_lum_ene);

lonRocket=saraRocket(:,4);
latRocket=saraRocket(:,3);
altRocket=saraRocket(:,2);
Up_pos = (altRocket.*1000 - z_init);
Nort_pos=(altRocket.*1000+r_planet).*deg2rad(latRocket-lat_OBS);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
West_pos=(altRocket.*1000+r_planet).*deg2rad(lonRocket-lon_OBS);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
r4=[Nort_pos West_pos Up_pos];
files(ii).azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
files(ii).elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));


figure(1)
subplot(2,1,1)
colororder(newcolors2)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),smooth(files(ii).vel(files(ii).alt<initial_altitude)/1000,20),'LineWidth',1.2)
hold on

subplot(2,1,2)
colororder(newcolors)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mag((files(ii).alt<initial_altitude),1),'LineWidth',1.2)
hold on
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mag((files(ii).alt<initial_altitude),2),'--','LineWidth',1.2)
ylim([-15 15])

figure(2)
colororder(newcolors)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).REAL_mag(:,1),'LineWidth',1.2)
hold on
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).REAL_mag(:,2),'--','LineWidth',1.2)
ylim([-10 10])

% figure(3)
% colororder(newcolors2) 
% % load('RocketSARA.mat')
% % plot3(saraRocket(:,4),saraRocket(:,3),saraRocket(:,2),'LineWidth',1.2); 
% hold on



% figure(4)
% skyplot(files(ii).azimuth,files(ii).elevation,'g.')
% hold on

%% METEOR

ii=2;

tnew=files(ii).time;
tnew.Format='mm:ss';
files(ii).time=tnew;
load('RocketSARA.mat')

Up_pos = (files(ii).alt*1000 - z_init);
Nort_pos=  ( 1000*(files(ii).gr_trk) .* sin( deg2rad(files(ii).head_ang )) + y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
West_pos=  ( 1000*(files(ii).gr_trk) .* cos( deg2rad(-files(ii).head_ang )) + x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
Distance = sqrt((West_pos).^2 + (Nort_pos).^2 + (Up_pos).^2 );  % distance respec to the observer

files(ii).REAL_lum_ene = 0.5 .* files(ii).vel.^2 .* files(ii).mass_rate.* lum_eff .* 1e10 ./ (Distance.^2);
files(ii).REAL_lum_ene(files(ii).REAL_lum_ene==0)=nan;
files(ii).REAL_mag = 6.8 - 1.086 * log(files(ii).REAL_lum_ene);

r4=[Nort_pos West_pos Up_pos];
files(ii).azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
files(ii).elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

rocketElevDist=deg2rad(abs(90-files(1).elevation(1:length(files(2).time))));
meteorElevDist=deg2rad(abs(90-files(2).elevation));
% rocketElevDist=deg2rad(files(1).elevation(1:length(files(2).time)));
% meteorElevDist=deg2rad(files(2).elevation);
meteorAzAngrocket=deg2rad(files(1).azimuth(1:length(files(2).time)))-deg2rad(files(2).azimuth); % distance
check_az=[(files(1).azimuth(1:length(files(2).time))) (files(2).azimuth)];

Separat_Ang=rad2deg(sqrt(rocketElevDist.^2+meteorElevDist.^2-2*rocketElevDist.*meteorElevDist.*cos(meteorAzAngrocket)));
    
allang=[rad2deg(rocketElevDist) rad2deg(meteorElevDist) rad2deg(meteorAzAngrocket) Separat_Ang];

if asteroidNum==10
    files(ii).REAL_mag(files(ii).sec<sec_asteroid,1)=nan;
    files(ii).REAL_mag(files(ii).sec<sec_asteroid,2)=nan;
    files(ii).mag(files(ii).sec<sec_asteroid,1)=nan;
    files(ii).mag(files(ii).sec<sec_asteroid,2)=nan;
end

figure(3)
sz = 5; % size of the points they can vary if use a vector
c=files(ii).REAL_mag(:,1);
c(isnan(c))=20;
scatter(files(2).time,Separat_Ang,sz,c,'filled')
grid on
hold on
scatter(files(2).time(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),Separat_Ang(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),25,c(files(ii).REAL_mag(:,1)==min(files(ii).REAL_mag(:,1))),'filled')
% plot(files(2).time(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),Separat_Ang(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),'y*')
colormap(winter)
set(gca,'yscale','log')
xlabel('Time mm:ss')
ylabel('Separation Angle [°] ') %separation angle \theta
title(titles, 'FontSize', fontSize)
h = colorbar('Direction','reverse');
ylabel(h, '\tau = 0.1 Magnitude')
fig3 = figure(3);
colormap(fig3, flipud(colormap(fig3)))
clim([-10 10])




% if files(2).init_size==5
% fig3 = figure(3);
% colormap(fig3, flipud(colormap(fig3)))
% end

xlim([files(3).time(1) files(2).time(end)])

figure(1)
subplot(2,1,1)
colororder(newcolors2)
plot(files(ii).time,files(ii).vel/1000,'b','LineWidth',1.2)
grid on
xlim([files(ii).time(1) files(ii).time(end)])

subplot(2,1,2)
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,1),'LineWidth',1.2)
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,2),'--','LineWidth',1.2)
xlim([files(ii).time(1) files(ii).time(end)])
% colororder(newcolors)
% time_consid=files(ii).time(files(ii).sec>sec_asteroid);
% MagMAX=files(ii).mag(files(ii).sec>sec_asteroid,1);
% Magmin=files(ii).mag(files(ii).sec>sec_asteroid,2);
% plot(time_consid(1:end-1),MagMAX(1:end-1),'LineWidth',1.2)
% plot(time_consid(1:end-1),Magmin(1:end-1),'--','LineWidth',1.2)

figure(2)
colororder(newcolors)
plot(files(ii).time(1:end-1),files(ii).REAL_mag(1:end-1,1),'LineWidth',1.2)
plot(files(ii).time(1:end-1),files(ii).REAL_mag(1:end-1,2),'--','LineWidth',1.2)

% time_consid=files(ii).time(files(ii).sec>sec_asteroid);
% RealMagMAX=files(ii).REAL_mag(files(ii).sec>sec_asteroid,1);
% RealMagmin=files(ii).REAL_mag(files(ii).sec>sec_asteroid,2);
% plot(time_consid(1:end-1),RealMagMAX(1:end-1),'LineWidth',1.2)
% plot(time_consid(1:end-1),RealMagmin(1:end-1),'--','LineWidth',1.2)




% figure(3)
% plot3(files(ii).lon,files(ii).lat,files(ii).alt,'LineWidth',1.2); 

% figure(4)
% skyplot(files(ii).azimuth,files(ii).elevation,'b.')

%% CAGE

ii=3;

tnew=files(ii).time;
tnew.Format='mm:ss';
files(ii).time=tnew;
load('RocketSARA.mat')

Up_pos = (files(ii).alt*1000 - z_init);
Nort_pos=  ( 1000*(files(ii).gr_trk) .* sin( deg2rad(files(ii).head_ang )) + y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
West_pos=  ( 1000*(files(ii).gr_trk) .* cos( deg2rad(-files(ii).head_ang )) + x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
Distance = sqrt((West_pos).^2 + (Nort_pos).^2 + (Up_pos).^2 );  % distance respec to the observer

files(ii).REAL_lum_ene = 0.5 .* files(ii).vel.^2 .* files(ii).mass_rate.* lum_eff .* 1e10 ./ (Distance.^2);
files(ii).REAL_lum_ene(files(ii).REAL_lum_ene==0)=nan;
files(ii).REAL_mag = 6.8 - 1.086 * log(files(ii).REAL_lum_ene);

r4=[Nort_pos West_pos Up_pos];
files(ii).azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
files(ii).elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));
    
figure(1)
subplot(2,1,1)
title(titles, 'FontSize', fontSize)
plot(files(ii).time,files(ii).vel/1000,'c','LineWidth',1.2)
xlabel('Time mm:ss')
ylabel('Velocity [km/s]')
hleg = legend({'Rocket','Meteor','Cage'},'Location','southwest');


subplot(2,1,2)
colororder(newcolors)
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,1),'LineWidth',1.2)
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Abs.Magnitude')
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,2),'--','LineWidth',1.2)
legend({'Rocket \tau = 0.1','Rocket \tau = 0.001','Meteor \tau = 0.1','Meteor \tau = 0.001','Cage \tau = 0.1','Cage \tau = 0.001'},'Location','southwest');

figure(2)
title(titles, 'FontSize', fontSize)
colororder(newcolors)
plot(files(ii).time(1:end-1),files(ii).REAL_mag(1:end-1,1),'LineWidth',1.2)
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Magnitude')
plot(files(ii).time(1:end-1),files(ii).REAL_mag(1:end-1,2),'--','LineWidth',1.2)
legend({'Rocket \tau = 0.1','Rocket \tau = 0.001','Meteor \tau = 0.1','Meteor \tau = 0.001','Cage \tau = 0.1','Cage \tau = 0.001'},'Location','northwest');
xlim([files(3).time(1) files(2).time(end)])


% figure(3)
% plot3(files(ii).lon,files(ii).lat,files(ii).alt,'LineWidth',1.2); 
% title(titles, 'FontSize', fontSize)
% hold on
% grid on
% xlabel('Longitude [degrees]','Fontsize',10);
% ylabel('Latitude [degrees]','Fontsize',10);
% zlabel('Altitude [km]','Fontsize',10);
% % view([0 0 90])


figure(4)
leggendsSKY(1)=skyplot(0,90,'k.');
% skyplot(files(ii).azimuth(files(ii).elevation>0),files(ii).elevation(files(ii).elevation>0),'k.')
hold on
ii=1;
% skyplot(files(ii).azimuth(files(ii).elevation>0),files(ii).elevation(files(ii).elevation>0),'k.')
yy = (90-files(ii).elevation(files(ii).elevation>0)).*cos(files(ii).azimuth(files(ii).elevation>0)/180*pi);
xx = (90-files(ii).elevation(files(ii).elevation>0)).*sin(files(ii).azimuth(files(ii).elevation>0)/180*pi);
leggendsSKY(2)=plot(xx,yy,'k.','DisplayName','Dark Flight');
ii=2;
% skyplot(files(ii).azimuth,files(ii).elevation,'k.')
yy = (90-files(ii).elevation).*cos(files(ii).azimuth/180*pi);
xx = (90-files(ii).elevation).*sin(files(ii).azimuth/180*pi);
leggendsSKY(3)=scatter(xx,yy,sz,c,'filled','DisplayName','Meteor');
leggendsSKY(4)=scatter(xx(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),yy(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),25,c(files(ii).REAL_mag(:,1)==min(files(ii).REAL_mag(:,1))),'filled','DisplayName','MAX Meteor');
colormap(winter)
% skyplot(files(ii).azimuth(~isnan(files(ii).REAL_mag(:,2))),files(ii).elevation(~isnan(files(ii).REAL_mag(:,2))),'r.')
% skyplot(files(ii).azimuth(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),files(ii).elevation(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),'y*')
ii=3;
% skyplot(files(ii).azimuth(~isnan(files(ii).REAL_mag(:,2))),files(ii).elevation(~isnan(files(ii).REAL_mag(:,2))),'c.')
yy = (90-files(ii).elevation).*cos(files(ii).azimuth/180*pi);
xx = (90-files(ii).elevation).*sin(files(ii).azimuth/180*pi);
leggendsSKY(5)=plot(xx,yy,'c.','DisplayName','Cage');
ii=1;
% skyplot(files(ii).azimuth(~isnan(files(ii).REAL_mag(:,2))),files(ii).elevation(~isnan(files(ii).REAL_mag(:,2))),'g.')
yy = (90-files(ii).elevation(~isnan(files(ii).REAL_mag(:,1)))).*cos(files(ii).azimuth(~isnan(files(ii).REAL_mag(:,1)))/180*pi);
xx = (90-files(ii).elevation(~isnan(files(ii).REAL_mag(:,1)))).*sin(files(ii).azimuth(~isnan(files(ii).REAL_mag(:,1)))/180*pi);
leggendsSKY(6)=plot(xx,yy,'g.','DisplayName','Rocket');

clim([-10 10])
fig4 = figure(4);
colormap(fig4, flipud(colormap(fig4)))
h = colorbar('Direction','reverse');

ylabel(h, '\tau = 0.1 Magnitude')
legend(leggendsSKY(2:6),'location','northwest');
% if files(2).init_size==5
% fig4 = figure(4);
% colormap(fig4, flipud(colormap(fig4)))
% end


title(titles, 'FontSize', fontSize)

%% Plot EARTH
% figure(3)
% 
% ii=2; 
% plot3(files(ii).lon(isnan(files(ii).REAL_mag(:,2))),files(ii).lat(isnan(files(ii).REAL_mag(:,2))),files(ii).alt(isnan(files(ii).REAL_mag(:,2))),'k.')
% hold on
% ii=1;
% plot3(saraRocket(isnan(files(ii).REAL_mag(:,2)),4),saraRocket(isnan(files(ii).REAL_mag(:,2)),3),saraRocket(isnan(files(ii).REAL_mag(:,2)),2),'k.')
% plot3(saraRocket(end-100:end,4),saraRocket(end-100:end,3),saraRocket(end-100:end,2),'k.')
% 
% plot3(lon_OBS,lat_OBS,z_init/1000,'*b','LineWidth',1)   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
% lat = load('coast_lat.dat');
% long = load('coast_long.dat');
% plot3(long,lat,zeros(1,length(lat)),'k')
% xlim([-30 30])
% ylim([-30 30])
% % zlim([0 500])
% set(gca,'Xtick',-180:30:180)
% set(gca,'Ytick',-90:30:90)
% % view([0 0 90])
% hleg = legend({'Rocket','Meteor','Cage'},'Location','southwest');
% 
% allChildren = get(hleg, 'Children');                % list of all objects on axes
% displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
% % Remove object associated with "data1" in legend
% delete(allChildren(strcmp(displayNames, 'data1')))

%% SAVE .jpg

FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\material_mag';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    f = figure(iFig);
        set(gca,'FontSize',20)
%     f.WindowState = 'maximized';
    saveas(gcf,fullfile(FolderName, [files(2).material,' Meteor ',num2str(files(2).init_size),' cm ',num2str(iFig),'_SCARAB_mag.jpg']))
end

% colnames = {'sec','DefAng','Mag 0.1','Mag 0.001'};
sec=files(2).sec;
DefAng=Separat_Ang;
Mag_MAX=files(2).REAL_mag(:,1);
Mag_min=files(2).REAL_mag(:,2);
txtFile = table(sec,DefAng,Mag_MAX,Mag_min);
writetable(txtFile,fullfile(FolderName, append(files(2).material,' Meteor ',num2str(files(2).init_size),' cm data_SCARAB_mag.txt')),'Delimiter',' ');
% txtFile=[files(2).sec Separat_Ang files(2).REAL_mag(:,1) files(2).REAL_mag(:,2)];
% save(txtFile,fullfile(FolderName, append(titles,'.dat')));

clear files fig3 fig4 c
close all
end