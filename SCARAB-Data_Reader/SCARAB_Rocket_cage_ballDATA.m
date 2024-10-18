clc
clear all
close all
%%
load('meteorRocketMass.mat')
load('meteorcageMass.mat')
load('meteorBallMass.mat')

newcolors = [0 1 0
             0 1 0
             0 0 1
             0 0 1
             0 1 1
             0 1 1];
         
newcolors2= [0 1 0
             0 0 1
             0 1 1];
         
files(1)=rocket;
files(2)=ball;
files(3)=cage;

r_planet = 6371000;
                                                                               
lat_OBS=6;
lon_OBS=1;

z_init= 10000;  %OBS_height(ii);%up 
x_init= (z_init+r_planet)*deg2rad(0-lon_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
y_init= (z_init+r_planet)*deg2rad(0-lat_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north

initial_altitude=130;
%% Rocket
ii=1;

figure(1)
colororder(newcolors2)
plot(files(ii).gr_trk(files(ii).alt<initial_altitude)-mean(files(ii).gr_trk(round(files(ii).alt)==initial_altitude)),files(ii).alt(files(ii).alt<initial_altitude),'LineWidth',1.2)
hold on

%% VEL PLOT
figure(2)
colororder(newcolors2)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).vel(files(ii).alt<initial_altitude)/1000,'LineWidth',1.2)
hold on

%% MASS PLOT
figure(3)
subplot(2,1,1)
colororder(newcolors2)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mass(files(ii).alt<initial_altitude)/files(ii).mass(1),'LineWidth',1.2)
hold on

subplot(2,1,2)
colororder(newcolors2)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mass_rate(files(ii).alt<initial_altitude),'LineWidth',1.2)
hold on

%% MAGNITUDE PLOT
figure(4)
colororder(newcolors)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mag((files(ii).alt<initial_altitude),1),'LineWidth',1.2)
hold on
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mag((files(ii).alt<initial_altitude),2),'--','LineWidth',1.2)

figure(5)
colororder(newcolors)
plot(files(ii).mag((files(ii).alt<initial_altitude),1),files(ii).alt(files(ii).alt<initial_altitude),'LineWidth',1.2)
hold on
plot(files(ii).mag((files(ii).alt<initial_altitude),2),files(ii).alt(files(ii).alt<initial_altitude),'--','LineWidth',1.2)

%% Inclination Ang

figure(6)
colororder(newcolors2)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),90+files(ii).flight_ang(files(ii).alt<initial_altitude),'LineWidth',1.2);
hold on

%% EARTH PLOT
figure(7)
colororder(newcolors2)
plot3(files(ii).lon(1:length(files(ii).alt(files(ii).alt<initial_altitude))),files(ii).lat(1:length(files(ii).alt(files(ii).alt<initial_altitude))),files(ii).alt(files(ii).alt<initial_altitude),'LineWidth',1.2); 
hold on
% view([0 0 90])

%% Azimuth Elev

Up_pos = (files(ii).alt(files(ii).alt<initial_altitude)*1000 - z_init);
Nort_pos=  ( 1000*(files(ii).gr_trk(files(ii).alt<initial_altitude)-mean(files(ii).gr_trk(round(files(ii).alt)==initial_altitude))) .* sin( deg2rad(files(ii).head_ang(1:length(files(ii).alt(files(ii).alt<initial_altitude))) )) + y_init );%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
West_pos=  ( 1000*(files(ii).gr_trk(files(ii).alt<initial_altitude)-mean(files(ii).gr_trk(round(files(ii).alt)==initial_altitude))) .* cos( deg2rad(-files(ii).head_ang(1:length(files(ii).alt(files(ii).alt<initial_altitude))) )) + x_init );%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
% Nort_pos=  (files(ii).alt(files(ii).alt<initial_altitude) +r_planet).*deg2rad(files(ii).lon(1:length(files(ii).alt(files(ii).alt<initial_altitude)))-lon_OBS); %( (files(ii).gr_trk(files(ii).alt<initial_altitude)-mean(files(ii).gr_trk(round(files(ii).alt)==initial_altitude))) .* cos( deg2rad(-files(ii).head_ang(files(ii).alt<initial_altitude) )) + x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
% West_pos=  (files(ii).alt(files(ii).alt<initial_altitude) +r_planet).*deg2rad(files(ii).lat(1:length(files(ii).alt(files(ii).alt<initial_altitude)))-lat_OBS); %( (files(ii).gr_trk(files(ii).alt<initial_altitude)-mean(files(ii).gr_trk(round(files(ii).alt)==initial_altitude))) .* sin( deg2rad( files(ii).head_ang(files(ii).alt<initial_altitude) )) + y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 

r4=[Nort_pos West_pos Up_pos];
files(ii).azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
files(ii).elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

figure(8)
skyplot(files(ii).azimuth,files(ii).elevation,'g.')
hold on
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')

%% Ball Meteor

for ii=[2 3]
%% ALTITUDE DOWNRANGE
figure(1)
plot(files(ii).gr_trk,files(ii).alt,'LineWidth',1.2)
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
grid on
hold on
ylim([0 initial_altitude])
ylabel('Altitude [km]')
xlabel('Dowrange [km]')
legend({'Rocket','Meteor','Cage'},'Location','southwest');

%% VEL PLOT
figure(2)
plot(files(ii).time,files(ii).vel/1000,'LineWidth',1.2)
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
grid on
hold on
ylabel('Velocity [km/s]')
xlabel('Time mm:ss')
legend({'Rocket','Meteor','Cage'},'Location','southwest');

%% MASS PLOT
figure(3)
subplot(2,1,1)
plot(files(ii).time,files(ii).mass/files(ii).mass(1),'LineWidth',1.2)
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Mass [-]')

subplot(2,1,2)
plot(files(ii).time(1:end-1),files(ii).mass_rate(1:end-1),'LineWidth',1.2)
grid on
hold on
xlabel('Time mm:ss')
ylabel('Mass loss rate [Kg/s]')
legend({'Rocket','Meteor','Cage'},'Location','southwest');

%% MAGNITUDE PLOT
figure(4)
colororder(newcolors)
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,1),'LineWidth',1.2)
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Abs.Magnitude')
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,2),'--','LineWidth',1.2)
legend({'Rocket \tau = 0.1','Rocket \tau = 0.001','Meteor \tau = 0.1','Meteor \tau = 0.001','Cage \tau = 0.1','Cage \tau = 0.001'},'Location','southwest');

figure(5)
colororder(newcolors)
plot(files(ii).mag(1:end-1,1),files(ii).alt(1:end-1),'LineWidth',1.2)
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
set(gca, 'XDir','reverse')
grid on
hold on
plot(files(ii).mag(1:end-1,2),files(ii).alt(1:end-1),'--','LineWidth',1.2)
ylim([0 100])
xlim([-10 10])
xlabel('Abs.Magnitude') 
ylabel('Altitude [km]')
legend({'Rocket \tau = 0.1','Rocket \tau = 0.001','Meteor \tau = 0.1','Meteor \tau = 0.001','Cage \tau = 0.1','Cage \tau = 0.001'},'Location','southwest');

%% Inclination Ang

figure(6)
plot(files(ii).time,90+files(ii).flight_ang,'LineWidth',1.2);
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
grid on
hold on
xlabel('Time mm:ss') 
ylabel('Zenith Angle [°]')
legend({'Rocket','Meteor','Cage'},'Location','southwest');

%% EARTH PLOT
figure(7)
plot3(files(ii).lon,files(ii).lat,files(ii).alt,'LineWidth',1.2); 
title('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
hold on
grid on
xlabel('Longitude [degrees]','Fontsize',10);
ylabel('Latitude [degrees]','Fontsize',10);
zlabel('Altitude [km]','Fontsize',10);
% view([0 0 90])

%% Azimuth Elevation

Up_pos = (files(ii).alt*1000 - z_init);
Nort_pos=  ( 1000*(files(ii).gr_trk) .* sin( deg2rad(files(ii).head_ang )) + y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
West_pos=  ( 1000*(files(ii).gr_trk) .* cos( deg2rad(-files(ii).head_ang )) + x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 

r4=[Nort_pos West_pos Up_pos];
files(ii).azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
files(ii).elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

end


figure(7)
plot3(lon_OBS,lat_OBS,z_init/1000,'*r','LineWidth',1)   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
lat = load('coast_lat.dat');
long = load('coast_long.dat');
plot3(long,lat,zeros(1,length(lat)),'k')
xlim([-30 30])
ylim([-30 30])
% zlim([0 500])
set(gca,'Xtick',-180:30:180)
set(gca,'Ytick',-90:30:90)
% view([0 0 90])
hleg = legend({'Rocket','Meteor','Cage'},'Location','southwest');

allChildren = get(hleg, 'Children');                % list of all objects on axes
displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
% Remove object associated with "data1" in legend
delete(allChildren(strcmp(displayNames, 'data1')))

%% Azimut Elev

figure(8)
ii=2;
skyplot(files(ii).azimuth,files(ii).elevation,'b.')
ii=3;
skyplot(files(ii).azimuth,files(ii).elevation,'c.')

%% SAVE .jpg

FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\Oct';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    figure(iFig)
    saveas(gcf,fullfile(FolderName, [num2str(iFig),'_SCARAB_IMC.jpg']))
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   savefig(fullfile(FolderName, [FigName '.fig']));
end
ii=1;
files(ii).gr_trk=files(ii).gr_trk(files(ii).alt<initial_altitude)-mean(files(ii).gr_trk(round(files(ii).alt)==initial_altitude));
files(ii).time=files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude));
files(ii).sec=files(ii).sec(files(ii).alt<initial_altitude)-mean(files(ii).sec(round(files(ii).alt)==initial_altitude));
files(ii).mag=[files(ii).mag((files(ii).alt<initial_altitude),1) files(ii).mag((files(ii).alt<initial_altitude),2)];

save('RocketMeteorCage.mat','files')

time_step=200;
figure
ii=1;
skyplot(files(ii).azimuth((round(files(ii).sec(:,1))<time_step)),files(ii).elevation((round(files(ii).sec(:,1))<time_step)),'g.')
hold on
ii=2;
skyplot(files(ii).azimuth((round(files(ii).sec(:,1))<time_step)),files(ii).elevation((round(files(ii).sec(:,1))<time_step)),'b.')
ii=3;
skyplot(files(ii).azimuth((round(files(ii).sec(:,1))<time_step)),files(ii).elevation((round(files(ii).sec(:,1))<time_step)),'c.')