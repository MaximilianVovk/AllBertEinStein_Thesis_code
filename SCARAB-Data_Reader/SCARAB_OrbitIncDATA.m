clc
clear all
close all
%%
load('meteorIronMass.mat')%[4 5 7]

newcolors = [0 0.4470 0.7410
             0 0.4470 0.7410
             0.8500 0.3250 0.098
             0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330
             0.3010 0.7450 0.9330
             0.6350 0.0780 0.1840
             0.6350 0.0780 0.1840];
         
for ii=[4 6 7 8 5]
%% ALTITUDE DOWNRANGE
figure(1)
plot(files(ii).gr_trk,files(ii).alt,'LineWidth',1.2)
disp('Iron Ball 20 cm 7.6 km/s initial Zenith Ang. 88°')
grid on
hold on
ylim([0 max(files(ii).alt)])
ylabel('Altitude [km]')
xlabel('Dowrange [km]')
hleg = legend({'0°','28°','90°','97.5°','180°'},'Location','southwest');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

%% VEL PLOT
figure(2)
plot(files(ii).time,files(ii).vel/1000,'LineWidth',1.2)
disp('Iron Ball 20 cm 7.6 km/s initial Zenith Ang. 88°')
grid on
hold on
ylabel('Velocity [km/s]')
xlabel('Time mm:ss')
hleg = legend({'0°','28°','90°','97.5°','180°'},'Location','southwest');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

%% MASS PLOT
figure(3)
% subplot(2,1,1)
plot(files(ii).time,files(ii).mass,'LineWidth',1.2)
disp('Iron Ball 20 cm 7.6 km/s initial Zenith Ang. 88°')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Mass [Kg]')

% subplot(2,1,2)
% plot(files(ii).time(1:end-1),files(ii).mass_rate(1:end-1),'LineWidth',1.2)
% grid on
% hold on
% xlabel('Time mm:ss')
% ylabel('Mass loss rate [Kg/s]')
hleg = legend({'0°','28°','90°','97.5°','180°'},'Location','southwest');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

%% MAGNITUDE PLOT
figure(4)
colororder(newcolors)
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,1),'LineWidth',1.2)
disp('Iron Ball 20 cm 7.6 km/s initial Zenith Ang. 88°')
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Abs.Magnitude')
plot(files(ii).time(1:end-1),files(ii).mag(1:end-1,2),'--','LineWidth',1.2)
hleg = legend({'0° MAX','0° min','28° MAX','28° min','90° MAX','90° min','97.5° MAX','97.5° min','180° MAX','180° min'},'Location','northeastoutside');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

figure(5)
colororder(newcolors)
plot(files(ii).mag(1:end-1,1),files(ii).alt(1:end-1),'LineWidth',1.2)
disp('Iron Ball 20 cm 7.6 km/s initial Zenith Ang. 88°')
set(gca, 'XDir','reverse')
grid on
hold on
plot(files(ii).mag(1:end-1,2),files(ii).alt(1:end-1),'--','LineWidth',1.2)
ylim([0 100])
xlim([-10 10])
xlabel('Abs.Magnitude') 
ylabel('Altitude [km]')
hleg = legend({'0° MAX','0° min','28° MAX','28° min','90° MAX','90° min','97.5° MAX','97.5° min','180° MAX','180° min'},'Location','northeastoutside');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

%% Inclination Ang

figure(6)
plot(files(ii).time,90+files(ii).flight_ang,'LineWidth',1.2);
disp('Iron Ball 20 cm 7.6 km/s initial Zenith Ang. 88°')
grid on
hold on
xlabel('Time mm:ss') 
ylabel('Zenith Angle [°]')
hleg = legend({'0°','28°','90°','97.5°','180°'},'Location','southwest');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

%% EARTH PLOT
figure(7)
plot3(files(ii).lon,files(ii).lat,files(ii).alt,'LineWidth',1.2); 
disp('Iron Ball 20 cm 7.6 km/s initial Zenith Ang. 88°')
hold on
grid on
xlabel('Longitude [degrees]','Fontsize',10);
ylabel('Latitude [degrees]','Fontsize',10);
zlabel('Altitude [km]','Fontsize',10);
view([0 0 90])
hleg = legend({'0°','28°','90°','97.5°','180°'},'Location','northwest');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

end

figure(7)
lat = load('coast_lat.dat');
long = load('coast_long.dat');
plot3(long,lat,zeros(1,length(lat)),'k')
xlim([-180 180])
ylim([-90 90])
zlim([0 500])
set(gca,'Xtick',-180:30:180)
set(gca,'Ytick',-90:30:90)
view([0 0 90])

allChildren = get(hleg, 'Children');                % list of all objects on axes
displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
% Remove object associated with "data1" in legend
delete(allChildren(strcmp(displayNames, 'data1')))

hleg = legend({'0°','28°','90°','97.5°','180°'},'Location','southwest');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Incl.')

%% SAVE .jpg

% FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\Oct';   % Your destination folder
FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\3 Chap'; 
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    figure(iFig)
    set(gca,'FontSize',20)
    saveas(gcf,fullfile(FolderName, [num2str(iFig),'_SCARAB_diffIncl.jpg']))
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   savefig(fullfile(FolderName, [FigName '.fig']));
end