f=figure;
colororder(newcolors2)
plot3(saraRocket(:,4),saraRocket(:,3),saraRocket(:,2),'LineWidth',1.2); 
hold on

ii=1;
colororder(newcolors2)
lonRocket=-abs(files(ii).lon(files(ii).alt<initial_altitude)-mean(files(ii).lon(round(files(ii).alt)==initial_altitude)));
lonRocket(lonRocket<-180)=-(360+lonRocket(lonRocket<-180));
latRocket=-(files(ii).lat(files(ii).alt<initial_altitude)-mean(files(ii).lat(round(files(ii).alt)==initial_altitude)));
plot3(lonRocket,latRocket,smooth(files(ii).alt(files(ii).alt<initial_altitude),20),'LineWidth',1.2); 
hold on
grid on
xlabel('Longitude [degrees]','Fontsize',10);
ylabel('Latitude [degrees]','Fontsize',10);
zlabel('Altitude [km]','Fontsize',10);
title('Meteor Rocket SCARAB Observers positions')
%Earth
xlim([-40 40])
ylim([-40 40])
% zlim([0 500])
set(gca,'Xtick',-180:30:180)
set(gca,'Ytick',-90:30:89)

f.WindowState = 'maximized';
lat = load('coast_lat.dat');
long = load('coast_long.dat');
plot3(long,lat,zeros(1,length(lat)),'k')
