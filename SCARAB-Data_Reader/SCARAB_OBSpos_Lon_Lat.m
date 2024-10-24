clc
clear all
close all
%% CHANGE FOR SIMULATION

sec_asteroid=300;%130
% sec_asteroid=180;%110 or 130 5
initial_altitude=130;
asteroidNum=10;

load('meteorRocketMass.mat')
load('meteorcageMass.mat')
load('meteorBallMass.mat')
load('meteorMissionMass.mat')


lum_eff = [0.1 0.001];

newcolors = [0 1 0
             0 1 0
             0 0 1
             0 0 1];
         
newcolors2= [0 1 0
             0 0 1
             0 1 1];

files(1)=rocket;
files(2)=ball;
% files(2)=files_mission(asteroidNum);
% files(3)=cage;

% lat_OBS=20;
% lon_OBS=-5;


figure(1)
colororder(newcolors2)
%     ii=1;
% lonRocket=-abs(files(ii).lon(files(ii).alt<initial_altitude)-mean(files(ii).lon(round(files(ii).alt)==initial_altitude)));
% lonRocket(lonRocket<-180)=-(360+lonRocket(lonRocket<-180));
% latRocket=-(files(ii).lat(files(ii).alt<initial_altitude)-mean(files(ii).lat(round(files(ii).alt)==initial_altitude)));
% plot3(lonRocket,latRocket,smooth(files(ii).alt(files(ii).alt<initial_altitude),20),'LineWidth',1.2); 
load('RocketSARA.mat')
plot3(saraRocket(:,4),saraRocket(:,3),saraRocket(:,2),'LineWidth',1.2); 
hold on
% view([0 0 90])
    ii=2;
plot3(files(ii).lon,files(ii).lat,files(ii).alt,'LineWidth',1.2); 

grid on
xlabel('Longitude [degrees]','Fontsize',10);
ylabel('Latitude [degrees]','Fontsize',10);
zlabel('Altitude [km]','Fontsize',10);
title('Meteor Rocket SCARAB Observers positions')
%Earth
xlim([-30 30])
ylim([-20 40])
% zlim([0 500])
set(gca,'Xtick',-180:30:180)
set(gca,'Ytick',-90:30:89)
% view([0 0 90])

FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\position_mag';


for lon_OBS=[2 1 0.5 0 -0.5 -1 -2]-5 %-5
    for lat_OBS=[2 1 0.5 0 -0.5 -1 -2]+23%[40 35 30 25 20 15 10 5 0] %21 rocket


        r_planet = 6371000;
                
        z_init= 10000;  %OBS_height(ii);%up 
        x_init= (z_init+r_planet)*deg2rad(0-lon_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
        y_init= (z_init+r_planet)*deg2rad(0-lat_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north
        
        figure(1)
        
        plot3(lon_OBS,lat_OBS,z_init/1000,'*r','LineWidth',1)   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
%         text(lon_OBS,lat_OBS,append(num2str(lat_OBS),'N ',num2str(lon_OBS),'W'),'FontSize',5)
        %% Rocket

        ii=1;

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
        DistanceRocket = sqrt((West_pos).^2 + (Nort_pos).^2 + (Up_pos).^2 );
        r4=[Nort_pos West_pos Up_pos];
        files(ii).azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
        files(ii).elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

        
        figure(2)
        colororder(newcolors)
        plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).REAL_mag(:,1),'g','LineWidth',1.2)
        hold on
        plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).REAL_mag(:,2),'g--','LineWidth',1.2)

        %% Meteor

        ii=2;

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

        figure(2)
        colororder(newcolors)
        time_consid=files(ii).time(files(ii).sec>sec_asteroid);
        RealMagMAX=files(ii).REAL_mag(files(ii).sec>sec_asteroid,1);
        RealMagmin=files(ii).REAL_mag(files(ii).sec>sec_asteroid,2);
        plot(time_consid(1:end-1),RealMagMAX(1:end-1),'b','LineWidth',1.2)
        plot(time_consid(1:end-1),RealMagmin(1:end-1),'b--','LineWidth',1.2)
        set(gca, 'YDir','reverse')
        grid on
        hold on
        xlabel('Time mm:ss')
        ylabel('Magnitude')
        ylim([-10 10])
        legend({'Rocket \tau = 0.1','Rocket \tau = 0.001','Meteor \tau = 0.1','Meteor \tau = 0.001'},'Location','southwest');
        
        titles=append(num2str(lat_OBS),'°N ',num2str(lon_OBS),'°W minimum Meteor Distance ',num2str(round(min(Distance)/1000)),' km');%files(ii).lon,files(ii).lat
        title(titles)
        subtitle(append('Rocket Landing Distance ',num2str(round(DistanceRocket(end)/1000)),' km Deg.Diff ',num2str(round(lat_OBS-latRocket(end),1)),'°N ',num2str(round(lon_OBS-lonRocket(end),1)),'°W'))
        
        figure(3)
        ii=1;
        skyplot(files(ii).azimuth(files(ii).elevation>0),files(ii).elevation(files(ii).elevation>0),'k.')
        hold on
        ii=2;
        skyplot(files(ii).azimuth(files(ii).elevation>0),files(ii).elevation(files(ii).elevation>0),'k.')
        skyplot(files(ii).azimuth(files(ii).sec>sec_asteroid),files(ii).elevation(files(ii).sec>sec_asteroid),'b.')
%         yy = (90-files(ii).elevation).*cos(files(ii).azimuth/180*pi);
%         xx = (90-files(ii).elevation).*sin(files(ii).azimuth/180*pi);
%         c=files(ii).REAL_lum_ene(:,2);
%         c(isnan(c))=0;
% %         if asteroidNum==10
%             c(files(ii).sec<sec_asteroid)=0;
% %         end
%         sz=5;
%         scatter(xx,yy,sz,c,'filled')
%         scatter(xx(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),yy(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),25,max(c),'filled')
%         colormap(winter)
% % skyplot(files(ii).azimuth(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),files(ii).elevation(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),'y*')
        ii=1;
        skyplot(files(ii).azimuth(~isnan(files(ii).REAL_mag(:,2))),files(ii).elevation(~isnan(files(ii).REAL_mag(:,2))),'g.')
        ii=2;
        g = figure(3);
%         colormap(g, flipud(colormap(g)))
        figure(3)
        saveas(gcf,fullfile(FolderName, ['Mag Azim Elev ',num2str(lat_OBS),'°N ',num2str(lon_OBS),'°W.jpg']))
        
        
        g=figure(2);
        axes('pos',[-.01 .395 .67 .51]) %[0.1300 0.1100 0.7750 0.8150]
        g.WindowState = 'maximized';
        imshow(fullfile(FolderName, 'The Border.jpg'));
        axes('pos',[.0 .4 .65 .5])
        imshow(fullfile(FolderName, ['Mag Azim Elev ',num2str(lat_OBS),'°N ',num2str(lon_OBS),'°W.jpg']));
        saveas(gcf,fullfile(FolderName, ['Mag Azim Elev ',num2str(lat_OBS),'°N ',num2str(lon_OBS),'°W.jpg']))
        
        close 2 3
    end
end

%% Plot EARTH
f = figure(1);
f.WindowState = 'maximized';
lat = load('coast_lat.dat');
long = load('coast_long.dat');
plot3(long,lat,zeros(1,length(lat)),'k')

ii=2; 
plot3(files(ii).lon(files(ii).sec<sec_asteroid),files(ii).lat(files(ii).sec<sec_asteroid),files(ii).alt(files(ii).sec<sec_asteroid),'k.')
plot3(saraRocket(~isnan(files(ii).REAL_mag(:,2)),4),saraRocket(~isnan(files(ii).REAL_mag(:,2)),3),saraRocket(~isnan(files(ii).REAL_mag(:,2)),2),'g.','LineWidth',1.2);
ii=1;
plot3(saraRocket(isnan(files(ii).REAL_mag(:,2)),4),saraRocket(isnan(files(ii).REAL_mag(:,2)),3),saraRocket(isnan(files(ii).REAL_mag(:,2)),2),'k.')
plot3(saraRocket(end-100:end,4),saraRocket(end-100:end,3),saraRocket(end-100:end,2),'k.')


hleg = legend({'Rocket','Meteor','Observer'},'Location','southwest');

allChildren = get(hleg, 'Children');                % list of all objects on axes
displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
% Remove object associated with "data1" in legend
delete(allChildren(strcmp(displayNames, 'data1')))

saveas(gcf,fullfile(FolderName, ['Meteor Rocket SCARAB Observers positions.jpg']))
