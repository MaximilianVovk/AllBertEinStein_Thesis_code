clc
clear all
close all

%% DATA

material_name=["iron" "basaltmax" "minbasalt" "ordchon" "carchon" "granite" "sandstone"];
material_density=[7870.0 2400.0 3100.0 3500.0 2800.0 2750.0 2000.0];
material_meltingheat=[272000.0 400000.0 506000.0 265000.0 265000.0 250000.0 680000.0];
material_color=["g" "b" "c" "r" "m" "k" "y"];

lum_eff = [0.1 0.001];

sec_asteroid=300;%130
initial_altitude=130;

lum_eff = [0.1 0.001];

r_planet = 6371000;
                                                                               
lat_OBS=23; %best 23 or closer to the second stage 22
lon_OBS=-5; %best -5

z_init= 10000;  %OBS_height(ii);%up 
x_init= (z_init+r_planet)*deg2rad(0-lon_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
y_init= (z_init+r_planet)*deg2rad(0-lat_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north

bmax5cm=[];
bmax10cm=[];
bmax20cm=[];

newcolors = [0 0.4470 0.7410
             0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330];
 
% r = [10 20 5]/2;
r=10;
siz = [2 1 3]/2;
for ii=1:3
    f = figure(ii);
    f.WindowState = 'maximized';
    subplot(2,4,8)

    [X,Y,Z] = sphere;
    s=surf(X * r(1),Y * r(1),Z * r(1));
    xlim([-20*siz(ii) 20*siz(ii)])
    ylim([-20*siz(ii) 20*siz(ii)])
    zlim([-40*siz(ii) 40*siz(ii)])

    col=colorbar('westoutside');
    colormap(flipud(winter))
%     colormap(f, flipud(colormap(f)));
    col.Label.String ='\tau = 0.1 Magnitude';
    col.Limits=([-10 10]);
    col.Direction='reverse';
    colormap(f, flipud(colormap(f)))

    hold on
    axis off  
end

%% Run
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
    meteorAzAngrocket=deg2rad(files(1).azimuth(1:length(files(2).time)))-deg2rad(files(2).azimuth); % distance
    check_az=[(files(1).azimuth(1:length(files(2).time))) (files(2).azimuth)];
    
    Separat_Ang=rad2deg(sqrt(rocketElevDist.^2+meteorElevDist.^2-2*rocketElevDist.*meteorElevDist.*cos(meteorAzAngrocket)));
        
    allang=[rad2deg(rocketElevDist) rad2deg(meteorElevDist) rad2deg(meteorAzAngrocket) Separat_Ang];
    
    if asteroidNum==10
        files(ii).REAL_mag(files(ii).sec<sec_asteroid,1)=nan;
        files(ii).REAL_mag(files(ii).sec<sec_asteroid,2)=nan;
    end
    
        %% subplot
        if files(ii).init_size==5
            g =figure(3);
            sgtitle('Diameter  5cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
        elseif files(ii).init_size==10
            g =figure(1);
            sgtitle('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
        elseif files(ii).init_size==20
            g =figure(2);
            sgtitle('Diameter 20cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
        end
        
        switch files(ii).material
            case "iron"
                subplot(2,4,1)
                colormap(g, flipud(colormap(g)))
            case "basaltmax"
                subplot(2,4,2)
                colormap(g, flipud(colormap(g)))
%                 subplot(1,4,1)
            case "minbasalt"
                subplot(2,4,3)
                colormap(g, flipud(colormap(g)))
            case "ordchon"
                subplot(2,4,4)
                                colormap(g, flipud(colormap(g)))
%                 subplot(1,4,2)
            case "carchon"
                subplot(2,4,5)
                                colormap(g, flipud(colormap(g)))
            case "granite"
                subplot(2,4,6)
%                 subplot(1,4,3)
                colormap(g, flipud(colormap(g)))
            case "sandstone"
                subplot(2,4,7)
                               colormap(g, flipud(colormap(g)))
        end

        sz = 5; % size of the points they can vary if use a vector
        c=files(ii).REAL_mag(:,1);
        c(isnan(c))=20;
        scatter(files(2).time,Separat_Ang,sz,c,'filled')
        grid on
        hold on
        scatter(files(2).time(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),Separat_Ang(files(ii).REAL_mag(:,2)==min(files(ii).REAL_mag(:,2))),25,c(files(ii).REAL_mag(:,1)==min(files(ii).REAL_mag(:,1))),'filled')

        colormap(flipud(winter))
        clim([-10 10])
        set(gca,'yscale','log')
        xlabel('Time mm:ss')
        ylabel('Separation Angle [°] ') %separation angle \theta
%         title(titles, 'FontSize', fontSize)
%         h = colorbar;
%         ylabel(h, '\tau = 0.1 Magnitude')
%         colormap(g, flipud(colormap(g)))
%         clim([-10 10])
        
        switch files(ii).material
            case "iron"
                title('Iron')
            case "basaltmax"
                title('Basalt Acidic')
            case "minbasalt"
                title('Basalt Basic')
            case "ordchon"
                title('Ordinary Chondrite')
            case "carchon"
                title('Carbonaceous Chondrite')
            case "granite"
                title('Granite')
            case "sandstone"
                title('Sandstone')
        end
        
        clear files g

end

FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\Oct';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
        if iFig==3
            size=5;
        elseif iFig==1
            size=10;
        elseif iFig==2
            size=20;
        end
    f = figure(iFig);
    saveas(gcf,fullfile(FolderName, [' Meteor ',num2str(size),' cm SCARAB distance material.jpg']))
end