clc
clear all
close all
%% SIMULATION data
load('FRIPONmeteors.mat')

for imeteor=3:11
% imeteor=11;
    for OBS=1:length(filesFRIPON(imeteor).OBS)
        %% PLOT
        figure%(1)
        plot(filesFRIPON(imeteor).gr_trk_code{1,OBS},filesFRIPON(imeteor).alt_code{1,OBS},'r')%,ast.s(jjt)/1000,ast.h(jjt)/1000)
        title(filesFRIPON(imeteor).OBS{1,OBS})
        grid on
        hold on
        plot(filesFRIPON(imeteor).gr_trk{1,OBS},filesFRIPON(imeteor).alt{1,OBS},'b.')
        ylabel('Altitude [km]')
        xlabel('Dowrange [km]')
            legend({'Theory','Meteor'},'Location','southwest')

        figure%(2)
        plot(filesFRIPON(imeteor).time_code{1,OBS},filesFRIPON(imeteor).vel_code{1,OBS},'r')
        title(filesFRIPON(imeteor).OBS{1,OBS})
        grid on
        hold on
        plot(filesFRIPON(imeteor).time{1,OBS},filesFRIPON(imeteor).vel{1,OBS},'b.')
        ylabel('Velocity [km/s]')
        xlabel('Time mm:ss')
            legend({'Theory','Meteor'},'Location','southwest')

        figure%(3)
        plot(filesFRIPON(imeteor).time_code{1,OBS}(2:end),filesFRIPON(imeteor).ABSmag_code{1,OBS}(2:end,1),'r')
        title(filesFRIPON(imeteor).OBS{1,OBS})
        set(gca, 'YDir','reverse')
        grid on
        hold on
        plot(filesFRIPON(imeteor).time_code{1,OBS}(2:end),filesFRIPON(imeteor).ABSmag_code{1,OBS}(2:end,2),'r--')
        plot(filesFRIPON(imeteor).time{1,OBS},filesFRIPON(imeteor).ABSmag{1,OBS},'b.')
        xlabel('Time mm:ss')
        ylabel('Abs.Magnitude')
        legend({'MAX Theory','min Theory','Meteor'},'Location','southwest')

        figure%(4)
        plot(filesFRIPON(imeteor).ABSmag_code{1,OBS}(2:end,1),filesFRIPON(imeteor).alt_code{1,OBS}(2:end),'r','Linewidth',1.2);
        set(gca, 'XDir','reverse')
        grid on
        hold on
        title(filesFRIPON(imeteor).OBS{1,OBS})
        plot(filesFRIPON(imeteor).ABSmag_code{1,OBS}(2:end,2),filesFRIPON(imeteor).alt_code{1,OBS}(2:end),'r--','Linewidth',1.2);
        plot(filesFRIPON(imeteor).ABSmag{1,OBS},filesFRIPON(imeteor).alt{1,OBS},'b.')
        ylim([0 100])
        xlim([-10 10])
        xlabel('Abs.Magnitude') 
        ylabel('Altitude [km]')
        grid on
        hold on
        legend({'MAX Theory','min Theory','Meteor'},'Location','southwest')

        figure%(5)
        plot(filesFRIPON(imeteor).time_code{1,OBS}(2:end),filesFRIPON(imeteor).mag_code{1,OBS}(2:end,1),'r')
        title(filesFRIPON(imeteor).OBS{1,OBS})
        set(gca, 'YDir','reverse')
        grid on
        hold on
        plot(filesFRIPON(imeteor).time_code{1,OBS}(2:end),filesFRIPON(imeteor).mag_code{1,OBS}(2:end,2),'r--')
        plot(filesFRIPON(imeteor).time{1,OBS},filesFRIPON(imeteor).mag{1,OBS},'b.')
        xlabel('Time mm:ss')
        ylabel('Magnitude')
        legend({'MAX Theory','min Theory','Meteor'},'Location','southwest')

        figure%(6)
        plot(filesFRIPON(imeteor).mag_code{1,OBS}(2:end,1),filesFRIPON(imeteor).alt_code{1,OBS}(2:end),'r','Linewidth',1.2);
        set(gca, 'XDir','reverse')
        grid on
        hold on
        title(filesFRIPON(imeteor).OBS{1,OBS})
        plot(filesFRIPON(imeteor).mag_code{1,OBS}(2:end,2),filesFRIPON(imeteor).alt_code{1,OBS}(2:end),'r--','Linewidth',1.2);
        plot(filesFRIPON(imeteor).mag{1,OBS},filesFRIPON(imeteor).alt{1,OBS},'b.')
        ylim([0 100])
        xlim([-10 10])
        xlabel('Magnitude') 
        ylabel('Altitude [km]')
        grid on
        hold on
        legend({'MAX Theory','min Theory','Meteor'},'Location','southwest')

        
        
    end  
%% PLOT lon lat
    
    figure%(7)
    titlename=[];
%     plot3(filesFRIPON(imeteor).lon_code{1,OBS},filesFRIPON(imeteor).lat_code{1,OBS},filesFRIPON(imeteor).alt_code{1,OBS},'.r','LineWidth',1.2)  
    hold on
%     plot3(filesFRIPON(imeteor).lon{1,OBS},filesFRIPON(imeteor).lat{1,OBS},filesFRIPON(imeteor).alt{1,OBS},'.b','LineWidth',1.2) 
    for OBS=1:length(filesFRIPON(imeteor).OBS)
        plot3(filesFRIPON(imeteor).lon_code{1,OBS},filesFRIPON(imeteor).lat_code{1,OBS},filesFRIPON(imeteor).alt_code{1,OBS},'.r','LineWidth',1.2)  

        plot3(filesFRIPON(imeteor).OBS_lon(1,OBS),filesFRIPON(imeteor).OBS_lat(1,OBS),filesFRIPON(imeteor).OBS_alt(1,OBS),'*m','LineWidth',1)   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
        titlename=[titlename filesFRIPON(imeteor).OBS{1,OBS} ' '];
        text(filesFRIPON(imeteor).OBS_lon(1,OBS),filesFRIPON(imeteor).OBS_lat(1,OBS),filesFRIPON(imeteor).OBS{1,OBS})
    end
    title( append(titlename) )
    grid on
    lat = load('coast_lat.dat');
    long = load('coast_long.dat');
    plot3(long,lat,zeros(1,length(lat)),'k')
    xlim([min(filesFRIPON(imeteor).lon{1,OBS})-5 max(filesFRIPON(imeteor).lon{1,OBS})+5])
    ylim([min(filesFRIPON(imeteor).lat{1,OBS})-5 max(filesFRIPON(imeteor).lat{1,OBS})+5])
    zlim([0 500])
    set(gca,'Xtick',-180:5:180)
    set(gca,'Ytick',-90:5:90)
    xlabel('Longitude [degrees]','Fontsize',10);
    ylabel('Latitude [degrees]','Fontsize',10);
    zlabel('Altitude [km]','Fontsize',10);
    view([0 0 90])
%     legend({'Theory','Meteor',filesFRIPON(imeteor).OBS{1,OBS}},'Location','southwest')
%     legend({'Theory',filesFRIPON(imeteor).OBS{1,OBS}},'Location','southwest')

    %% PLOT SKY
    figure
    for OBS=1:length(filesFRIPON(imeteor).OBS)
        skyplot(filesFRIPON(imeteor).azimuth_code{1,OBS}(2:end),filesFRIPON(imeteor).elevation_code{1,OBS}(2:end),'.r');
        hold on
        
        % transform data to Cartesian coordinates.
        yy = (90-filesFRIPON(imeteor).elevation_code{1,OBS}(end)).*cos(filesFRIPON(imeteor).azimuth_code{1,OBS}(end)/180*pi);
        xx = (90-filesFRIPON(imeteor).elevation_code{1,OBS}(end)).*sin(filesFRIPON(imeteor).azimuth_code{1,OBS}(end)/180*pi);
        
        text(xx,yy,filesFRIPON(imeteor).OBS{1,OBS},'FontSize',8)

    end
    title( append(titlename) )
end

%% SAVE .jpg

FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\FRIPON_Image';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    figure(iFig)
    saveas(gcf,fullfile(FolderName, ['FRIPON_',filesFRIPON(imeteor).OBS{1,1},'_',filesFRIPON(imeteor).OBS{1,2},'_',num2str(iFig),'_Meteors.jpg']))
end