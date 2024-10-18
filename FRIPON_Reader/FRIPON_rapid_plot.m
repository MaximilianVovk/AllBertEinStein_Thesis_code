clc
clear all
close all
%% SIMULATION
load('FRIPONmeteors.mat')
run('config.m')

imeteor=10; % 3<ii<11
OBS=1; %number of the observatory
rhoii=1; %1 more mass 6 less mass
     %1.20000000000000	0.0350000000000000	1000000	1500
minMass=0.000000001;
    %% VARIABLE DATA
    if isempty(filesFRIPON(imeteor).A)
        ast.A(t) = 1.2;    % shape factor (1.2 for sphere). %A=S/V^2/3 (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);%
%             more mv t
        ast.heat_transfer = 0.035; %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
%             less mv t
        ast.heat_of_ablation = 1000000; %2 106 8 106  %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6
%STONY&IRON heat of fusion 265000 heat of melt 6000000
    else
        ast.A(t)=filesFRIPON(imeteor).A; %0.3
        ast.heat_transfer=filesFRIPON(imeteor).heat_tran;
        ast.heat_of_ablation=filesFRIPON(imeteor).heat_ablat;
    end
%     
%     ast.A(t) = 1.6; %1.2 0.3 1.5
%     ast.heat_transfer = 0.2; %0.002 0.0175 0.5
%     ast.heat_of_ablation = 4000000; % 100000 2000000 melt - 8000000 vap 8000000

    %% DATA
    ii=imeteor;
                                        filesFRIPON(imeteor).A=ast.A(t);
                                        filesFRIPON(imeteor).heat_tran=ast.heat_transfer;
                                        filesFRIPON(imeteor).heat_ablat=ast.heat_of_ablation;
                                        
    ast.init_v = mean(filesFRIPON(imeteor).vel{1,OBS})*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                                        ast.v(1) = ast.init_v; %10.519535315283466*1000; for 7OBS1 %10.519535315283466*1000;%
                 filesFRIPON(imeteor).init_vel(1,OBS)=ast.init_v/1000;
                 VELOCITY=ast.init_v/1000
        
        switch imeteor
            case 3
                ast.init_zenith = 91-mean(filesFRIPON(imeteor).incl_ang{1,OBS});
                minMass=1.46344918961e-151;
            case 4
                ast.init_zenith = 90-mean(filesFRIPON(imeteor).incl_ang{1,OBS});
                minMass=2.17414763837e-12;
            case 5
                ast.init_zenith = 101-mean(filesFRIPON(imeteor).incl_ang{1,OBS});
                minMass=6.57663478005e-08;
            case 6
                ast.init_zenith = 96-mean(filesFRIPON(imeteor).incl_ang{1,OBS});
                minMass=3.42778744896e-06; %1.82544893732e-07;
            case 7
                ast.init_zenith = 94-mean(filesFRIPON(imeteor).incl_ang{1,OBS}); %94-67.731045774778140  95.5-67.731045774778140;%
                minMass=0.14630148555;
            case 8
                 ast.init_zenith = 94-mean(filesFRIPON(imeteor).incl_ang{1,OBS});     
                 minMass=7.42247835131e-11; %2.15207952198e-11;
            case 9
                 ast.init_zenith = 99-mean(filesFRIPON(imeteor).incl_ang{1,OBS}); 
                 minMass=0.0198691859557;
            case 10
                ast.init_zenith = 98-mean(filesFRIPON(imeteor).incl_ang{1,OBS});
                minMass=2.90350304292e-196;
            case 11
                ast.init_zenith = 93-mean(filesFRIPON(imeteor).incl_ang{1,OBS}); %zenith angle 98/10 93/11 94-95/7 99/9 94/8 96/6 100/5 90/4 91/3
                minMass=8.35715115021e-26;
            otherwise
                error('not in mapped')
        end
                                        ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith)); % 95-67.731045774778140; for 7OBS1
                                        ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));

    ast.init_pos = r_planet + filesFRIPON(imeteor).alt{1,OBS}(1)*1000;
                                        ast.h(1) = filesFRIPON(imeteor).alt{1,OBS}(1)*1000;
    ast.inc_orbit=mean(filesFRIPON(imeteor).head_ang{1,OBS}(1));%files(ii).init_inc;%90
    ast.lat(1)=filesFRIPON(imeteor).lat{1,OBS}(1);
    ast.lon(1)=filesFRIPON(imeteor).lon{1,OBS}(1);

    ast.z_init= filesFRIPON(imeteor).OBS_alt(1,OBS);%up                                        
    ast.x_init= (ast.z_init+r_planet)*deg2rad(ast.lon(1)-filesFRIPON(imeteor).OBS_lon(1,OBS));%west
    ast.y_init= (ast.z_init+r_planet)*deg2rad(ast.lat(1)-filesFRIPON(imeteor).OBS_lat(1,OBS));%north
    
    ast.luminous_energy(1:2)= [nan nan];
    ast.luminous_energy_real(1:2)= [nan nan];
    ast.lum_eff = [0.1 0.001];%0.001;%0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006 lum_eff = [0.1 0.001];
    ast.drag_coeff = 0.5;%0.5; % due to aircap actual Cd=2drag_coef

    ast.Up_pos(1)=(ast.h(1) - ast.z_init);
    ast.Nort_pos(1)=(ast.h(t)+r_planet)*deg2rad(ast.lon(t)-filesFRIPON(imeteor).OBS_lon(1,OBS));%(ast.h(1)+r_planet)*deg2rad(ast.lon(1)-filesFRIPON(ii).OBS_lon(1,OBS)); %(0) + ast.x_init;%
    ast.West_pos(1)=(ast.h(t)+r_planet)*deg2rad(ast.lat(t)-filesFRIPON(imeteor).OBS_lat(1,OBS));%(ast.h(1)+r_planet)*deg2rad(ast.lat(1)-filesFRIPON(ii).OBS_lat(1,OBS)); %(0) + ast.y_init;%
    
    if isempty(filesFRIPON(imeteor).rho_code)
        ast.r(1) = filesFRIPON(imeteor).size(rhoii); %in meters
        ast.rho = filesFRIPON(imeteor).rho(rhoii);   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
        ast.init_mass = filesFRIPON(imeteor).mass(rhoii);
        
        filesFRIPON(imeteor).rho_code=ast.rho;
        filesFRIPON(imeteor).mass_code=ast.init_mass;  
        filesFRIPON(imeteor).size_code=ast.r;
        
    else      
        ast.r(1) = filesFRIPON(imeteor).size_code; %in meters
        ast.rho = filesFRIPON(imeteor).rho_code;   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
        ast.init_mass = filesFRIPON(imeteor).mass_code;

    end
                                            ast.m(1)= ast.init_mass;
                                            ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki

    %% CODE

    mvmax = 9.0;  %# faintest magnitude to be printed
    tic
    t=1; jjt=0; ast.mv(1) = nan; ast.luminous_energy(1)=0;rho_aria(1)=rho_air;
    h_dist(1)=sqrt((ast.s(1)*cos(deg2rad(ast.inc_orbit)) - ast.x_init)^2 + (ast.s(1)*sin(deg2rad(ast.inc_orbit)) - ast.y_init)^2 + (ast.h(1) - ast.z_init)^2 );  % distance respec to the observer
ast.h(1) = filesFRIPON(imeteor).alt{1,OBS}(1)*1000;

    ast.v_side(t)=0; ast.side(t)=0; init_inc=ast.inc_orbit; ast.AnghSide(t)=0;
    %     # ----------------------------------------------------------------------------------------------
    %     # do the computation
    %     # ----------------------------------------------------------------------------------------------
% ast.mv(1)=0;
while ast.h(t) > 0 && (ast.m(t) > minMass)

    [rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
        rho_aria(t+1)=rho_air;

        mean_free_path=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
        Kn(t)=mean_free_path/(ast.r(t)*2);
        Cd(t)=Drag_coef(Kn(t),ast.A(t))/2;
        ast.drag_coeff=Cd(t); 

        %om_earth=0.00007292;
        om_earth=2*pi/(23.934472*3600);
        ast.lat(t+1)=-rad2deg( (ast.s(t)*sin(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) + ast.lat(1);
        ast.lon(t+1)=-rad2deg( (ast.s(t)*cos(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) + ast.lon(1);

        % calculate rates of change
        a2 = ast.drag_coeff * ast.A(t) * rho_air * ast.v(t)^2 / (ast.m(t) * ast.rho^2)^0.33333; %ast.drag_coeff * 1/4 * 4*pi*ast.r^2 * rho_air * ast.v(t)^2 
        gv = g0_planet / (1 + ast.h(t) / r_planet)^2;

        ah = -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
        av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
        aSide = -a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.AnghSide(t))) + 2*om_earth * ast.v_horizontal(t)* sin(deg2rad(ast.lat(t))) - 2*om_earth * ast.v_vertical(t) * cos(deg2rad(ast.lat(t))) * sin(deg2rad(ast.inc_orbit));
    
        ast.zd(t+1) = rad2deg(atan2(ast.v_vertical(t), ast.v_horizontal(t)));
        ast.AnghSide(t+1) = rad2deg(atan2(ast.side(t),ast.s(t)));
        ast.inc_orbit=init_inc- ast.AnghSide(t+1);
        control_inc(t)=ast.inc_orbit;
    
        ast.s(t+1) = ast.s(t) + ast.v_horizontal(t) * dt;
        ast.h(t+1) = ast.h(t) + ast.v_vertical(t) * dt;
        ast.side(t+1) = ast.side(t) + ast.v_side(t) * dt;
    
        ast.Up_pos(t+1)=(ast.h(t+1) - ast.z_init);
        ast.Nort_pos(t+1)=  (ast.s(t+1)*cos(deg2rad(-ast.inc_orbit)) + ast.side(t+1)*cos(deg2rad(ast.inc_orbit)) + ast.x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
        ast.West_pos(t+1)=  (ast.s(t+1)*sin(deg2rad(ast.inc_orbit)) + ast.side(t+1)*sin(deg2rad(ast.inc_orbit)) + ast.y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
        h_dist(t+1) = sqrt((ast.West_pos(t+1))^2 + (ast.Nort_pos(t+1))^2 + (ast.Up_pos(t+1))^2 );  % distance respec to the observer

        ast.v_vertical(t+1) = ast.v_vertical(t) + av * dt;
        ast.v_side(t+1) = ast.v_side(t) + aSide * dt;
        ast.v_horizontal(t+1) = ast.v_horizontal(t) + ah * dt;
        ast.v(t+1) = sqrt((ast.v_horizontal(t))^2 + ast.v_vertical(t)^2 + ast.v_side(t)); %-ast.Earth_atm_speed*ast.v(t)/ast.v(1)

    %% temperature and ablation fragmentatation

        ast.m_sec(t+1)=ast.A(t)*ast.heat_transfer/(2*ast.heat_of_ablation)*(ast.m(t)/ast.rho)^(2/3)*rho_air*ast.v(t)^3;%*abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+1);%abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+0.05);
        ast.m(t+1)=ast.m(t)-(ast.m_sec(t+1)) * dt;
        ast.r(t+1)=((ast.m(t)/ast.rho)*3/(4*pi))^(1/3);

        ast.A(t+1)=ast.A(t);         

    %% luminosity eff
        ast.luminous_energy(t+1,1:2) = 0.5 * ast.v(t+1)^2 * ast.m_sec(t+1) .* ast.lum_eff ;%* 1e10 / (h_dist(t+1) * h_dist(t+1));
        ast.mv(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy(t+1,:));  %# magnitude before 2.086 but in -1.086 kampbel-brown & koshny

        ast.luminous_energy_real(t+1,1:2) = 0.5 * ast.v(t+1)^2 * ast.m_sec(t+1) .* ast.lum_eff * 1e10 / (h_dist(t+1) * h_dist(t+1));
        ast.mv_real(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy_real(t+1,:));
        %% cases
        t = t+1;

         if toc > 20
             % stop simulation if it takes too much
             disp('too much time')
            break
         end
    end
    if ast.m(end) < 0.001
        ast.m(end)=0;
    end

    thetaa = mod((ast.s / r_planet),(2 * pi));
    timeSec=[0 1*dt:dt:(t-1)*dt]+filesFRIPON(imeteor).sec{1,OBS}(1);
    sec = seconds(timeSec);
    sec.Format = 'mm:ss.SSS';

    %% SAVE DATA
    
    r4=[ast.West_pos' ast.Nort_pos' ast.Up_pos'];
    azimuth_ast=rad2deg(atan2(r4(:,2),r4(:,1)));
    elevation_ast=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));
    
    
    filesFRIPON(imeteor).sec_code{1,OBS} = timeSec';
    filesFRIPON(imeteor).time_code{1,OBS} = sec';
    filesFRIPON(imeteor).alt_code{1,OBS} = ast.h'/1000;
    filesFRIPON(imeteor).gr_trk_code{1,OBS} = ast.s'/1000; % downrange_code
    filesFRIPON(imeteor).zenithang_code{1,OBS} = ast.zd';
    filesFRIPON(imeteor).vel_code{1,OBS} = ast.v'/1000;
    filesFRIPON(imeteor).distOBS_code{1,OBS} = h_dist'/1000;
    filesFRIPON(imeteor).mag_code{1,OBS} = ast.mv_real;
    filesFRIPON(imeteor).ABSmag_code{1,OBS} = ast.mv;
    filesFRIPON(imeteor).A_code{1,OBS} = ast.A';
    filesFRIPON(imeteor).azimuth_code{1,OBS} = azimuth_ast;
    filesFRIPON(imeteor).elevation_code{1,OBS} = elevation_ast;
    filesFRIPON(imeteor).lon_code{1,OBS} = ast.lon';
    filesFRIPON(imeteor).lat_code{1,OBS} = ast.lat';

    %% PLOT

    figure(1)
    plot(filesFRIPON(imeteor).gr_trk_code{1,OBS},filesFRIPON(imeteor).alt_code{1,OBS},'r')%,ast.s(jjt)/1000,ast.h(jjt)/1000)
title(filesFRIPON(imeteor).OBS{1,OBS})
    grid on
    hold on
    plot(filesFRIPON(imeteor).gr_trk{1,OBS},filesFRIPON(imeteor).alt{1,OBS},'b.')
    ylabel('Altitude [km]')
    xlabel('Dowrange [km]')
        legend({'Theory','Meteor'},'Location','southwest')

    figure(2)
    plot(filesFRIPON(imeteor).time_code{1,OBS},filesFRIPON(imeteor).vel_code{1,OBS},'r')
title(filesFRIPON(imeteor).OBS{1,OBS})
    grid on
    hold on
    plot(filesFRIPON(imeteor).time{1,OBS},filesFRIPON(imeteor).vel{1,OBS},'b.')
    ylabel('Velocity [km/s]')
    xlabel('Time mm:ss')
        legend({'Theory','Meteor'},'Location','southwest')

        figure(3)
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
        
    figure(4)
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
        
    figure(5)
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

    figure(6)
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
        
figure(7)
    plot3(filesFRIPON(imeteor).lon_code{1,OBS},filesFRIPON(imeteor).lat_code{1,OBS},filesFRIPON(imeteor).alt_code{1,OBS},'.r','LineWidth',1.2)  
    hold on
    plot3(filesFRIPON(imeteor).lon{1,OBS},filesFRIPON(imeteor).lat{1,OBS},filesFRIPON(imeteor).alt{1,OBS},'.b','LineWidth',1.2) 
    plot3(filesFRIPON(imeteor).OBS_lon(1,OBS),filesFRIPON(imeteor).OBS_lat(1,OBS),filesFRIPON(imeteor).OBS_alt(1,OBS),'*m','LineWidth',1)   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
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
        legend({'Theory','Meteor',filesFRIPON(imeteor).OBS{1,OBS}},'Location','southwest')


        

%% PLOT SKY

figure
skyplot(filesFRIPON(imeteor).azimuth_code{1,OBS}(2:end),filesFRIPON(imeteor).elevation_code{1,OBS}(2:end),'.r')
hold on
title(filesFRIPON(imeteor).OBS{1,OBS})

save('FRIPONmeteors.mat','filesFRIPON')

figure(6)