clc
clear all
close all
%% INPUT

% fragm_r=[3.5 2.1 1.8 5 3.5 2.1 1.8 5]./100; %3.5 1.7 2 3
% fragm_A= [1.2 3 1.2 3 1.2 3 1.2 3]; %1.2 7.5 1.2 4
% fragm_transf= [0.45 0.35 0.8 0.6 0.45 0.35 0.8 0.6]; %0.6 0.35 0.8 0.65

% fragm_r=[3.5 1 1.5 3 3.5 2.1 1.8 5]./100; %3.5 1.7 2 3
% fragm_A= [1.2 3 2 3 1.2 3 1.2 3]; %1.2 7.5 1.2 4
% fragm_transf= [0.45 0.35 0.8 0.6 0.45 0.35 0.8 0.6]; %0.6 0.35 0.8 0.65

% fragm_r=[3.5 3 2.5 5.5 3.5 3 2.5 5.5]./100; %3.5 1.7 2 3
% fragm_A= [1.2 3 2.2 3 1.2 3 2.2 3]; %1.2 7.5 1.2 4
% fragm_transf= [0.45 0.8 0.8 0.95 0.45 0.8 0.8 0.95]; %0.6 0.35 0.8 0.65

fragm_r=[3.5 6 7 11 3.5 6 7 11]./100; %3.5 1.7 2 3
fragm_A= [1.2 6 7 6.5 1.2 6 7 6.5]; %1.2 7.5 1.2 4
fragm_transf= [0.4 0.6 0.5 0.65 0.4 0.6 0.5 0.65]; %0.6 0.35 0.8 0.65

fragm_rho= [2813 2813 2813 2813 2813 2813 2813 2813];
fragm_melt= [400000 400000 400000 400000 400000 400000 400000 400000];
OBS_lat= [19.823108 19.823108 19.823108 19.823108 20.707521 20.707521 20.707521 20.707521]; %Maunakea and Haleaana
OBS_lon= [-155.473366 -155.473366 -155.473366 -155.473366 -156.256904 -156.256904 -156.256904 -156.256904];
OBS_height= [4207 4207 4207 4207 3055 3055 3055 3055]; % mountain altitude

%% DATA

newcolors2= [0.3010 0.7450 0.9330
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880];

load('LongMarch3B_Fragments.mat')

for ii=1:8
    run('config.m')
    
    if LongMarch3B(ii).sec(1) ~= 0
        ast.init_v = (LongMarch3B(ii+4).vel(1))*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                                            ast.v(1) = ast.init_v;
        ast.init_zenith = 90-(360-LongMarch3B(ii+4).ang(1)); %zenith angle
                                            ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
                                            ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));

        ast.init_pos = r_planet + LongMarch3B(ii+4).alt(1)*1000;
        
        ast.inc_orbit=mean(LongMarch3B(ii+4).head_ang);%files(ii).init_inc;%90
        ast.lat(1)=LongMarch3B(ii+4).lat(1);
        ast.lon(1)=LongMarch3B(ii+4).lon(1);

        ast.init_pos = r_planet + LongMarch3B(ii+4).alt(1)*1000;
                                                ast.h(1) = ast.init_pos - r_planet;
    else
        ast.init_v = (LongMarch3B(ii).vel(1))*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                                            ast.v(1) = ast.init_v;
        ast.init_zenith = 90-(360-LongMarch3B(ii).ang(1)); %zenith angle
                                            ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
                                            ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));

        ast.init_pos = r_planet + LongMarch3B(ii).alt(1)*1000;
        
        ast.inc_orbit=mean(LongMarch3B(ii).head_ang);%files(ii).init_inc;%90
        ast.lat(1)=LongMarch3B(ii).lat(1);
        ast.lon(1)=LongMarch3B(ii).lon(1);

        ast.init_pos = r_planet + LongMarch3B(ii).alt(1)*1000;
                                                ast.h(1) = ast.init_pos - r_planet;
    end
    
    ast.z_init= OBS_height(ii);%up                                        
    ast.x_init= (ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
    ast.y_init= (ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north
    
    ast.luminous_energy(1:2)= [nan nan];
    ast.luminous_energy_real(1:2)= [nan nan];
    ast.lum_eff = [0.1 0.001];%0.001;%0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006 lum_eff = [0.1 0.001];
    ast.drag_coeff = 0.5;%0.5; % due to aircap actual Cd=2drag_coef

    ast.Up_pos(1)=(ast.h(1) - ast.z_init);
    ast.Nort_pos(1)=(ast.h(1)+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii)); %(0) + ast.x_init;%
    ast.West_pos(1)=(ast.h(1)+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii)); %(0) + ast.y_init;%
    %% VARIABLE DATA
    ast.r(1) = fragm_r(ii); %in meters
    ast.rho = fragm_rho(ii);   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
    ast.A = fragm_A(ii);    % shape factor (1.2 for sphere). %A=S/V^2/3 (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);%
    ast.heat_transfer = fragm_transf(ii); %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
    ast.heat_of_ablation = fragm_melt(ii);   %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6

    if ast.A~=1.2
        L=(pi * (ast.r(1))^2)^(3/2)/(4/3 * pi * (ast.r(1))^2 * ast.A^(3/2));
        ast.init_mass = 4/3 * pi * ast.r^2 * L * ast.rho;%4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
    else
        ast.init_mass = 4/3 * pi * ast.r^3 * ast.rho;%4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
    end
                                            ast.m(1)=ast.init_mass;
                                            massss=ast.init_mass
                                            ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki

    %% CODE

    mvmax = 9.0;  %# faintest magnitude to be printed
    tic
    t=1; jj=1; jjt=0; ast.mv(1) = nan; ast.luminous_energy(1)=0;rho_aria(1)=rho_air;
    h_dist(1)=sqrt((ast.s(1)*cos(deg2rad(ast.inc_orbit)) - ast.x_init)^2 + (ast.s(1)*sin(deg2rad(ast.inc_orbit)) - ast.y_init)^2 + (ast.h(1) - ast.z_init)^2 );  % distance respec to the observer

    ast.v_side(t)=0; ast.side(t)=0; init_inc=ast.inc_orbit; ast.AnghSide(t)=0;
    %     # ----------------------------------------------------------------------------------------------
    %     # do the computation
    %     # ----------------------------------------------------------------------------------------------
    while ast.h(t) > 0 && ast.m(t) > ast.init_mass*0.01
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

        if (ii==1 || ii==5)
            ah =  -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
            av =  4-gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
        elseif (ii==2 || ii==6)
            ah =  -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
            av =  11-gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect

        elseif (ii==3 || ii==7)
            ah =  -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
            av =  9-gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
        elseif (ii==4 || ii==8)
            ah =  -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
            av =  11-gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect

        else
            ah =  -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
            av =  -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
            
        end
        aSide = -a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.AnghSide(t))) + 2*om_earth * ast.v_horizontal(t)* sin(deg2rad(ast.lat(t))) - 2*om_earth * ast.v_vertical(t) * cos(deg2rad(ast.lat(t))) * sin(deg2rad(ast.inc_orbit));
        
%         ah2 = -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t));
%         av2 = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t));


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
        
        if (ii==1 || ii==5) && t*dt>30
            ast.A(t+1)=ast.A(t)+0.05;
        end
        
%         if (ii==2 || ii==6) && (ii==3 || ii==7) && (ii==4 || ii==8)
%             R=ast.r(t+1);
%             L=L;
%             ast.A(t+1)=(pi * R^2)/(4/3 * pi * R^2*L)^(2/3);
%         end
        
%         if (ii==4 || ii==8)
%             ast.A(t+1)=ast.A(1)+sin((t)/40);
%         end
        

    %% luminosity eff
        ast.luminous_energy(t+1,1:2) = 0.5 * ast.v(t+1)^2 * ast.m_sec(t+1) .* ast.lum_eff ;%* 1e10 / (h_dist(t+1) * h_dist(t+1));

        if ast.luminous_energy(t+1,2) < 0.1
            ast.mv(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy(t+1,:));  %# magnitude before 2.086 but in -1.086 kampbel-brown & koshny
        else
            ast.mv(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy(t+1,:));  %# magnitude before 2.086 but in -1.086 kampbel-brown & koshny
            if ast.mv(t+1) < mvmax
            jjt(jj)=t+1;
            jj=jj+1;
            end
        end

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
%     sec = seconds(1*dt:dt:t*dt);
    sec = seconds([0 1*dt:dt:(t-1)*dt]);
    sec.Format = 'mm:ss';

    %% SAVE DATA
    timeSec=1*dt:dt:t*dt;
    r4=[ast.West_pos' ast.Nort_pos' ast.Up_pos'];
    azimuth_ast=rad2deg(atan2(r4(:,2),r4(:,1)));
    elevation_ast=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));
    
    LongMarch3B(ii).sec_code = timeSec';
    LongMarch3B(ii).time_code = sec';
    LongMarch3B(ii).alt_code = ast.h'/1000;
    LongMarch3B(ii).gr_trk_code = ast.s'/1000; % downrange_code
    LongMarch3B(ii).zenithang_code = ast.zd';
    LongMarch3B(ii).vel_code = ast.v'/1000;
    LongMarch3B(ii).distOBS_code = h_dist'/1000;
    LongMarch3B(ii).mag_code = ast.mv_real;
    LongMarch3B(ii).ABSmag_code = ast.mv;
    LongMarch3B(ii).A_code = ast.A';
    LongMarch3B(ii).azimuth_code = azimuth_ast;
    LongMarch3B(ii).elevation_code = elevation_ast;
    LongMarch3B(ii).lon_code = ast.lon';
    LongMarch3B(ii).lat_code = ast.lat';

    %% PLOT
    numfig=7;
    
    numTimes=0;
    if ii>=5
        numTimes=1;
    end
    
    figure(1+(numTimes*numfig))
    subplot(2,2,ii-(numTimes*4))
    plot(LongMarch3B(ii).gr_trk_code,LongMarch3B(ii).alt_code,'r')%,ast.s(jjt)/1000,ast.h(jjt)/1000)
    title(['Fragment 0' num2str(ii-(numTimes*4))])
    grid on
    hold on
    xlim([0 310])
%     xlim([0 320])
%     ylim([LongMarch3B(ii).alt_code(end)-1 LongMarch3B(ii).alt_code(1)])
    plot(LongMarch3B(ii).gr_trk,LongMarch3B(ii).alt,'b.')
    ylabel('Altitude / km')
    xlabel('Downrange / km')
%     xlim([LongMarch3B(ii).gr_trk_code(1) (LongMarch3B(ii).gr_trk_code(end)+100)])
    if ii==3 || ii==7
        legend({'Theory','CZ-3B R/B'},'Location','southwest')
    end
    if ii>=5        
        sgtitle('$Haleakal\bar{a}$ Observatory','interpreter','latex')
    else
        sgtitle('$Maunakea$ Observatory','interpreter','latex')
    end

    figure(2+(numTimes*numfig))
    subplot(2,2,ii-(numTimes*4))
    plot(LongMarch3B(ii).time_code,LongMarch3B(ii).vel_code,'r')
    title(['Fragment 0' num2str(ii-(numTimes*4))])
    grid on
    hold on
    xlim([0 LongMarch3B(3).time(end)])
    plot(LongMarch3B(ii).time,LongMarch3B(ii).vel,'b.')
    ylabel('Velocity [km/s]')
    xlabel('Time mm:ss')
    if ii==3 || ii==7
        legend({'Theory','CZ-3B R/B'},'Location','southwest')
    end
    if ii>=5        
        sgtitle('$Haleakal\bar{a}$ Observatory','interpreter','latex')
    else
        sgtitle('$Maunakea$ Observatory','interpreter','latex')
    end

    figure(3+(numTimes*numfig))
    subplot(2,2,ii-(numTimes*4))
    plot(LongMarch3B(ii).time_code(2:end),LongMarch3B(ii).ABSmag_code(2:end,1),'r')
    title(['Fragment 0' num2str(ii-(numTimes*4))])
    set(gca, 'YDir','reverse')
    grid on
    hold on
    xlim([0 LongMarch3B(3).time(end)])
    plot(LongMarch3B(ii).time_code(2:end),LongMarch3B(ii).ABSmag_code(2:end,2),'r--')
    plot(LongMarch3B(ii).time,LongMarch3B(ii).ABSmag,'b.')
    xlabel('Time mm:ss')
    ylabel('Abs.Magnitude')
    if ii==1 || ii==5
        legend({'\tau = 0.1','\tau = 0.001','CZ-3B R/B'},'Location','northwest')
    end
    if ii>=5        
        sgtitle('$Haleakal\bar{a}$ Observatory','interpreter','latex')
    else
        sgtitle('$Maunakea$ Observatory','interpreter','latex')
    end

    figure(4+(numTimes*numfig))
    subplot(2,2,ii-(numTimes*4))
    plot(LongMarch3B(ii).ABSmag_code(2:end,1),LongMarch3B(ii).alt_code(2:end),'r','Linewidth',1.2);
    set(gca, 'XDir','reverse')
    grid on
    hold on
    title(['Fragment 0' num2str(ii-(numTimes*4))])
    plot(LongMarch3B(ii).ABSmag_code(2:end,2),LongMarch3B(ii).alt_code(2:end),'r--','Linewidth',1.2);
    plot(LongMarch3B(ii).ABSmag,LongMarch3B(ii).alt,'b.')
    ylim([0 100])
    xlim([-10 10])
    xlabel('Abs.Magnitude') 
    ylabel('Altitude [km]')
    grid on
    hold on
    if ii==3 || ii==7
        legend({'\tau = 0.1 Theory','\tau = 0.001 Theory','CZ-3B R/B'},'Location','southwest')
    end
    if ii>=5        
        sgtitle('$Haleakal\bar{a}$ Observatory','interpreter','latex')
    else
        sgtitle('$Maunakea$ Observatory','interpreter','latex')
    end

    figure(5+(numTimes*numfig))
    subplot(2,2,ii-(numTimes*4))
    plot(LongMarch3B(ii).time_code(2:end),LongMarch3B(ii).mag_code(2:end,1),'r')
    title(['Fragment 0' num2str(ii-(numTimes*4))])
    set(gca, 'YDir','reverse')
    grid on
    hold on
    xlim([0 LongMarch3B(3).time(end)])
    plot(LongMarch3B(ii).time_code(2:end),LongMarch3B(ii).mag_code(2:end,2),'r--')
    plot(LongMarch3B(ii).time,LongMarch3B(ii).mag,'b.')
    xlabel('Time mm:ss')
    ylabel('Magnitude')
    if ii==1 || ii==5
        legend({'\tau = 0.1 Theory','\tau = 0.001 Theory','CZ-3B R/B'},'Location','northwest')
    end
    if ii>=5        
        sgtitle('$Haleakal\bar{a}$ Observatory','interpreter','latex')
    else
        sgtitle('$Maunakea$ Observatory','interpreter','latex')
    end
    
    figure(6+(numTimes*numfig))
    subplot(2,2,ii-(numTimes*4))
    plot(LongMarch3B(ii).mag_code(2:end,1),LongMarch3B(ii).alt_code(2:end),'r','Linewidth',1.2);
    set(gca, 'XDir','reverse')
    grid on
    hold on
    title(['Fragment 0' num2str(ii-(numTimes*4))])
    plot(LongMarch3B(ii).mag_code(2:end,2),LongMarch3B(ii).alt_code(2:end),'r--','Linewidth',1.2);
    plot(LongMarch3B(ii).mag,LongMarch3B(ii).alt,'b.') 
    ylim([0 100])
    xlim([-10 10])
    xlabel('Magnitude') 
    ylabel('Altitude [km]')
    grid on
    hold on
    if ii==3 || ii==7
        legend({'\tau = 0.1 Theory','\tau = 0.001 Theory','CZ-3B R/B'},'Location','southwest')
    end
    if ii>=5        
        sgtitle('$Haleakal\bar{a}$ Observatory','interpreter','latex')
    else
        sgtitle('$Maunakea$ Observatory','interpreter','latex')
    end

%     figure(7+(numTimes*numfig))
%     subplot(2,2,ii-(numTimes*4))
if (numTimes==1)
    figure(7)
    colororder(newcolors2)
%     plot3(LongMarch3B(ii).lon_code,LongMarch3B(ii).lat_code,LongMarch3B(ii).alt_code,'.r','LineWidth',1.2)
%     
    hold on
    le(ii-4)=plot3([LongMarch3B(ii).lon;LongMarch3B(ii-4).lon],[LongMarch3B(ii).lat;LongMarch3B(ii-4).lat],[LongMarch3B(ii).alt;LongMarch3B(ii-4).alt],'.','LineWidth',1.2,'DisplayName',append('Fragment 0',num2str(ii-4)) ) ;
    grid on


    if ii==8
        le(ii-4+1)=plot3(OBS_lon(ii),OBS_lat(ii),OBS_height(ii)/1000,'*m','LineWidth',1);   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
        le(ii-4+2)=plot3(OBS_lon(ii-4),OBS_lat(ii-4),OBS_height(ii)/1000,'*m','LineWidth',1);   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west

        text(OBS_lon(ii)-0.1,OBS_lat(ii)-0.3,'$Haleakal\bar{a}$ Observatory','interpreter','latex')
        text(OBS_lon(ii)-0.1,OBS_lat(ii-4)-0.3,'$Maunakea$ Observatory','interpreter','latex')

        lat = load('coast_lat.dat');
        long = load('coast_long.dat');
        le(ii-4+3)=plot3(long,lat,zeros(1,length(lat)),'k');
        xlim([-163 -153])
        ylim([15 25])
    %     zlim([0 500])
        set(gca,'Xtick',-180:5:180)
        set(gca,'Ytick',-90:5:90)
        xlabel('Longitude [degrees]','Fontsize',10);
        ylabel('Latitude [degrees]','Fontsize',10);
        zlabel('Altitude [km]','Fontsize',10);
    
        view([0 0 90])

        hleg = legend([le(1) le(2) le(3) le(4)],'location','southwest');

%         allChildren = get(hleg, 'Children');                % list of all objects on axes
%         displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
%         % Remove object associated with "data1" in legend
%         delete(allChildren(strcmp(displayNames, 'data1')))
    end
end
%     if ii==4 || ii==8
%         legend({'Theory','CZ-3B R/B','Observatory'},'Location','southwest')
%     end
%     if ii>=5        
%         sgtitle('$Haleakal\bar{a}$ Observatory','interpreter','latex')
%     else
%         sgtitle('$Maunakea$ Observatory','interpreter','latex')
%     end
    
    if ii==8
    else
        clear t ast
    end

end

%% PLOT SKY

figure
leggendsSKY(1)=skyplot(0,90,'k.');
hold on
colororder(newcolors2)
for ii=1:4

%     skyplot(LongMarch3B(ii).azimuth+180,LongMarch3B(ii).elevationORaltitude,'.')
%     hold on

    yy = (90-LongMarch3B(ii).elevationORaltitude).*cos((LongMarch3B(ii).azimuth+180)/180*pi);
    xx = (90-LongMarch3B(ii).elevationORaltitude).*sin((LongMarch3B(ii).azimuth+180)/180*pi);
    leggendsSKY(ii+1)=plot(xx,yy,'.','DisplayName', append('Fragment 0',num2str(ii) ) );

%     skyplot(LongMarch3B(ii).azimuth_code,LongMarch3B(ii).elevation_code,'.')
%     hold on
end
% legend(leggendsSKY(2:ii+1),'location','southwest');
% title('$Maunakea$ Observatory','interpreter','latex')
text(mean(xx),mean(yy)-5,'$Maunakea$ Observatory','interpreter','latex')

% figure
% leggendsSKY2(1)=skyplot(0,90,'k.');
% hold on
for ii=5:8
%     skyplot(LongMarch3B(ii).azimuth+180,LongMarch3B(ii).elevationORaltitude,'.')
%     hold on

    yy = (90-LongMarch3B(ii).elevationORaltitude).*cos((LongMarch3B(ii).azimuth+180)/180*pi);
    xx = (90-LongMarch3B(ii).elevationORaltitude).*sin((LongMarch3B(ii).azimuth+180)/180*pi);
    leggendsSKY(ii+1)=plot(xx,yy,'.','DisplayName', append('Fragment 0',num2str(ii-4) ) );

%     skyplot(LongMarch3B(ii).azimuth_code,LongMarch3B(ii).elevation_code,'.')
%     hold on
end
text(mean(xx),mean(yy)-10,'$Haleakal\bar{a}$ Observatory','interpreter','latex')
legend(leggendsSKY(6:ii+1),'location','southwest'); %,'NumColumns',2
% title('$Haleakal\bar{a}$ Observatory','interpreter','latex')

figure(7)
legend(leggendsSKY(6:ii+1),'location','southwest'); %,'NumColumns',2
    
save('LongMarch3B_Fragments.mat','LongMarch3B')
%% SAVE .jpg

FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\AMOS results';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    figure(iFig)
    saveas(gcf,fullfile(FolderName, ['AMOS_',num2str(iFig+13),'_Meteors.jpg']))
end