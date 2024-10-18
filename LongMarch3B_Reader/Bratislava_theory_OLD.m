clc
clear all
close all
%% DATA
ii=4;
for ii=1:4
run('config.m')
% load('meteorMass.mat')

load('sat1_6.3kms.mat')
load('sat2_6.7kms.mat')
load('sat3_7.1kms.mat')
load('sat4_7.5kms.mat')

allSat_fragm{1}=sat1;
allSat_fragm{2}=sat2;
allSat_fragm{3}=sat3;
allSat_fragm{4}=sat4;

    ast.init_v = allSat_fragm{ii}(1,5)*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                                            ast.v(1) = ast.init_v;
    ast.init_zenith = 90-(360-allSat_fragm{ii}(1,7));    %zenith angle
                                            ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
                                            ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));
    ast.init_pos = r_planet + allSat_fragm{ii}(1,2)*1000;
                                            ast.h(1) = ast.init_pos - r_planet;
obs.x_init = -100000;
obs.y_init = 100000;
obs.z_init = 10000.0;

    if ii==1
        ast.r(1) = 3.5/100; %in meters
    elseif ii==2
        ast.r(1) = 1.7/100; %in meters
    elseif ii==3
        ast.r(1) = 2/100; %in meters
    elseif ii==4
        ast.r(1) = 3/100; %in meters
    end
    ast.rho = 2813;%7870;   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
    ast.A = 1.2;    % shape factor (1.2 for sphere). %A=S/V^2/3 (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);%
    if ii==2
        ast.A = 7.5; % shape factor
        L=(pi * (ast.r(1))^2)^(3/2)/(4/3 * pi * (ast.r(1))^2 * ast.A^(3/2));
    elseif ii==4
        ast.A = 4; % shape factor for rotation 10
        L=(pi * (ast.r(1))^2)^(3/2)/(4/3 * pi * (ast.r(1))^2 * ast.A^(3/2));
%         R=ast.r(1);
%         Afin=(pi * R^2)/(4/3 * pi * R^2*L)^(2/3);
    end
    
    ast.drag_coeff = 0.5;%0.5; % due to aircap actual Cd=2drag_coef
%     heatCoefGauss=Gauss_HeatCoef(1:0.1:100,ast.init_v/1000,ast.r(1)*100,"iron"); %iron ordchon granite
%     heat_transfer_m=max(heatCoefGauss);
    ast.heat_transfer = 0.6;%heat_transfer_m*2; %0.6; %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
    if ii==2
        ast.heat_transfer = 0.35; % heat_transfer
    end
    if ii==3
        ast.heat_transfer = 0.8; % heat_transfer
    end
    if ii==4
        ast.heat_transfer = 0.65; %rotation0.57% heat_transfer
    end

    ast.heat_of_ablation = 400000;%272000;   %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6

                                            ast.init_mass = 4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
                                            ast.m(1)=ast.init_mass;
    ast.lum_eff = [0.1 0.001];%0.001;%0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006 lum_eff = [0.1 0.001];
                                            ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki
    
    ast.dir_ang= -45; %yaw to determine its position 
    %% CODE

    mvmax = 9.0;  %# faintest magnitude to be printed

    tic
    t=1;
    jj=1;
    jjt=0;
    h_dist(1)=sqrt((cos(deg2rad(ast.dir_ang))*0 - obs.x_init)^2 + (sin(deg2rad(ast.dir_ang))*0 - obs.y_init)^2 + (ast.h(1) - obs.z_init)^2 );  % distance respec to the observer
    ast.mv(1) = nan;
    ast.luminous_energy(1)=0;
    rho_aria(1)=rho_air;

    % ast.heat_of_ablation = ;
    % ast.heat_transfer = mean();

    %     # ----------------------------------------------------------------------------------------------
    %     # do the computation
    %     # ----------------------------------------------------------------------------------------------
    while ast.h(t) > 0 && ast.m(t) > 0.001
        [rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
        rho_aria(t+1)=rho_air;

        mean_free_path=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
        Kn(t)=mean_free_path/(ast.r(t)*2);
        Cd(t)=Drag_coef(Kn(t),ast.A(t))/2;
%         if ii==2 
%             Cd(t)=(Drag_coef(Kn(t),ast.A(t)))/2; % shape factor
%         elseif ii==4
%             Cd(t)=(Drag_coef(Kn(t),ast.A(t)))/2;
%         end
        ast.drag_coeff=Cd(t);

    %     ast.heat_transfer=mean(heatCoefGauss([1:0.1:100]*10==round(ast.h(t)/100)));
    %     if isnan(ast.heat_transfer)
    %         ast.heat_transfer=0;
    %     end

        % calculate rates of change
        a2 = ast.drag_coeff * ast.A(t) * rho_air * ast.v(t)^2 / (ast.m(t) * ast.rho^2)^0.33333; %ast.drag_coeff * 1/4 * 4*pi*ast.r^2 * rho_air * ast.v(t)^2
        gv = g0_planet / (1 + ast.h(t) / r_planet)^2;
        av = -gv - a2 * ast.v_vertical(t) / ast.v(t) +ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t));
        ah = -a2 * ast.v_horizontal(t) / ast.v(t) - ast.v_vertical(t) * ast.v_horizontal(t) / (r_planet + ast.h(t));
        ast.zd(t+1) = rad2deg(atan2(ast.v_vertical(t), ast.v_horizontal(t)));

        ast.s(t+1) = ast.s(t) + ast.v_horizontal(t) * dt * r_planet / (r_planet + ast.h(t));
        ast.h(t+1) = ast.h(t) + ast.v_vertical(t) * dt;
        h_dist(t+1) = sqrt((cos(deg2rad(ast.dir_ang))*ast.s(t+1) - obs.x_init)^2 + (sin(deg2rad(ast.dir_ang))*ast.s(t+1) - obs.y_init)^2 + (ast.h(t+1) - obs.z_init)^2 );  % distance respec to the observer

        ast.v_vertical(t+1) = ast.v_vertical(t) + av * dt;
        ast.v_horizontal(t+1) = ast.v_horizontal(t) + ah * dt;
        ast.v(t+1) = sqrt(ast.v_horizontal(t)^2 + ast.v_vertical(t)^2);

    %% temperature and ablation fragmentatation

        ast.m_sec(t+1)=ast.A(t)*ast.heat_transfer/(2*ast.heat_of_ablation)*(ast.m(t)/ast.rho)^(2/3)*rho_air*ast.v(t)^3;%*abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+1);%abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+0.05);
        ast.m(t+1)=ast.m(t)-(ast.m_sec(t+1)) * dt;
        ast.r(t+1)=((ast.m(t)/ast.rho)*3/(4*pi))^(1/3);

        ast.A(t+1)=ast.A(t); 
        
        if ii==1
            R=ast.r(1);
            L=((ast.m(t)/ast.rho)*3/(4*pi*R^2));
            ast.A(t+1)=(pi * R^2)/(4/3 * pi * R^2*L)^(2/3);
        end
        
        if ii==2
            R=ast.r(t+1);
            L=L;
            ast.A(t+1)=(pi * R^2)/(4/3 * pi * R^2*L)^(2/3);
        end
        
        if ii==4
            R=ast.r(t+1);
            L=L;
            ast.A(t+1)=(pi * R^2)/(4/3 * pi * R^2*L)^(2/3);
%             if ast.A(t)>=Afin
%                 reduce=1;
%             elseif ast.A(t)<0.2
%                 reduce=0;
%             end
%             
%             if reduce==1
%                 ast.A(t+1)=ast.A(t)-0.1;
%             else
%                 ast.A(t+1)=ast.A(t)+0.1;
%             end
%             Afin=(pi * R^2)/(4/3 * pi * R^2*L)^(2/3);
        end

    %% luminosity eff
%     ast.h_dist(t+1)=1000*mean(allSat_fragm{ii}(( round(allSat_fragm{ii}(:,2)) == round(ast.h(t)/1000) ),8));    %     ast.heat_transfer=mean(heatCoefGauss([1:0.1:100]*10==round(ast.h(t)/100)));
%     if isnan(ast.h_dist(t+1))
%         ast.h_dist(t+1)=ast.h_dist(t);
%     end
        ast.luminous_energy(t+1,1:2) = 0.5 * ast.v(t+1)^2 * ast.m_sec(t+1) .* ast.lum_eff ;%* 1e10 / (h_dist(t+1) * h_dist(t+1));
        ast.luminous_energy_OBS(t+1,1:2) = 0.5 * ast.v(t+1)^2 * ast.m_sec(t+1) .* ast.lum_eff * 1e10 / (h_dist(t+1) * h_dist(t+1));

    %                     lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff ;
    %                     lum_eneAvg(lum_eneAvg==0)=nan;
    %                     magAvg = 6.8 - 1.086 * log(lum_eneAvg);
        if ast.luminous_energy(t+1,2) < 0.1
            ast.mv(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy(t+1,:));  %# magnitude before 2.086 but in -1.086 kampbel-brown & koshny
            ast.mv_OBS(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy_OBS(t+1,:));
        else
            ast.mv(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy(t+1,:));  %# magnitude before 2.086 but in -1.086 kampbel-brown & koshny
            ast.mv_OBS(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy_OBS(t+1,:));
            if ast.mv(t+1) < mvmax
            jjt(jj)=t+1;
            jj=jj+1;
            end
        end
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
    sec = seconds(1*dt:dt:t*dt);
    sec.Format = 'mm:ss'; 
    sec_obs=seconds(allSat_fragm{ii}(:,1));
    sec_obs.Format='mm:ss';

    %% PLOT
    figure(1)
    subplot(2,2,ii)
    plot(ast.s/1000,(ast.h)/1000,'r')%,ast.s(jjt)/1000,ast.h(jjt)/1000)
    title(['Fragment ' num2str(ii)])
    grid on
    hold on
    plot(allSat_fragm{ii}(:,4),allSat_fragm{ii}(:,2),'b.')
    %plot(obs.x_init/1000,obs.z_init/1000,'xr','linewidth',2)
    ylabel('Altitude [km]')
    xlabel('Dowrange [km]')
    xlim([0 (ast.s(end)+100)/1000])
    legend({'Theory','CZ-3B R/B'},'Location','southwest')

    figure(2)
    subplot(2,2,ii)
    plot(sec,ast.v/1000,'r')
    title(['Fragment ' num2str(ii)])
    grid on
    hold on
    plot(sec_obs,allSat_fragm{ii}(:,5),'b.')
    % plot(sec(jjt),ast.v(jjt)/1000)
    ylabel('Velocity [km/s]')
    xlabel('Time mm:ss')
    legend({'Theory','CZ-3B R/B'},'Location','southwest')

    figure(3)
    subplot(2,2,ii)
    % plot(sec,ast.mv)
    % plot(sec,ast.mv)
    %plot([sec'; flip(sec'); sec(1)],[ast.mv(:,2); flip(ast.mv(:,1)); ast.mv(2,2)],'b')
    plot(sec(2:end)',ast.mv(2:end,1),'r')
    title(['Fragment ' num2str(ii)])
    set(gca, 'YDir','reverse')
    grid on
    hold on
    plot(sec(2:end)',ast.mv(2:end,2),'r--')
    plot(sec_obs,allSat_fragm{ii}(:,3),'b.')
    % plot([files(ii).time(1:end-1); flip(files(ii).time(1:end-1)); files(ii).time(1)],[files(ii).mag(1:end-1,2); flip(files(ii).mag(1:end-1,1)); files(ii).mag(1,2)],'r')
    xlabel('Time mm:ss')
    ylabel('Magnitude')
    % plot(sec(jjt),ast.mv(jjt))
    legend({'MAX Theory','min Theory','CZ-3B R/B'},'Location','southwest')

    figure(4)
    subplot(2,2,ii)
    plot(ast.mv(2:end,1),ast.h(2:end)/1000,'r','Linewidth',1.2);
    % plot([ast.mv(:,2); flip(ast.mv(:,1)); ast.mv(1,2)],[ast.h/1000 flip(ast.h/1000) ast.h(1)/1000],'b')
    % title('magnitude')
    set(gca, 'XDir','reverse')
    grid on
    hold on
    title(['Fragment ' num2str(ii)])
    plot(ast.mv(2:end,2),ast.h(2:end)/1000,'r--','Linewidth',1.2);
    plot(allSat_fragm{ii}(:,3),allSat_fragm{ii}(:,2),'b.')
    % plot(ast.mv(jjt),ast.h(jjt)/1000)
    % plot([files(ii).mag(1:end-1,2); flip(files(ii).mag(1:end-1,1)); files(ii).mag(2,2)],[files(ii).alt(1:end-1); flip(files(ii).alt(1:end-1)); files(ii).alt(1)],'r')
    ylim([0 100])
    xlim([-10 10])
    xlabel('Magnitude') 
    ylabel('Altitude [km]')
    % set(gca, 'XDir','reverse')
    grid on
    hold on
    legend({'MAX Theory','min Theory','CZ-3B R/B'},'Location','southwest')
    
%     figure
%     plot(ast.h_dist(2:end)/1000,ast.mv_OBS(2:end,1),'r')
%     title(['Fragment ' num2str(ii)])
%     set(gca, 'YDir','reverse')
%     grid on
%     hold on
%     plot(ast.h_dist(2:end)/1000,ast.mv_OBS(2:end,2),'r--')
%     plot(allSat_fragm{ii}(:,8),allSat_fragm{ii}(:,6),'b.')
%     xlabel('Distance [km]')
%     ylabel('Magnitude')
%     legend({'MAX Theory','min Theory','CZ-3B R/B'},'Location','southwest')
% 
    disp(ast.m(1))
    disp(ast.m(end))
    
%     x_mot=cos(deg2rad(ast.dir_ang)).*ast.s - obs.x_init;
%     y_mot=sin(deg2rad(ast.dir_ang)).*ast.s - obs.y_init;
%     z_mot=abs(ast.h-obs.z_init);
%     elev=rad2deg(atan(z_mot./(vecnorm([x_mot' y_mot'],2,2))' ));
%     azim=rad2deg(atan(x_mot./y_mot));
    figure(5)
    skyplot(allSat_fragm{ii}(:,9)+180,allSat_fragm{ii}(:,10),'.')
    hold on
%     skyplot(azim,elev,'.')

    clear all
end
