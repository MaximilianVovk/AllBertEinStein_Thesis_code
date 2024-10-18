clc
clear all
close all

%% DATA
fileroot='C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\SCARAB\';
files= dir(fullfile([append(fileroot,'Meteors')]));
%files=dir(fullfile('C:\Users\maxivy\Documents\DRAMA\TEST folder','*.txt'));


material_name=["iron" "basaltmax" "minbasalt" "ordchon" "carchon" "granite" "sandstone"];
material_density=[7870.0 2400.0 3100.0 3500.0 2800.0 2750.0 2000.0];
material_meltingheat=[272000.0 400000.0 506000.0 265000.0 265000.0 250000.0 680000.0];
material_color=["g" "b" "c" "r" "m" "k" "y"];
lum_eff = [0.1 0.001];

bmax5cm=[];
bmax10cm=[];
bmax20cm=[];  

newcolors = [0 0.4470 0.7410
             0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330];
 
r=10;
siz = [2 1 3]/2;
for ii=1:3
    f = figure(ii);
    f.WindowState = 'maximized';
%     subplot(2,4,8)


    f = figure(ii);
    f.WindowState = 'maximized';
    subplot(3,3,9)
%     subplot(2,4,8)
%     subplot(1,4,4)
    plot3(1,1,1,'b','Linewidth',1.2)
    hold on
    plot3(1,1,1,'b--','Linewidth',1.2)
    plot3(1,1,1,'r')
    plot3(1,1,1,'r--')
    [X,Y,Z] = sphere;
    s=surf(X * 1,Y * 1,Z * 1);
    xlim([-20*siz(ii) 20*siz(ii)])
    ylim([-20*siz(ii) 20*siz(ii)])
    zlim([-20*siz(ii) 20*siz(ii)])
%     xlim([-20 20])
%     ylim([-20 20])
%     zlim([-50 50])
    colormap([1 1 0]);
    legend('SCARAB Abs.Magnitude \tau = 0.1','SCARAB Abs.Magnitude \tau = 0.001','Theoretical Abs.Magnitude \tau = 0.1','Theoretical Abs.Magnitude \tau = 0.001','Initial Meteor Shape','location','north');
%     legend('SCARAB Abs.Magnitude \tau = 0.1','SCARAB Abs.Magnitude \tau = 0.001','Initial Meteor Shape','location','north');
    hold on
    axis off 

    subplot(3,3,8)
    [X,Y,Z] = sphere;
    s=surf(X * r(1),Y * r(1),Z * r(1));
    xlim([-40*siz(ii) 40*siz(ii)])
    ylim([-40*siz(ii) 40*siz(ii)])
    zlim([-20*siz(ii) 20*siz(ii)])
    colormap([1 1 0]);
    hold on
    axis off  
end

%% Run
matrix_holder=[];
load('meteorMass.mat')
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'1.83','match')) 
        %if ~isempty(regexp(files(ii).material,'basaltmax','match')) || ~isempty(regexp(files(ii).material,'ordchon','match')) || ~isempty(regexp(files(ii).material,'granite','match'))
%         run('config.m')
%         ast.init_v = files(ii).init_vel*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
%                                                 ast.v(1) = ast.init_v;
%         ast.init_zenith = 90-files(ii).init_ang;    %zenith angle
% 
%         % % %inclination to heading angle
%         i=[0.5751561E-01 28.24824 89.86656 97.39682 179.9425];
%         headANG=[0 30 93.5 101 180];
%         [p,S] = polyfit(i,headANG,3);
%         iMAX=[0 28 90 97.5 180];
%         head_ang = polyval(p,files(ii).init_inc);
%         % % ast.heading_ang=head_ang;
% 
%         ast.inc_orbit=head_ang;%files(ii).init_inc;%90
%         ast.lat(1)=0;
%         ast.lon(1)=0;
%                                                 ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
%                                                 ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));
%         ast.init_pos = r_planet + files(ii).init_alt*1000;
%                                                 ast.h(1) = ast.init_pos - r_planet;
% 
%         ast.r(1) = files(ii).init_size/100/2; %in meters
%         ast.rho = files(ii).rho;   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
%         ast.A = 1.2;    % shape factor (1.2 for sphere). %A=S/V^2/3 (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);%
%         ast.drag_coeff = 0.5;%0.5; % due to aircap actual Cd=2drag_coef
%         ast.heat_transfer = mean(files(ii).heat_transf_coef(files(ii).heat_transf_coef>0)); %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
%         ast.heat_of_ablation = files(ii).Hmelt;   %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6
% 
%         ast.init_mass = files(ii).mass(1);%4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
%                                                 ast.m(1)=ast.init_mass;
%         ast.lum_eff = [0.1 0.001];%0.001;%0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006 lum_eff = [0.1 0.001];
%                                                 ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki
% 
%                                                 ast.luminous_energy(1:2)= [nan nan];
%                                                 ast.luminous_energy_real(1:2)= [nan nan];
% 
%         ast.z_init= 10000;  %OBS_height(ii);%up                                        
%         ast.x_init= -9.3601e+04; %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
%         ast.y_init= -1.0177e+06; %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north
% 
%         ast.Up_pos(1)=(ast.h(1) - ast.z_init);
%         ast.Nort_pos(1)=(0) + ast.x_init; %(0) + ast.x_init;%
%         ast.West_pos(1)=(0) + ast.y_init; %(0) + ast.y_init;%
% 
%         %% CODE
% 
%         mvmax = 9.0;  %# faintest magnitude to be printed
%         tic
%         t=1; jj=1; jjt=0; ast.mv(1) = nan; ast.luminous_energy(1)=0;rho_aria(1)=rho_air;
%         h_dist(1)=sqrt((0 - obs.x_init)^2 + (ast.h(1) - obs.z_init)^2 + obs.y_init^2);  % distance respec to the observer
% 
%         %     # ----------------------------------------------------------------------------------------------
%         %     # do the computation
%         %     # ----------------------------------------------------------------------------------------------
%         while ast.h(t) > 0 && ast.m(t) > 0.001
%             [rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
%             rho_aria(t+1)=rho_air;
% 
%             mean_free_path(t)=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
%             Kn(t)=mean_free_path(t)/(ast.r(t)*2);
%             Cd(t)=Drag_coef(Kn(t),ast.A(t))/2;
%             ast.drag_coeff=Cd(t);
% 
% % %             ast.heat_transfer=mean(files(ii).heat_transf_coef(find(round(files(ii).alt)==round(ast.h(t)/1000))));
%         % %     ast.heat_transfer=mean(files(ii).coef_gauss(find(round(1:100)==round(ast.h(t)/1000))));
%             if isnan(ast.heat_transfer)
%                 ast.heat_transfer=0;
%             end
%             if isempty(ast.heat_transfer)
%                 ast.heat_transfer=0;
%             end
%             heatTransfer(t)=ast.heat_transfer;
% 
%             %om_earth=0.00007292;
%             om_earth=2*pi/(23.934472*3600);
%             ast.lat(t+1)=rad2deg( (ast.s(t)*sin(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lat(1);
%             ast.lon(t+1)=rad2deg( (ast.s(t)*cos(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lon(1);
% 
%             % calculate rates of change
%             a2 = ast.drag_coeff * ast.A(t) * rho_air * ast.v(t)^2 / (ast.m(t) * ast.rho^2)^0.33333; %ast.drag_coeff * 1/4 * 4*pi*ast.r^2 * rho_air * ast.v(t)^2 
%             gv = g0_planet / (1 + ast.h(t) / r_planet)^2;
% 
%         %     aW = -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.inc_orbit)) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) * cos(deg2rad(ast.inc_orbit)) + om_earth*sin(deg2rad(ast.lat(t)))*ast.v_horizontal(t) * sin(deg2rad(ast.inc_orbit));
%         %     aN = -a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.inc_orbit)) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) * sin(deg2rad(ast.inc_orbit)) + om_earth*sin(deg2rad(ast.lat(t)))*ast.v_horizontal(t) * cos(deg2rad(ast.inc_orbit));
%             ah = -a2 * ast.v_horizontal(t) / ast.v(t) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) + om_earth*sin(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t);
%             av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit));% Eötvös effect
% 
%             ast.zd(t+1) = rad2deg(atan2(ast.v_vertical(t), ast.v_horizontal(t)));
% 
%             ast.s(t+1) = ast.s(t) + ast.v_horizontal(t) * dt;
%             ast.h(t+1) = ast.h(t) + ast.v_vertical(t) * dt;
% 
%             ast.Up_pos(t+1)=(ast.h(t+1) - ast.z_init);
%             ast.Nort_pos(t+1)=  (ast.s(t+1)*cos(deg2rad(-ast.inc_orbit)) + ast.x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
%             ast.West_pos(t+1)=  (ast.s(t+1)*sin(deg2rad(ast.inc_orbit)) + ast.y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
%             h_dist(t+1) = sqrt((ast.West_pos(t+1))^2 + (ast.Nort_pos(t+1))^2 + (ast.Up_pos(t+1))^2 );  % distance respec to the observer
% 
%             ast.v_vertical(t+1) = ast.v_vertical(t) + av * dt;
%         %     ast.v_horizontal(t+1) = sqrt((ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit)) + aW * dt)^2 + (ast.v_horizontal(t)*sin(deg2rad(ast.inc_orbit)) + aN * dt)^2);
%             ast.v_horizontal(t+1) = ast.v_horizontal(t) + ah * dt;
%             ast.v(t+1) = sqrt((ast.v_horizontal(t))^2 + ast.v_vertical(t)^2); %-ast.Earth_atm_speed*ast.v(t)/ast.v(1)
% 
%         %% temperature and ablation fragmentatation
% 
%             ast.m_sec(t+1)=ast.A(t)*ast.heat_transfer/(2*ast.heat_of_ablation)*(ast.m(t)/ast.rho)^(2/3)*rho_air*ast.v(t)^3;%*abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+1);%abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+0.05);
%             ast.m(t+1)=ast.m(t)-(ast.m_sec(t+1)) * dt;
%             ast.r(t+1)=((ast.m(t)/ast.rho)*3/(4*pi))^(1/3);
% 
%             ast.A(t+1)=ast.A(t); 
% 
%         %% luminosity eff
%             ast.luminous_energy(t+1,1:2) = 0.5 * ast.v(t+1)^2 * ast.m_sec(t+1) .* ast.lum_eff ;%* 1e10 / (h_dist(t+1) * h_dist(t+1));
% 
%             if ast.luminous_energy(t+1,2) < 0.1
%                 ast.mv(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy(t+1,:));  %# magnitude before 2.086 but in -1.086 kampbel-brown & koshny
%             else
%                 ast.mv(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy(t+1,:));  %# magnitude before 2.086 but in -1.086 kampbel-brown & koshny
%                 if ast.mv(t+1) < mvmax
%                 jjt(jj)=t+1;
%                 jj=jj+1;
%                 end
%             end
% 
%             ast.luminous_energy_real(t+1,1:2) = 0.5 * ast.v(t+1)^2 * ast.m_sec(t+1) .* ast.lum_eff * 1e10 / (h_dist(t+1) * h_dist(t+1));
%             ast.mv_real(t+1,1:2) = 6.8 - 1.086 .* log(ast.luminous_energy_real(t+1,:));
%         %% cases
%             t = t+1;
% 
%              if toc > 20
%                  % stop simulation if it takes too much
%                  disp('too much time')
%                 break
%              end
%         end
%         
%         files(ii).AbsMagMean_Code= ast.mv;
%         files(ii).h_Code= ast.h;
%         files(ii).sec_Code=[0 1*dt:dt:(t-1)*dt];
        
        %% subplot
        if files(ii).init_size==5
            figure(3)
            sgtitle('Diameter  5cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
        elseif files(ii).init_size==10
            figure(1)
            sgtitle('Diameter 10cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
        elseif files(ii).init_size==20
            figure(2)
            sgtitle('Diameter 20cm 130 km 7.97 km/s Zenith Ang.=88.17° Orbit inc.=97.5°')
        end
        
        switch files(ii).material
            case "iron"
%                 subplot(2,4,1)
                subplot(3,3,1)
            case "basaltmax"
                subplot(3,3,2)
%                 subplot(2,4,2)
%                 subplot(1,4,1)
            case "minbasalt"
                subplot(3,3,3)
%                 subplot(2,4,3)
            case "ordchon"
                subplot(3,3,4)
%                 subplot(2,4,4)
%                 subplot(1,4,2)
            case "carchon"
                subplot(3,3,5)
%                 subplot(2,4,5)
            case "granite"
                subplot(3,3,6)
%                 subplot(2,4,6)
%                 subplot(1,4,3)
            case "sandstone"
                subplot(3,3,7)
%                 subplot(2,4,7)
        end
        
        plot(files(ii).mag(1:end-1,1),files(ii).alt(1:end-1),'b','Linewidth',1.2)
        ylim([0 100])
        xlim([-10 10])
        xlabel('Abs.Magnitude [-]') 
        ylabel('Altitude / km')
        set(gca, 'XDir','reverse')
        grid on
        hold on
        plot(files(ii).mag(1:end-1,2),files(ii).alt(1:end-1),'b--','Linewidth',1.2)
        plot(files(ii).AbsMagMean_Code(1:end-1,1),files(ii).h_Code(1:end-1)/1000,'r')
        plot(files(ii).AbsMagMean_Code(1:end-1,2),files(ii).h_Code(1:end-1)/1000,'r--')
        
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
        
        clear ast t
%         end
    end
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
    saveas(gcf,fullfile(FolderName, [' Meteor ',num2str(size),' cm SCARAB_mag #2.jpg']))
end

%save('meteorMass.mat','files')
%%
% [X,Y,Z] = sphere;
% r=1;
%  s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [0 0 1]);
%  hold on
%   s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [0 1 1]);
%    s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [0 1 0]);
%     s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [1 0 1]);
%      s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [34 139 34]/255);
%       s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [1 1 0]);
%        s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [1 0 0]);
%         s=surf(X * r(1),Y * r(1),Z * r(1), 'FaceColor', [0 0 0]);
% hleg = legend('500 km','450 km','400 km','350 km','300 km','250 km','200 km','170 km');
% htitle = get(hleg,'Title');
% set(htitle,'String','Detachment altitudes')
% 
%         
