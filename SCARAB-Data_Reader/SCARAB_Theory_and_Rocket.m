clc
clear all
close all
%% CHANGE FOR SIMULATION

sec_asteroid=0;%130
% sec_asteroid=180;%110 or 130 5
initial_altitude=130;
asteroidNum=10;

load('meteorRocketMass.mat')
load('meteorcageMass.mat')
load('meteorBallMass.mat')
load('meteorMissionMass.mat')

files(1)=rocket;
files(2)=files_mission(asteroidNum);

titles=append(files(2).material,' Meteor ',num2str(files(2).init_size),' cm');

%% METEOR DATA
ast.r(1) = 0.1; %in meters
ast.rho = 2800;   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
ast.A = 1.2;    % shape factor (1.2 for sphere). %A=S/V^2/3 (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);%
ast.drag_coeff = 0.5;%0.5; % due to aircap actual Cd=2drag_coef
ast.heat_transfer = 0.02; %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
ast.heat_of_ablation = 500000;   %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6

%% DATA

load('meteorRocketMass.mat')
load('meteorcageMass.mat')
load('meteorBallMass.mat')
load('meteorMissionMass.mat')

lum_eff = [0.1 0.001];

newcolors = [0 1 0
             0 1 0
             0 0 1
             0 0 1
             0 1 1
             0 1 1];
         
newcolors2= [0 1 0
             0 0 1
             0 1 1];

r_planet = 6371000;
                                                                               
lat_OBS=20;
lon_OBS=-5;

z_init= 10000;  %OBS_height(ii);%up 
x_init= (z_init+r_planet)*deg2rad(0-lon_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
y_init= (z_init+r_planet)*deg2rad(0-lat_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north

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
r4=[Nort_pos West_pos Up_pos];
files(ii).azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
files(ii).elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

figure(1)
subplot(2,1,1)
colororder(newcolors2)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),smooth(files(ii).vel(files(ii).alt<initial_altitude)/1000,20),'LineWidth',1.2)
hold on

subplot(2,1,2)
colororder(newcolors)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mag((files(ii).alt<initial_altitude),1),'LineWidth',1.2)
hold on
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).mag((files(ii).alt<initial_altitude),2),'--','LineWidth',1.2)


figure(2)
colororder(newcolors)
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).REAL_mag(:,1),'LineWidth',1.2)
hold on
plot(files(ii).time(files(ii).alt<initial_altitude)-mean(files(ii).time(round(files(ii).alt)==initial_altitude)),files(ii).REAL_mag(:,2),'--','LineWidth',1.2)

figure(3)
colororder(newcolors2)
% plot3(lonRocket,latRocket,files(ii).alt(files(ii).alt<initial_altitude),'LineWidth',1.2); 
% load('RocketSARA.mat')
plot3(saraRocket(:,4),saraRocket(:,3),saraRocket(:,2),'LineWidth',1.2); 
hold on
% view([0 0 90])


figure(4)
skyplot(files(ii).azimuth,files(ii).elevation,'g.')
hold on

%% asteroid SCARAB

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
    
figure(1)
subplot(2,1,1)
colororder(newcolors2)
plot(files(ii).time,files(ii).vel/1000,'b','LineWidth',1.2)
grid on
xlim([files(ii).time(1) files(ii).time(end)])

subplot(2,1,2)
colororder(newcolors)
time_consid=files(ii).time(files(ii).sec>sec_asteroid);
MagMAX=files(ii).mag(files(ii).sec>sec_asteroid,1);
Magmin=files(ii).mag(files(ii).sec>sec_asteroid,2);
plot(time_consid(1:end-1),MagMAX(1:end-1),'LineWidth',1.2)
plot(time_consid(1:end-1),Magmin(1:end-1),'--','LineWidth',1.2)
xlim([files(ii).time(1) files(ii).time(end)])

figure(2)
colororder(newcolors)
time_consid=files(ii).time(files(ii).sec>sec_asteroid);
RealMagMAX=files(ii).REAL_mag(files(ii).sec>sec_asteroid,1);
RealMagmin=files(ii).REAL_mag(files(ii).sec>sec_asteroid,2);
plot(time_consid(1:end-1),RealMagMAX(1:end-1),'LineWidth',1.2)
plot(time_consid(1:end-1),RealMagmin(1:end-1),'--','LineWidth',1.2)

figure(3)
plot3(files(ii).lon,files(ii).lat,files(ii).alt,'LineWidth',1.2); 

figure(4)
skyplot(files(ii).azimuth,files(ii).elevation,'b.')

%% METEOR

R_specif=287;
r_planet = 6371000;
g0_planet = 9.81;
dt = 0.025;
t_max = 300;

ast.s=0;

%% DATA
ii=1;

ast.init_v = mean(files(ii).vel(round(files(1).alt)==initial_altitude));%files(ii).vel(2);  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                                        ast.v(1) = ast.init_v;
ast.init_zenith = 90+mean(files(ii).flight_ang(round(files(ii).alt)==initial_altitude));    %zenith angle

ast.inc_orbit=(files(ii).head_ang(1));%111 %files(ii).init_inc;%90 %-mean(files(ii).head_ang(round(files(ii).alt)==initial_altitude-100))
ast.lat(1)=files(ii).lat(1);
ast.lon(1)=files(ii).lon(1);
                                        ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
                                        ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));
ast.init_pos = r_planet + initial_altitude*1000;
                                        ast.h(1) = ast.init_pos - r_planet;

ast.r(1)=files(2).init_size/100;
ast.init_mass = 4/3 * pi * ast.r(1)^3 * ast.rho;%4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
                                        ast.m(1)=ast.init_mass;
ast.lum_eff = [0.1 0.001];%0.001;%0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006 lum_eff = [0.1 0.001];
                                        ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki

                                        ast.luminous_energy(1:2)= [nan nan];
                                        ast.luminous_energy_real(1:2)= [nan nan];

ast.z_init= z_init;  %OBS_height(ii);%up                                                                                
ast.x_init= x_init;
ast.y_init= y_init;

ast.Up_pos(1)=(ast.h(1) - ast.z_init);
ast.Nort_pos(1)=(0) + ast.x_init; %(0) + ast.x_init;%
ast.West_pos(1)=(0) + ast.y_init; %(0) + ast.y_init;%

om_earth=2*pi/(23.934472*3600);
%% CODE

mvmax = 9.0;  %# faintest magnitude to be printed
tic
t=1;
[rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
jj=1; jjt=0; ast.mv(1) = nan; ast.luminous_energy(1)=0;rho_aria(1)=rho_air;
h_dist(1)=sqrt((0 - ast.x_init)^2 + (ast.h(1) - ast.z_init)^2 + ast.y_init^2);  % distance respec to the observer

ii=2;
%     # ----------------------------------------------------------------------------------------------
%     # do the computation
%     # ----------------------------------------------------------------------------------------------
while ast.h(t) > 0 && ast.m(t) > 0.001 %ast.m(1)*0.0001 %0.0000001
    [rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
    rho_aria(t+1)=rho_air;

    mean_free_path(t)=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
    Kn(t)=mean_free_path(t)/(ast.r(t)*2);
    Cd(t)=Drag_coef(Kn(t),ast.A(t))/2;
    ast.drag_coeff=Cd(t);
    
    ast.heat_transfer=mean(files(ii).heat_transf_coef(find(round(files(ii).alt)==round(ast.h(t)/1000))));
    %ast.heat_transfer=mean(files(ii).heat_transf_coef(files(ii).heat_transf_coef>0));
    if isnan(ast.heat_transfer)
        ast.heat_transfer=0;
    end
    if isempty(ast.heat_transfer)
        ast.heat_transfer=0;
    end
    heatTransfer(t)=ast.heat_transfer;

    %om_earth=0.00007292;
    ast.lat(t+1)=rad2deg( (ast.s(t)*sin(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) + ast.lat(1);
    ast.lon(t+1)=rad2deg( (ast.s(t)*cos(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) + ast.lon(1);

    % calculate rates of change
    a2 = ast.drag_coeff * ast.A(t) * rho_air * ast.v(t)^2 / (ast.m(t) * ast.rho^2)^0.33333; %ast.drag_coeff * 1/4 * 4*pi*ast.r^2 * rho_air * ast.v(t)^2 
    gv = g0_planet / (1 + ast.h(t) / r_planet)^2;

%     aW = -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.inc_orbit)) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) * cos(deg2rad(ast.inc_orbit)) + om_earth*sin(deg2rad(ast.lat(t)))*ast.v_horizontal(t) * sin(deg2rad(ast.inc_orbit));
%     aN = -a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.inc_orbit)) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) * sin(deg2rad(ast.inc_orbit)) + om_earth*sin(deg2rad(ast.lat(t)))*ast.v_horizontal(t) * cos(deg2rad(ast.inc_orbit));
    ah = -a2 * ast.v_horizontal(t) / ast.v(t) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) + om_earth*sin(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t);
    av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit));% Eötvös effect

    ast.zd(t+1) = rad2deg(atan2(ast.v_vertical(t), ast.v_horizontal(t)));

    ast.s(t+1) = ast.s(t) + ast.v_horizontal(t) * dt;
    ast.h(t+1) = ast.h(t) + ast.v_vertical(t) * dt;

    ast.Up_pos(t+1)=(ast.h(t+1) - ast.z_init);
    ast.Nort_pos(t+1)=  (ast.s(t+1)*cos(deg2rad(-ast.inc_orbit)) + ast.x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
    ast.West_pos(t+1)=  (ast.s(t+1)*sin(deg2rad(ast.inc_orbit)) + ast.y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
    h_dist(t+1) = sqrt((ast.West_pos(t+1))^2 + (ast.Nort_pos(t+1))^2 + (ast.Up_pos(t+1))^2 );  % distance respec to the observer

    ast.v_vertical(t+1) = ast.v_vertical(t) + av * dt;
%     ast.v_horizontal(t+1) = sqrt((ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit)) + aW * dt)^2 + (ast.v_horizontal(t)*sin(deg2rad(ast.inc_orbit)) + aN * dt)^2);
    ast.v_horizontal(t+1) = ast.v_horizontal(t) + ah * dt;
    ast.v(t+1) = sqrt((ast.v_horizontal(t))^2 + ast.v_vertical(t)^2); %-ast.Earth_atm_speed*ast.v(t)/ast.v(1)

%% temperature and ablation fragmentatation

    ast.m_sec(t+1)=ast.A(t)*ast.heat_transfer/(2*ast.heat_of_ablation)*(ast.m(t)/ast.rho)^(2/3)*rho_air*ast.v(t)^3;%*abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+1);%abs(atan(0.04*(ast.T(t)-ast.melting_T))/pi+0.05);
    ast.m(t+1)=ast.m(t)-(ast.m_sec(t+1)) * dt;
    ast.r(t+1)=((ast.m(t)/ast.rho)*3/(4*pi))^(1/3);

    ast.A(t+1)=ast.A(t); 

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
% sec = seconds(1*dt:dt:t*dt);
sec = seconds([0 1*dt:dt:(t-1)*dt]);
sec.Format = 'mm:ss'; 

ast.time=sec;

r4=[ast.West_pos' ast.Nort_pos' ast.Up_pos'];
ast.azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
ast.elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

%% METEOR

figure(1)
title(titles)
subplot(2,1,1)
plot(files(ii).time,files(ii).vel/1000,'c','LineWidth',1.2)
xlabel('Time mm:ss')
ylabel('Velocity [km/s]')
legend({'Rocket','Meteor','Theory'},'Location','southwest');

subplot(2,1,2)
colororder(newcolors)
plot(ast.time(1:end-1),ast.mv(1:end-1,1),'LineWidth',1.2)
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Abs.Magnitude')
plot(ast.time(1:end-1),ast.mv(1:end-1,2),'--','LineWidth',1.2)
legend({'Rocket \tau = 0.1','Rocket \tau = 0.001','Meteor \tau = 0.1','Meteor \tau = 0.001','Theory \tau = 0.1','Theory \tau = 0.001'},'Location','southwest');


figure(2)
title(titles)
colororder(newcolors)
plot(ast.time(1:end-1),ast.mv_real(1:end-1,1),'LineWidth',1.2)
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Magnitude')
plot(ast.time(1:end-1),ast.mv_real(1:end-1,2),'--','LineWidth',1.2)
legend({'Rocket \tau = 0.1','Rocket \tau = 0.001','Meteor \tau = 0.1','Meteor \tau = 0.001','Theory \tau = 0.1','Theory \tau = 0.001'},'Location','southwest');


figure(3)
if sin(deg2rad(ast.inc_orbit))>0
    ast.lon(ast.lat>90) = 180+ast.lon(ast.lat>90);
    ast.lon(ast.lat<-90) = 180+ast.lon(ast.lat<-90);
else
    ast.lon(ast.lat>90) = ast.lon(ast.lat>90)-180;
    ast.lon(ast.lat<-90) = ast.lon(ast.lat<-90)-180;
end
ast.lon(ast.lon>180) = -360+ast.lon(ast.lon>180);
ast.lon(ast.lon<-180) = 360+ast.lon(ast.lon<-180);
ast.lat(ast.lat>90) = 180-ast.lat(ast.lat>90);
ast.lat(ast.lat<-90) = -(180+ast.lat(ast.lat<-90)); 
plot3(ast.lon,ast.lat,ast.h/1000,'LineWidth',1.2); 
title(titles)
hold on
grid on
xlabel('Longitude [degrees]','Fontsize',10);
ylabel('Latitude [degrees]','Fontsize',10);
zlabel('Altitude [km]','Fontsize',10);
% view([0 0 90])


figure(4)
skyplot(ast.azimuth,ast.elevation,'c.')
title(titles)

%% EARTH
figure(3)
plot3(lon_OBS,lat_OBS,z_init/1000,'*r','LineWidth',1)   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
lat = load('coast_lat.dat');
long = load('coast_long.dat');
grid on
plot3(long,lat,zeros(1,length(lat)),'k')
xlim([-40 40])
ylim([-40 40])
% xlim([-180 180])
% ylim([-90 90])
% zlim([0 500])
set(gca,'Xtick',-180:30:180)
set(gca,'Ytick',-90:30:90)
% view([0 0 90])
hleg = legend({'Rocket','Meteor','Theory'},'Location','southwest');

allChildren = get(hleg, 'Children');                % list of all objects on axes
displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
% Remove object associated with "data1" in legend
delete(allChildren(strcmp(displayNames, 'data1')))

% figure
% plot3(files(1).lon,files(1).lat,files(1).alt,'.g','LineWidth',1.2); 
% hold on
% plot3(files(1).lon(files(1).alt<initial_altitude),files(1).lat(files(1).alt<initial_altitude),files(1).alt(files(1).alt<initial_altitude),'.r','LineWidth',1.2); 
% plot3(files(2).lon(files(2).alt<initial_altitude),files(2).lat(files(2).alt<initial_altitude),files(2).alt(files(2).alt<initial_altitude),'.r','LineWidth',1.2); 
% plot3(long,lat,zeros(1,length(lat)),'k')
% set(gca,'Xtick',-180:30:180)
% set(gca,'Ytick',-90:30:90)
% grid on
