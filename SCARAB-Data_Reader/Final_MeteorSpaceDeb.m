clc
clear all
close all
%% COST
run('config.m')

namePlot='Space Debris';

R_specif=287;
r_planet = 6371000;
g0_planet = 9.81;
dt = 0.025;
t_max = 300;

ast.s=0;

%% DATA

ast.init_v = 7.970283806*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                                        ast.v(1) = ast.init_v;
ast.init_zenith = 90-1.836319289;%90-1.836319289;    %zenith angle

ast.inc_orbit=90;%files(ii).init_inc;%90
ast.lat(1)=0;
ast.lon(1)=0;
                                        ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
                                        ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));
ast.init_pos = r_planet + 130.5925573*1000;
                                        ast.h(1) = ast.init_pos - r_planet;

ast.r(1) = 0.1; %in meters
ast.rho = 2800.0;   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
ast.A = 1.2;    % shape factor (1.2 for sphere). %A=S/V^2/3 (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);%
ast.drag_coeff = 0.5;%0.5; % due to aircap actual Cd=2drag_coef
ast.heat_transfer = 0.02; %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
ast.heat_of_ablation = 500000;   %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6

ast.init_mass = 4/3 * pi * ast.r(1)^3 * ast.rho;%4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
                                        ast.m(1)=ast.init_mass;
ast.lum_eff = [0.1 0.001];%0.001;%0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006 lum_eff = [0.1 0.001];
                                        ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki

                                        ast.luminous_energy(1:2)= [nan nan];
                                        ast.luminous_energy_real(1:2)= [nan nan];

ast.z_init= 10000;  %OBS_height(ii);%up                                                                                
lat_OBS=2;
lon_OBS=0;
ast.x_init= (ast.z_init+r_planet)*deg2rad(ast.lon(1)-lon_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
ast.y_init= (ast.z_init+r_planet)*deg2rad(ast.lat(1)-lat_OBS); %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north

% ast.x_init= -9.3601e+04; %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
% ast.y_init= -1.0177e+06; %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north
% lon_OBS=ast.lon(1)-rad2deg(ast.x_init/(ast.z_init+r_planet));
% lat_OBS=ast.lat(1)-rad2deg(ast.y_init/(ast.z_init+r_planet));

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

mean_free_path(t)=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
Kn(t)=mean_free_path(t)/(ast.r(t)*2);
Cd(t)=Drag_coef(Kn(t),ast.A(t))/2;

ast.v_side(t)=0; ast.side(t)=0; init_inc=ast.inc_orbit; ast.AnghSide(t)=0;
%     # ----------------------------------------------------------------------------------------------
%     # do the computation
%     # ----------------------------------------------------------------------------------------------
while ast.h(t) > 0 && ast.m(t) > 0.001
%     rho_air=mean(files(ii).rho_air(find(round(files(ii).alt)==round(ast.h(t)/1000))));
%     rho_aria(t+1)=rho_air;
%     if isnan(rho_air)
        [rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
        rho_aria(t+1)=rho_air;
%     end
    
    mean_free_path(t)=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
%     Kn(t+1)=mean(files(ii).Kn(find(round(files(ii).alt)==round(ast.h(t)/1000))));
%     if isnan(Kn(t+1))
        Kn(t+1)=mean_free_path(t)/(ast.r(t)*2);
%     end
    Cd(t+1)=Drag_coef(Kn(t+1))/2;
    ast.drag_coeff=Cd(t+1);
    
    
%       ast.heat_transfer=mean(files(ii).heat_transf_coef(find(round(files(ii).alt)==round(ast.h(t)/1000))));
% %     ast.heat_transfer=mean(files(ii).coef_gauss(find(round(1:100)==round(ast.h(t)/1000))));
    if isnan(ast.heat_transfer)
        ast.heat_transfer=0;
    end
    if isempty(ast.heat_transfer)
        ast.heat_transfer=0;
    end
    heatTransfer(t)=ast.heat_transfer;

    %om_earth=0.00007292;
    om_earth=2*pi/(23.934472*3600);
    ast.lat(t+1)=rad2deg( (ast.s(t)*sin(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lat(1);
    ast.lon(t+1)=rad2deg( (ast.s(t)*cos(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lon(1);

    % calculate rates of change
    a2 = ast.drag_coeff * ast.A(t) * rho_air * ast.v(t)^2 / (ast.m(t) * ast.rho^2)^0.33333; %ast.drag_coeff * 1/4 * 4*pi*ast.r^2 * rho_air * ast.v(t)^2 
    gv = g0_planet / (1 + ast.h(t) / r_planet)^2;

    ah = -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
    av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
    aSide = -a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.AnghSide(t))) + 2*om_earth * ast.v_horizontal(t)* sin(deg2rad(ast.lat(t))) - 2*om_earth * ast.v_vertical(t) * cos(deg2rad(ast.lat(t))) * sin(deg2rad(ast.inc_orbit));
%     ah = -a2 * ast.v_horizontal(t) / ast.v(t) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) + om_earth*sin(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t);
%     av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit));% Eötvös effect


    ast.zd(t+1) = rad2deg(atan2(ast.v_vertical(t), ast.v_horizontal(t)));
    ast.AnghSide(t+1) = rad2deg(atan2(ast.side(t),ast.s(t)));
    ast.inc_orbit=init_inc- ast.AnghSide(t+1);
    control_inc(t)=ast.inc_orbit;

    ast.s(t+1) = ast.s(t) + ast.v_horizontal(t) * dt * r_planet / (r_planet + ast.h(t));
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

r4=[ast.West_pos' ast.Nort_pos' ast.Up_pos'];
ast.azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
ast.elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

%% ALTITUDE DOWNRANGE
figure
plot(ast.s/1000,(ast.h)/1000,'r')
title(namePlot)
grid on
hold on
ylabel('Altitude [km]')
xlabel('Dowrange [km]')

%% VEL PLOT
figure
plot(sec,ast.v/1000,'r')
title(namePlot)
grid on
hold on
ylabel('Velocity [km/s]')
xlabel('Time mm:ss')

%% MASS PLOT
figure
subplot(2,1,1)
plot(sec,ast.m,'r')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Mass [Kg]')

subplot(2,1,2)
plot(sec,ast.m_sec,'r')
sgtitle(namePlot)
grid on
hold on
xlabel('Time mm:ss')
ylabel('Mass loss rate [Kg/s]')

%% MAGNITUDE PLOT
figure
subplot(2,1,1)
plot(sec(2:end)',ast.mv(2:end,1),'r')
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Abs.Magnitude')
plot(sec(2:end)',ast.mv(2:end,2),'r--')
legend({'MAX Theory','min Theory'},'Location','southwest')

subplot(2,1,2)
plot(sec(2:end)',ast.mv_real(2:end,1),'r')
set(gca, 'YDir','reverse')
grid on
hold on
xlabel('Time mm:ss')
ylabel('Magnitude')
plot(sec(2:end)',ast.mv_real(2:end,2),'r--')
legend({'MAX Theory','min Theory'},'Location','southwest')

sgtitle(namePlot)

figure
subplot(2,1,1)
plot(ast.mv(:,1),ast.h/1000,'r');
set(gca, 'XDir','reverse')
grid on
hold on
ylim([0 100])
xlim([-10 10])
plot(ast.mv(:,2),ast.h/1000,'r--');
legend({'MAX Theory','min Theory'},'Location','southwest')
xlabel('Abs.Magnitude') 
ylabel('Altitude [km]')

subplot(2,1,2)
plot(ast.mv_real(:,1),ast.h/1000,'r');
set(gca, 'XDir','reverse')
grid on
hold on
ylim([0 100])
xlim([-10 10])
plot(ast.mv_real(:,2),ast.h/1000,'r--');
legend({'MAX Theory','min Theory'},'Location','southwest')
xlabel('Magnitude') 
ylabel('Altitude [km]')

sgtitle(namePlot)

%% EARTH PLOT
figure

% if sin(deg2rad(ast.inc_orbit))>0
%     ast.lon(ast.lat>90) = 180+ast.lon(ast.lat>90);
%     ast.lon(ast.lat<-90) = 180+ast.lon(ast.lat<-90);
% else
%     ast.lon(ast.lat>90) = ast.lon(ast.lat>90)-180;
%     ast.lon(ast.lat<-90) = ast.lon(ast.lat<-90)-180;
% end
% ast.lon(ast.lon>180) = -360+ast.lon(ast.lon>180);
% ast.lon(ast.lon<-180) = 360+ast.lon(ast.lon<-180);
% ast.lat(ast.lat>90) = 180-ast.lat(ast.lat>90);
% ast.lat(ast.lat<-90) = -(180+ast.lat(ast.lat<-90));        
plot3(ast.lon,ast.lat,ast.h/1000,'.r','LineWidth',1.2)   
hold on
plot3(lon_OBS,lat_OBS,ast.z_init/1000,'*b','LineWidth',1)   %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
grid on
lat = load('coast_lat.dat');
long = load('coast_long.dat');
plot3(long,lat,zeros(1,length(lat)),'k')
xlim([-180 180])
ylim([-90 90])
zlim([0 500])
set(gca,'Xtick',-180:30:180)
set(gca,'Ytick',-90:30:90)
xlabel('Longitude [degrees]','Fontsize',10);
ylabel('Latitude [degrees]','Fontsize',10);
zlabel('Altitude [km]','Fontsize',10);
view([0 0 90])
title(namePlot)

%% Azimuth Elevation
figure
skyplot(ast.azimuth(round(ast.azimuth)~=round(ast.azimuth(end))),ast.elevation(round(ast.azimuth)~=round(ast.azimuth(end))),'r.')
hold on
%skyplot(ast.azimuth(ast.elevation>5),ast.elevation(ast.elevation>5),'g.')
hold on
title(namePlot)

%% Knusseld number and Drag
figure
subplot(2,1,1)
semilogy(sec,1*ones(1,length(Kn)),'k','LineWidth',1.2);
grid on
hold on
semilogy(sec,0.01*ones(1,length(Kn)),'k','LineWidth',1.2);
semilogy(sec,Kn,'r');

subplot(2,1,2)
plot(sec,Cd*2,'r');
grid on
ylabel('C_D Drag Coef.')
xlabel('Time mm:ss')

figure
semilogx(Kn,(ast.h)/1000,'w');
hold on
semilogx(1*ones(1,length(Kn)),(ast.h)/1000,'k','LineWidth',1.2);
grid on
semilogx(0.01*ones(1,length(Kn)),(ast.h)/1000,'k','LineWidth',1.2);

xlabel('Kn Number')
ylabel('Altitude [km]')
ylim([0 130])
xlim([10^-8 1000])

%% Pressure
ast.stagn_pres=rho_aria.*ast.v.^2;

CFD_stagP=[9.9e4 3.1e4 5.5e3];
CFD_HuskinP=[9.9e4/2 3.1e4/2 5.5e3/2];
% CFD_HuskinP=[1e4 2e3 5e2];
CFD_alt=[40 55 70];

figure
semilogx(ast.stagn_pres,ast.h/1000,'b')%,sec(jjt),ast.dynamic_pres(jjt))
title('Fragmentation')
ylabel('Altitude [km]') 
xlabel('Pressure [Pa]')
grid on
hold on
ylim([0 120])
xlim([min(ast.stagn_pres) ast.ultim_compres*10])
semilogx(CFD_stagP,CFD_alt,'b*')
%plot(sec,ones(1,length(sec(jjt)))*ast.ultim_compres,'b',sec(jjt),ones(1,length(sec(jjt)))*ast.ultim_tensile,'r')
plot(ones(1,length(sec))*ast.ultim_compres,ast.h/1000,'r','linewidth',1.2)
legend('stagnation presure','CFD stagnation presure','ultimate compressive strength')

figure
semilogx(ast.stagn_pres*0.5,ast.h/1000,'b')%,sec(jjt),ast.dynamic_pres(jjt)) %ast.stagn_pres*0.365
title('Husking')
ylabel('Altitude [km]') 
xlabel('Pressure [Pa]')
grid on
hold on
ylim([0 120])
xlim([min(ast.stagn_pres) ast.ultim_compres*10])
semilogx(CFD_HuskinP,CFD_alt,'b*')
%plot(sec,ones(1,length(sec(jjt)))*ast.ultim_compres,'b',sec(jjt),ones(1,length(sec(jjt)))*ast.ultim_tensile,'r')
plot(ones(1,length(sec))*ast.ultim_tensile,ast.h/1000,'r','linewidth',1.2)
legend('marginal zone presure','CFD marginal zone presure','ultimate shear strength')

%% Spraing

figure
%plot(sec,7.6*ast.surface_tension./(drag_coef.*rho_aria.*ast.v.^2))
plot(7.6*ast.surface_tension./(Cd.*rho_aria.*ast.v.^2),ast.h/1000,'b')
title('Radius for fusion')
hold on
plot(ones(1,length(sec))*ast.r(1),ast.h/1000,'r','linewidth',1.2)
ylabel('Altitude [km]')  
xlabel('Melting Radius [m]')
xlim([0 0.5])
grid on
ylim([0 120])


%% Vapour
heat_stag=1.7415e-4.*(rho_aria./ast.r).^0.5.*ast.v.^3; %J/s/m^2 experimental sutton graves formula
ast.T=(heat_stag/(ast.emissivity.*Stefan_Boltzmann_cost)).^(1/4);
if ast.rho < 3800 %ROCKY asteroid
    ast.vap_pres=10.^(10.6-13500./ast.T).*0.1;
else %IRON asteroid
    ast.vap_pres=10.^(10.607-16120./ast.T).*0.1;
end

figure; 
plot(ast.stagn_pres.*Cd,ast.h/1000,'r','linewidth',1.2);
hold on
plot(ast.vap_pres,ast.h/1000,'b'); 
grid on
legend('Stagnation pressure','Vapour pressure')
ylabel('Altitude [km]')  
xlabel('Pressure [Pa]')
% title('ps>pe for evaporation')
title('Stony Meteor')
ylim([0 120])

figure; 
plot(ast.T,ast.h/1000,'b');
hold on
% plot(ast.vap_pres,ast.h/1000,'r','linewidth',1.2); 
grid on
% legend('Stagnation pressure','Vapour pressure')
ylabel('Altitude [km]')  
xlabel('Temperature [K]')
% title('ps>pe for evaporation')
title('Stony Meteor')
ylim([0 120])

%% MASS
disp(append('Initial mass : ',num2str(ast.m(1)),' kg'))
disp(append('Theoretical final mass : ',num2str(ast.m(end)),' kg'))



% FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\Oct';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    figure(iFig)
%     saveas(gcf,fullfile(FolderName, [num2str(iFig),'_Theory_diffInclFlightAng.jpg']))
% set(gca,'FontSize',20)
end
