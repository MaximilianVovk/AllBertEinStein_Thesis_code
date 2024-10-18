clc
clear all
close all
%% DATA


vel=[7.5 9 11.2 15 30 72];

radius=1;

head_ang = [0 45 90 135 180];

zenith_angs = [45 65 79.5];

for incl_num=1:length(zenith_angs)
for ivel=1:length(vel)
    for iang=1:length(head_ang)
        run('config.m')
ast.init_v = vel(ivel)*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                                        ast.v(1) = ast.init_v;
ast.init_zenith = zenith_angs(incl_num);    %zenith angle

ast.inc_orbit=head_ang(iang);%files(ii).init_inc;%90
ast.lat(1)=0;
ast.lon(1)=0;
                                        ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
                                        ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));
ast.init_pos = r_planet + 130*1000;
                                        ast.h(1) = ast.init_pos - r_planet;

ast.r(1) = radius; %in meters
ast.rho = 2800.0;   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
ast.A = 1.2;    % shape factor (1.2 for sphere). %A=S/V^2/3 (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);%
ast.drag_coeff = 0.5;%0.5; % due to aircap actual Cd=2drag_coef
ast.heat_transfer = 0.00002; %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
ast.heat_of_ablation = 500000;   %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6

ast.init_mass = 4/3 * pi * ast.r^3 * ast.rho;%4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
                                        ast.m(1)=ast.init_mass;
ast.lum_eff = [0.1 0.001];%0.001;%0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006 lum_eff = [0.1 0.001];
                                        ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki

                                        ast.luminous_energy(1:2)= [nan nan];
                                        ast.luminous_energy_real(1:2)= [nan nan];

ast.z_init= 10000;  %OBS_height(ii);%up                                        
ast.x_init= -9.3601e+04; %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
ast.y_init= -1.0177e+06; %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north

ast.Up_pos(1)=(ast.h(1) - ast.z_init);
ast.Nort_pos(1)=(0) + ast.x_init; %(0) + ast.x_init;%
ast.West_pos(1)=(0) + ast.y_init; %(0) + ast.y_init;%

%% CODE

mvmax = 9.0;  %# faintest magnitude to be printed
tic
t=1; jj=1; jjt=0; ast.mv(1) = nan; ast.luminous_energy(1)=0;rho_aria(1)=rho_air;
h_dist(1)=sqrt((0 - obs.x_init)^2 + (ast.h(1) - obs.z_init)^2 + obs.y_init^2);  % distance respec to the observer

%     # ----------------------------------------------------------------------------------------------
%     # do the computation
%     # ----------------------------------------------------------------------------------------------
        

        dragH{ivel,iang}(t)= 0;
        EotvosH{ivel,iang}(t)= 0;
        CoriolisH{ivel,iang}(t)= 0;

        dragV{ivel,iang}(t)=0;
        EotvosV{ivel,iang}(t)= 0;
        CentrifV{ivel,iang}(t)= 0;
        GraV{ivel,iang}(t)= 0;

        dragSide{ivel,iang}(t)= 0;
        EotvosS{ivel,iang}(t)=0;
        
ast.v_side(t)=0; ast.side(t)=0; init_inc=ast.inc_orbit; ast.AnghSide(t)=0;
while ast.h(t) > 0 && ast.m(t) > 0.001
    [rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
    rho_aria(t+1)=rho_air;

    mean_free_path=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
    Kn(t)=mean_free_path/(ast.r(t)*2);
    Cd(t)=Drag_coef(Kn(t))/2;
    ast.drag_coeff=Cd(t);

    %om_earth=0.00007292;
    om_earth=2*pi/(23.934472*3600);
    ast.lat(t+1)=rad2deg( (ast.s(t)*sin(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lat(1);
    ast.lon(t+1)=rad2deg( (ast.s(t)*cos(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lon(1);
    
    % calculate rates of change
    a2 = ast.drag_coeff * ast.A(t) * rho_air * ast.v(t)^2 / (ast.m(t) * ast.rho^2)^0.33333; %ast.drag_coeff * 1/4 * 4*pi*ast.r^2 * rho_air * ast.v(t)^2 
    gv = g0_planet / (1 + ast.h(t) / r_planet)^2;

    dragH{ivel,iang}(t+1)= -a2 * ast.v_horizontal(t) / ast.v(t)*cos(deg2rad(ast.AnghSide(t)));
    EotvosH{ivel,iang}(t+1)= 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
%     EotvosH{ivel,iang}(t+1)= abs(2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit)))  +   abs(2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t))));
    CoriolisH{ivel,iang}(t+1)= - ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t));
%     ah = dragH{ivel,iang}(t+1) + CoriolisH{ivel,iang}(t+1) + EotvosH{ivel,iang}(t+1);
%     ah = -a2 * ast.v_horizontal(t) / ast.v(t) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) + om_earth*sin(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t);

    dragV{ivel,iang}(t+1)=- a2 * ast.v_vertical(t) / ast.v(t);
    EotvosV{ivel,iang}(t+1)=2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))-2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
%     EotvosV{ivel,iang}(t+1)=  abs(2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))) + abs(2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit)));% Eötvös effect
    CentrifV{ivel,iang}(t+1)= ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t));
    GraV{ivel,iang}(t+1)=-gv; %Eötvös effect

    dragSide{ivel,iang}(t+1)=-a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.AnghSide(t)));
    EotvosS{ivel,iang}(t+1)=2*om_earth * ast.v_horizontal(t)* sin(deg2rad(ast.lat(t))) - 2*om_earth * ast.v_vertical(t) * cos(deg2rad(ast.lat(t))) * sin(deg2rad(ast.inc_orbit));
%     EotvosS{ivel,iang}(t+1)=abs(2*om_earth * ast.v_horizontal(t)* sin(deg2rad(ast.lat(t)))) + abs(2*om_earth * ast.v_vertical(t) * cos(deg2rad(ast.lat(t))) * sin(deg2rad(ast.inc_orbit))); %Eötvös effect
%     av = GraV{ivel,iang}(t+1) + dragV{ivel,iang}(t+1) + CentrifV{ivel,iang}(t+1) + EotvosV{ivel,iang}(t+1);
%     av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit));% Eötvös effect

    ah = -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
    av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
    aSide = -a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.AnghSide(t))) + 2*om_earth * ast.v_horizontal(t)* sin(deg2rad(ast.lat(t))) - 2*om_earth * ast.v_vertical(t) * cos(deg2rad(ast.lat(t))) * sin(deg2rad(ast.inc_orbit));
%     ah = -a2 * ast.v_horizontal(t) / ast.v(t) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) + om_earth*sin(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t);
%     av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit));% Eötvös effect


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
sec = seconds(1*dt:dt:t*dt);
sec.Format = 'mm:ss'; 

sec_CASES{ivel,iang}=sec;
alt_CASES{ivel,iang}=ast.h;
space_CASES{ivel,iang}=ast.s;

lat_CASES{ivel,iang}=ast.lat;
lon_CASES{ivel,iang}=ast.lon;

r4=[ast.West_pos' ast.Nort_pos' ast.Up_pos'];
ast.azimuth=rad2deg(atan2(r4(:,2),r4(:,1)));
ast.elevation=rad2deg(atan2(r4(:,3),sqrt(r4(:,1).^2 + r4(:,2).^2)));

clear ast t
end
end

%% PLOT Bar

y = [1 2 3 4 5 6]; % cases
z = [nan nan nan nan nan nan nan; %material
     nan nan nan nan nan nan nan; 
     nan nan nan nan nan nan nan;
     nan nan nan nan nan nan nan;
     nan nan nan nan nan nan nan;
     nan nan nan nan nan nan nan];
 
% for iang=1:length(head_ang)-2
% for iang=2:2
%     for ivel=1:length(vel)
%         for cases=1:7
%            switch cases
%                 case 1
%                     zacc{iang}(ivel,1)=sum(abs(EotvosH{ivel,iang}));
%                 case 2
%                     zacc{iang}(ivel,2)=sum(abs(EotvosV{ivel,iang}));
%                 case 3
%                     zacc{iang}(ivel,3)=sum(abs(GraV{ivel,iang}));
%                 case 4
%                     zacc{iang}(ivel,4)=sum(abs(CentrifV{ivel,iang}));
%                 case 5
%                     zacc{iang}(ivel,5)=sum(abs(CoriolisH{ivel,iang}));
%                 case 6
%                     zacc{iang}(ivel,6)=sum(abs(dragV{ivel,iang}));
%                 case 7
%                     zacc{iang}(ivel,7)=sum(abs(dragH{ivel,iang}));
%            end
%            
%         end
% %         totAccsum=sum(accTOT{iang}(ivel,:));
% %         accTOT{iang}(ivel,:)=100*accTOT{iang}(ivel,:)./totAccsum;
%     end
%     figure(iang)
%     bar3(y,zacc{iang})
%     zlabel('Tot.Acc. [m/s^2]')
% %     bar3(accTOT{iang},0.5,'stacked')
% %     zlabel('Tot.Acc. percent [%]')
%     xticklabels({'Eötvös Horizontal','Eötvös Vertical','Gravity Acc.','Centrifugal Vertical','Coriolis Horizontal','Drag Vertical','Drag Horizontal'})
% %     legend({'Eötvös Horizontal','Eötvös Vertical','Gravity Acc.','Centrifugal Vertical','Coriolis Horizontal','Drag Vertical','Drag Horizontal'})
% %     zlim([0 100])
%     yticklabels({'7.5 km/s','9 km/s','11.2 km/s','15 km/s','30 km/s','72 km/s'})
%     %vel=[7.5 9 11.2 15 30 72];
%     disp('Total Acceleration Terms of 1 m Meteor')
%     title(['Heading Angle ',num2str(head_ang(iang)),'°'])
% end
%  
% for iang=1:length(head_ang)-2
%     for ivel=1:length(vel)
%         for cases=1:5
%            switch cases
%                 case 1
%                     accTOT{iang}(ivel,1)=sum(abs(EotvosH{ivel,iang}));
%                 case 2
%                     accTOT{iang}(ivel,2)=sum(abs(EotvosV{ivel,iang}));
%                 case 3
%                     accTOT{iang}(ivel,3)=sum(abs(GraV{ivel,iang}));
%                 case 4
%                     accTOT{iang}(ivel,4)=sum(abs(CentrifV{ivel,iang}));
%                 case 5
%                     accTOT{iang}(ivel,5)=sum(abs(CoriolisH{ivel,iang}));
%                 case 6
%                     accTOT{iang}(ivel,6)=sum(abs(dragV{ivel,iang}));
%                 case 7
%                     accTOT{iang}(ivel,7)=sum(abs(dragH{ivel,iang}));
%            end
%            
%         end
%         totAccsum=sum(accTOT{iang}(ivel,:));
%         accTOT{iang}(ivel,:)=100*accTOT{iang}(ivel,:)./totAccsum;
%     end
%     figure(iang+3)
% %     bar3(y,accTOT{iang})
% %     zlabel('Tot.Acc. [m/s^2]')
%     bar3(accTOT{iang},0.5,'stacked')
%     zlabel('Tot.Acc. percent [%]')
% %     xticklabels({'Eötvös Horizontal','Eötvös Vertical','Gravity Acc.','Centrifugal Vertical','Coriolis Horizontal','Drag Vertical','Drag Horizontal'})
%     legend({'Eötvös Horizontal','Eötvös Vertical','Gravity Acc.','Centrifugal Vertical','Coriolis Horizontal','Drag Vertical','Drag Horizontal'})
%     zlim([0 100])
%     yticklabels({'7.5 km/s','9 km/s','11.2 km/s','15 km/s','30 km/s','72 km/s'})
%     %vel=[7.5 9 11.2 15 30 72];
%     disp('Total Acceleration Terms of 1 m Meteor')
%     title(['Heading Angle ',num2str(head_ang(iang)),'°'])
% end
for iang=1:length(head_ang)
    for ivel=1:length(vel)
        TOTEarth(ivel,iang)=sum(abs(EotvosH{ivel,iang}))+sum(abs(EotvosS{ivel,iang}));
        TOTEotvos(ivel,iang)=sum(abs(EotvosV{ivel,iang}));
        TOTGrav(ivel,iang)=sum(abs(GraV{ivel,iang}));
        TOTCentr(ivel,iang)=sum(abs(CentrifV{ivel,iang}));
        TOTCoriol(ivel,iang)=sum(abs(CoriolisH{ivel,iang}));
        TOTDrag(ivel,iang)=sum(abs(dragV{ivel,iang}))+sum(abs(dragH{ivel,iang}))+sum(abs(dragSide{ivel,iang}));
    end
end
% for iang=2:2
    iang=2;
    for ivel=1:length(vel)
        for cases=1:6
           switch cases
                case 1
                    zacc{iang}(ivel,1)=TOTEarth(ivel,iang);
                case 2
                    zacc{iang}(ivel,2)=TOTEotvos(ivel,iang);
                case 3
                    zacc{iang}(ivel,3)=TOTGrav(ivel,iang);
                case 4
                    zacc{iang}(ivel,4)=TOTCoriol(ivel,iang);
                case 5
                    zacc{iang}(ivel,5)=TOTCentr(ivel,iang);
                case 6
                    zacc{iang}(ivel,6)=TOTDrag(ivel,iang);
           end
           
        end
    end
%     figure(iang)
figure
    bar3(y,zacc{iang})
    zlabel('Tot.Acc. [m/s^2]')
    xticklabels({'Earth Coriolis','Eötvös Effect','Gravity Acc.','Coriolis Acc.','Centrifugal Acc.','Drag Acc.'})
    yticklabels({'7.5 km/s','9 km/s','11.2 km/s','15 km/s','30 km/s','72 km/s'})
    disp('Total Acceleration Terms for a 1 m Meteor')
    title(['Head.Ang. ',num2str(head_ang(iang)),'° Zen.Ang. ',num2str(zenith_angs(incl_num)),'°'])
    zlim([0 4*10^6])
% end
 
for iang=1:length(head_ang)-2
    for ivel=1:length(vel)
        for cases=1:6
               
           switch cases
                case 1
                    accTOT{iang}(ivel,1)=TOTEarth(ivel,iang);
                case 2
                    accTOT{iang}(ivel,2)=TOTEotvos(ivel,iang);
                case 3
                    accTOT{iang}(ivel,3)=TOTGrav(ivel,iang);
                case 4
                    accTOT{iang}(ivel,4)=TOTCoriol(ivel,iang);
                case 5
                    accTOT{iang}(ivel,5)=TOTCentr(ivel,iang);
                case 6
%                     accTOT{iang}(ivel,6)=TOTDrag(ivel,iang);
                    accTOT{iang}(ivel,6)=0;
           end
               
           
        end
        totAccsum=sum(accTOT{iang}(ivel,1:4))+TOTDrag(ivel,iang);
        accTOT{iang}(ivel,:)=100*accTOT{iang}(ivel,:)./totAccsum;
    end
    figure
    bar3(accTOT{iang},0.5,'stacked')
%     bar(accTOT{iang},'stacked')
    grid on
    zlabel('Tot.Acc. [%]')
    legend({'Earth Coriolis','Eötvös Effect','Gravity Acc.','Coriolis Acc.','Centrifugal Acc.','Drag Acc.'})
%     xticklabels({'Earth Coriolis','Eötvös Effect','Gravity Acc.','Coriolis Acc.','Centrifugal Acc.','Drag Acc.'})
%     zlim([0 sum(accTOT{1}(1,:))+1])
    zlim([0 25])
    yticklabels({'7.5 km/s','9 km/s','11.2 km/s','15 km/s','30 km/s','72 km/s'})
    disp('Accelerations without showing Drag for a 1 m Meteor')
%     view([-90 0 0])
    title(['Head.Ang. ',num2str(head_ang(iang)),'° Zen.Ang. ',num2str(zenith_angs(incl_num)),'°'])
end

% colorVel=[0 0.4470 0.7410;...
%     0.8500 0.3250 0.0980;...
%     0.9290 0.6940 0.1250;...
%     0.4940 0.1840 0.5560;...
%     0.4660 0.6740 0.1880;...
%     0.3010 0.7450 0.9330;...
%     0.6350 0.0780 0.1840];%"MarkerEdgeColor",[0.6350 0.0780 0.1840]
% 
% % positz=[1 2 3 4 5 6 7];
% positz=[7 6 5 4 3 2 1];

% figure
% lat = load('coast_lat.dat');
% long = load('coast_long.dat');
% plot3(long,lat,zeros(1,length(lat)),'k')
% grid on
% hold on
% % xlim([-180 180])
% % ylim([-90 90])
% % xlim([-10 10])
% xlim([-max(lat_CASES{6,3}) max(lat_CASES{6,3})])
% ylim([min(lon_CASES{6,5}) max(lon_CASES{6,1})])
% % zlim([0 10])
% set(gca,'Xtick',-180:30:180)
% set(gca,'Ytick',-90:30:90)
% xlabel('Longitude [degrees]','Fontsize',10);
% ylabel('Latitude [degrees]','Fontsize',10);
% zlabel('Altitude [km]','Fontsize',10);
% view([0 0 90])
% for ivel=1:length(vel)
%     for iang=1:length(head_ang)
%         plot3(lon_CASES{ivel,iang},lat_CASES{ivel,iang},positz(ivel)*ones(1,length(lon_CASES{ivel,iang})),'.',"MarkerEdgeColor",colorVel(ivel,:),'LineWidth',1.2)   
% %         plot3(lon_CASES{ivel,iang},lat_CASES{ivel,iang},alt_CASES{ivel,iang},'.',"MarkerEdgeColor",colorVel(ivel,:),'LineWidth',1.2)  
%     end
% end


for iang=1:length(head_ang)
    for ivel=1:length(vel)
        spaceTOT(ivel,iang)=space_CASES{ivel,iang}(end)/1000;
    end
end
figure
% bar3(spaceTOT,1,'grouped')
bar(spaceTOT,'grouped')
% zlabel('Downrange [km]')
ylabel('Downrange [km]')
grid on
hleg = legend({'0°','45°','90°','135°','180°'},'location','northeastoutside');
htitle = get(hleg,'Title');
set(htitle,'String','Orbit Inclinations')
% xticklabels({'0°','45°','90°','135°','180°'})
% yticklabels({'7.5 km/s','9 km/s','11.2 km/s','15 km/s','30 km/s','72 km/s'})
xticklabels({'7.5 km/s','9 km/s','11.2 km/s','15 km/s','30 km/s','72 km/s'})
disp('Total Downrange for a 1 m Meteor')
title(['Zenith Angle ',num2str(zenith_angs(incl_num)),'°'])

clearvars -except vel radius head_ang zenith_angs
end

%% SAVE .jpg

% FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\Oct';   % Your destination folder
FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\3 Chap'; 
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    figure(iFig)
%     set(gca,'FontSize',15)
    saveas(gcf,fullfile(FolderName, [num2str(iFig),'_Theory_diffInclFlightAng.jpg']))
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   savefig(fullfile(FolderName, [FigName '.fig']));
end