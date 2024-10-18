clc
clear all
close all
%% DATA

% % load('meteorMass.mat')
% load('SACARBMissionKnMa.mat')
% for ii=1:length(files)
%     if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))   
%         run('config.m')
% 
%         ast.init_v = files(ii).init_vel*1000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
%                                         ast.v(1) = ast.init_v;
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
%         ast.init_lat=0;
%         ast.init_lon=0;
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
%         %% CODE
% 
%         mvmax = 9.0;  %# faintest magnitude to be printed
%         tic
%         t=1; jj=1; jjt=0; ast.mv(1) = nan; ast.luminous_energy(1)=0;rho_aria(1)=rho_air;
%         h_dist(1)=sqrt((0 - obs.x_init)^2 + (ast.h(1) - obs.z_init)^2 + obs.y_init^2);  % distance respec to the observer
% 
%         ast.lat(t)=0;
%         ast.lon(t)=0;
% ast.z_init= 10000;  %OBS_height(ii);%up                                        
% ast.x_init= -9.3601e+04; %(ast.z_init+r_planet)*deg2rad(ast.lon(1)-OBS_lon(ii));%west
% ast.y_init= -1.0177e+06; %(ast.z_init+r_planet)*deg2rad(ast.lat(1)-OBS_lat(ii));%north
% ast.v_side(t)=0; ast.side(t)=0; init_inc=ast.inc_orbit;
% 
% 
% ast.v_side(t)=0; ast.side(t)=0; init_inc=ast.inc_orbit; ast.AnghSide(t)=0;
% %     # ----------------------------------------------------------------------------------------------
% %     # do the computation
% %     # ----------------------------------------------------------------------------------------------
% while ast.h(t) > 0 && ast.m(t) > 0.001
%     rho_air=mean(files(ii).rho_air(find(round(files(ii).alt)==round(ast.h(t)/1000))));
% %     rho_air=mean(files(ii).rho_air(find(round(files(ii).sec)==round(t*dt))));
%     rho_aria(t+1)=rho_air;
%     if isnan(rho_air)
%         [rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
%         rho_aria(t+1)=rho_air;
%     end
%     
%     mean_free_path(t)=nu_air*rho_air/P_air*sqrt(pi*R_specif*T_air/2);
%     Kn(t+1)=mean(files(ii).Kn(find(round(files(ii).alt)==round(ast.h(t)/1000))));
% %     Kn(t+1)=mean(files(ii).Kn(find(round(files(ii).sec)==round(t*dt))));
%     if isnan(Kn(t+1))
%         Kn(t+1)=mean_free_path(t)/(ast.r(t)*2);
%     end
%     Cd(t+1)=Drag_coef(Kn(t+1))/2;
%     ast.drag_coeff=Cd(t+1);
%     
%     ast.heat_transfer=mean(files(ii).heat_transf_coef(find(round(files(ii).sec)==round(t*dt))));
% %        ast.heat_transfer=mean(files(ii).heat_transf_coef(find(round(files(ii).alt)==round(ast.h(t)/1000))));
% % %     ast.heat_transfer=mean(files(ii).coef_gauss(find(round(1:100)==round(ast.h(t)/1000))));
%     if isnan(ast.heat_transfer)
%         ast.heat_transfer=0;
%     end
%     if isempty(ast.heat_transfer)
%         ast.heat_transfer=0;
%     end
%     heatTransfer(t)=ast.heat_transfer;
% 
%     %om_earth=0.00007292;
%     om_earth=2*pi/(23.934472*3600);
%     ast.lat(t+1)=rad2deg( (ast.s(t)*sin(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lat(1);
%     ast.lon(t+1)=rad2deg( (ast.s(t)*cos(deg2rad(ast.inc_orbit))) / (ast.h(t)+r_planet) ) - ast.lon(1);
% 
%     % calculate rates of change
%     a2 = ast.drag_coeff * ast.A(t) * rho_air * ast.v(t)^2 / (ast.m(t) * ast.rho^2)^0.33333; %ast.drag_coeff * 1/4 * 4*pi*ast.r^2 * rho_air * ast.v(t)^2 
%     gv = g0_planet / (1 + ast.h(t) / r_planet)^2;
% 
%     ah = -a2 * ast.v_horizontal(t) / ast.v(t) * cos(deg2rad(ast.AnghSide(t))) - ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) - 2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_vertical(t) * cos(deg2rad(ast.inc_orbit))- 2* om_earth*ast.v_side(t)*sin(deg2rad(ast.lat(t)));
%     av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit))+2* om_earth*cos(deg2rad(ast.lat(t)))*ast.v_side(t) * sin(deg2rad(ast.inc_orbit));% Eötvös effect
%     aSide = -a2 * ast.v_horizontal(t) / ast.v(t) * sin(deg2rad(ast.AnghSide(t))) + 2*om_earth * ast.v_horizontal(t)* sin(deg2rad(ast.lat(t))) - 2*om_earth * ast.v_vertical(t) * cos(deg2rad(ast.lat(t))) * sin(deg2rad(ast.inc_orbit));
% %     ah = -a2 * ast.v_horizontal(t) / ast.v(t) - 2 * ast.v_vertical(t) * ast.v_horizontal(t)/ (r_planet + ast.h(t)) + om_earth*sin(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t);
% %     av = -gv - a2 * ast.v_vertical(t) / ast.v(t) + ast.v_horizontal(t) * ast.v_horizontal(t) / (r_planet + ast.h(t)) + 2*om_earth*cos(deg2rad(ast.lat(t+1)))*ast.v_horizontal(t)*cos(deg2rad(ast.inc_orbit));% Eötvös effect
% 
%     ast.zd(t+1) = rad2deg(atan2(ast.v_vertical(t), ast.v_horizontal(t)));
%     ast.AnghSide(t+1) = rad2deg(atan2(ast.side(t),ast.s(t)));
%     ast.inc_orbit=init_inc- ast.AnghSide(t+1);
%     control_inc(t)=ast.inc_orbit;
% 
%     ast.s(t+1) = ast.s(t) + ast.v_horizontal(t) * dt; %* r_planet / (r_planet + ast.h(t)) ;
%     ast.h(t+1) = ast.h(t) + ast.v_vertical(t) * dt;
%     ast.side(t+1) = ast.side(t) + ast.v_side(t) * dt;
% 
%     ast.Up_pos(t+1)=(ast.h(t+1) - ast.z_init);
%     ast.Nort_pos(t+1)=  (ast.s(t+1)*cos(deg2rad(-ast.inc_orbit)) + ast.side(t+1)*cos(deg2rad(ast.inc_orbit)) + ast.x_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lon(t+1)-OBS_lon(ii));% 
%     ast.West_pos(t+1)=  (ast.s(t+1)*sin(deg2rad(ast.inc_orbit)) + ast.side(t+1)*sin(deg2rad(ast.inc_orbit)) + ast.y_init);%(ast.h(t+1)+r_planet)*deg2rad(ast.lat(t+1)-OBS_lat(ii));% 
%     h_dist(t+1) = sqrt((ast.West_pos(t+1))^2 + (ast.Nort_pos(t+1))^2 + (ast.Up_pos(t+1))^2 );  % distance respec to the observer
% 
%     ast.v_vertical(t+1) = ast.v_vertical(t) + av * dt;
%     ast.v_side(t+1) = ast.v_side(t) + aSide * dt;
%     ast.v_horizontal(t+1) = ast.v_horizontal(t) + ah * dt;
%     ast.v(t+1) = sqrt((ast.v_horizontal(t))^2 + ast.v_vertical(t)^2 + ast.v_side(t)); %-ast.Earth_atm_speed*ast.v(t)/ast.v(1)
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
%             [~,col,~] = find(1:100==round(ast.h(t+1)/1000));
%             if ~isempty(col)
%                 magAst(col)=ast.mv(t+1,1);
%                 magFile(col)=mean(files(ii).mag((col==round(files(ii).alt(1:end-1))),1));
%             end
%             
% 
%             
%             %% cases
%             t = t+1;
% 
%              if toc > 20
%                  % stop simulation if it takes too much
%                  disp('too much time')
%                 break
%              end
%         end
%         if ast.m(end) < 0.001
%             ast.m(end)=0;
%         end
% 
% % %         sec = seconds(1*dt:dt:t*dt);
%         sec = seconds([0 1*dt:dt:(t-1)*dt]);
%         sec.Format = 'mm:ss';
%         
%         files(ii).mag_100km=magFile;
%         files(ii).mag_theory_100km=magAst;
%         files(ii).RMS=sqrt( mean( ( magFile(~isnan(magFile))-magAst(~isnan(magFile)) ).^2 ) );
% 
%         clear ast t
%     end
% end
% 
% save('meteorMassmeanNODrag.mat','files')

%% What to PLOT?

%  load('meteorMassPrecise.mat') %BEST almost same copy below 1
% load('meteorMassMEAN.mat') %better mean of all
% load('meteorMassMEANwithZEROs.mat') %WORSE mean only without zeros
% load('meteorMassGaussRMS.mat') %vary depend if H is simil Gaussian
% load('meteorMassOldRMS.mat') %LEO is way too high
% load('meteorMassNO2.mat') %LEO is way too high
load('meteorMassmeanNO2.mat') %LEO is way too high

%% PLOT

y = [1 2 3 4 5 6]; % cases
z = [nan nan nan nan nan nan nan; %material
     nan nan nan nan nan nan nan; 
     nan nan nan nan nan nan nan;
     nan nan nan nan nan nan nan;
     nan nan nan nan nan nan nan;
     nan nan nan nan nan nan nan];
zRMS{1}=z;
zRMS{2}=z;
zRMS{3}=z;

casesNum=1;

for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))
        if files(ii).init_size==5
            sizeNum=1;
        elseif files(ii).init_size==10
            sizeNum=2;
        elseif files(ii).init_size==20
            sizeNum=3;
        end
        
        if casesNum==7
            casesNum=1;
        end
        
        switch files(ii).material
            case "iron"
                materNum=1;
            case "basaltmax"
                materNum=2;
            case "minbasalt"
                materNum=3;
            case "ordchon"
                materNum=4;
            case "carchon"
                materNum=5;
            case "granite"
                materNum=6;
            case "sandstone"
                materNum=7;
        end
        zRMS{sizeNum}(casesNum,materNum)=files(ii).mass(end);
        casesNum=casesNum+1;
    end
end

for ii=1:3
    figure(ii)
    bar3(y,zRMS{ii})

    zRMS{ii}

%     zlabel('RMS')
    zlabel('Final Mass / kg')
%     zlim([0 4])
    xticklabels({'Iron','Acidic Basalt','Basic Basalt','Ordinary Chondrite','Carbonaceous Chondrite','Granite','Sandstone'})
%     if ii==1
%         xticklabels({'Iron 0.51 kg','Acidic Basalt 0.15 kg','Basic Basalt 0.20 kg','Ordinary Chondrite 0.22 kg','Carbonaceous Chondrite 0.18 kg','Granite 0.17 kg','Sandstone 0.13 kg'})
%     elseif ii==2
%         xticklabels({'Iron 4.07 kg','Acidic Basalt 1.24 kg','Basic Basalt 1.60 kg','Ordinary Chondrite 1.81 kg','Carbonaceous Chondrite 1.44 kg','Granite 1.42 kg','Sandstone 1.03 kg'})
%     elseif ii==3
%         xticklabels({'Iron 32.57 kg','Acidic Basalt 9.93 kg','Basic Basalt 12.83 kg','Ordinary Chondrite 14.48 kg','Carbonaceous Chondrite 11.59 kg','Granite 11.38 kg','Sandstone 8.27 kg'})
%     end
    yticklabels({'7.6 km/s 2° LEO','7.97 km/s 1.83° SSO','8.5 km/s 5° SSO','8 km/s 5° SSO','9 km/s 10° SSO','9 km/s 5° SSO'})
%     if ii==1
%         subtitle('RMS SCARAB & Theory')
%         title('Meteor of 5 cm ')
%     elseif ii==2
%         subtitle('RMS SCARAB & Theory')
%         title('Meteor of 10 cm')
%     elseif ii==3
%         subtitle('RMS SCARAB & Theory')
%         title('Meteor of 20 cm ')
%     end
    if ii==1
        title('Final mass SCARAB simulations 5 cm')
        subtitle('Initial mass : Iron 0.51 kg | Acidic Basalt 0.15 kg | Basic Basalt 0.20 kg | Ordinary Chondrite 0.22 kg | Carbonaceous Chondrite 0.18 kg | Granite 0.17 kg Sandstone 0.13 kg')
    elseif ii==2
        title('Final mass SCARAB simulations 10 cm')
        subtitle('Initial mass : Iron 4.07 kg | Acidic Basalt 1.24 kg | Basic Basalt 1.60 kg | Ordinary Chondrite 1.81 kg | Carbonaceous Chondrite 1.44 kg | Granite 1.42 kg | Sandstone 1.03 kg')
    elseif ii==3
        title('Final mass SCARAB simulations 20 cm')
        subtitle('Initial mass : Iron 32.57 kg | Acidic Basalt 9.93 kg | Basic Basalt 12.83 kg | Ordinary Chondrite 14.48 kg | Carbonaceous Chondrite 11.59 kg | Granite 11.38 kg | Sandstone 8.27 kg')
    end
end

%% SAVE IMAGE
% FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\3 Chap';   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
% %     fh=figure(iFig);
% %     fh.WindowState = 'maximized';
% %     set(gca,'FontSize',20)
%     saveas(gcf,fullfile(FolderName, ['RMSD_OLD_',num2str(iFig),'_Meteors.jpg']))
% end