clc
clear all
close all

%%

load('sat1_6.3kms.mat')
load('sat2_6.7kms.mat')
load('sat3_7.1kms.mat')
load('sat4_7.5kms.mat')

allSat_fragm{1}=sat1;
allSat_fragm{2}=sat2;
allSat_fragm{3}=sat3;
allSat_fragm{4}=sat4;

material_name=["iron" "basaltmax" "minbasalt" "ordchon" "carchon" "granite" "sandstone"];
material_density=[7870.0 2400.0 3100.0 3500.0 2800.0 2750.0 2000.0];
material_meltingheat=[272000.0 400000.0 506000.0 265000.0 265000.0 250000.0 680000.0];
material_color=["g" "b" "c" "r" "m" "k" "y"];

altitude=1:100;
lum_eff = [0.1 0.001];

%% Saves
% % "iron" 1 poly22
% % files(ii).init_alt~=120 && files(ii).init_ang~=1.83
% velAvg=9E3;
% r=0.01;
% H_melt=272000.0;

% % all
% iron
% velAvg=8E3;

% % all
% basaltmax
% velAvg=8E3;
% r=0.05;

% % no 120
% material="minbasalt";
% velAvg=8E3;
% r=0.1;

% % all
% material="iron";
% velAvg=6E3;
% r=0.05;
% H_melt=material_meltingheat(find(material_name==material));

% material="iron";
% velAvg=7.5E3;
% r=0.01;
% H_melt=material_meltingheat(find(material_name==material));
% colorMat=material_color(find(material_name==material));

%% plots
material=["iron" "iron" "iron" "iron"];
% velAvgVEC=[6.3E3 6.7E3 7.1E3 7.5E3];
% rVEC=[0.037 0.015 0.02 0.01];

velAvgVEC=[allSat_fragm{1}(1,5) allSat_fragm{2}(1,5) allSat_fragm{3}(1,5) allSat_fragm{4}(1,5)]*1000;
rVEC=[0.035 0.017 0.02 0.01];

figure
for ii=1:4    
    subplot(2,2,ii)
    plot(allSat_fragm{ii}(:,3),allSat_fragm{ii}(:,2),'.b')
    hold on
    
    velAvg=velAvgVEC(ii);
    r=rVEC(ii);
    H_melt=400000;%material_meltingheat(find(material_name==material(ii)));
    colorMat=material_color(find(material_name==material(ii)));
    
    [rho_air,~,~,~,~,~,~] = atmos(altitude*1000);
    heatCoef=Gauss_HeatCoef(altitude,velAvg/1000,r*100,material);
    mass_rateAvg = heatCoef.*((pi*r^2) .* rho_air' .* (velAvg).^3) ./ (2*H_melt);
    lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff;
    lum_eneAvg(lum_eneAvg==0)=nan;
    magAvg = 6.8 - 1.086 * log(lum_eneAvg);
    plot(magAvg(:,1),altitude,'r','Linewidth',1.2)
    plot(magAvg(:,2),altitude,'r--','Linewidth',1.2)
    
    ylim([0 100])
    xlim([-10 10])
    xlabel('Magnitude') 
    ylabel('Altitude [km]')
    set(gca, 'XDir','reverse')
    grid on
    title(['Fragment ' num2str(ii)])
    sgtitle('AMOS Long March 3B Re-entry') %CZ-3B R/B
end

material=["ordchon" "ordchon" "ordchon" "ordchon"];
% velAvgVEC=[6.3E3 6.7E3 7.1E3 7.5E3];
velAvgVEC=[allSat_fragm{1}(1,5) allSat_fragm{2}(1,5) allSat_fragm{3}(1,5) allSat_fragm{4}(1,5)]*1000;

rVEC=[0.06 0.03 0.025 0.01];

figure
for ii=1:4    
    subplot(2,2,ii)
    plot(allSat_fragm{ii}(:,3),allSat_fragm{ii}(:,2),'.b')
    hold on
    
    velAvg=velAvgVEC(ii);
    r=rVEC(ii);
    H_melt=material_meltingheat(find(material_name==material(ii)));
    colorMat=material_color(find(material_name==material(ii)));
    
    [rho_air,~,~,~,~,~,~] = atmos(altitude*1000);
    heatCoef=Gauss_HeatCoef(altitude,velAvg/1000,r*100,material);
    mass_rateAvg = heatCoef.*((pi*r^2) .* rho_air' .* (velAvg).^3) ./ (2*H_melt);
    lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff;
    lum_eneAvg(lum_eneAvg==0)=nan;
    magAvg = 6.8 - 1.086 * log(lum_eneAvg);
    plot(magAvg(:,1),altitude,'g','Linewidth',1.2)
    plot(magAvg(:,2),altitude,'g--','Linewidth',1.2)
    
    ylim([0 100])
    xlim([-10 10])
    xlabel('Magnitude') 
    ylabel('Altitude [km]')
    set(gca, 'XDir','reverse')
    grid on
    title(['Fragment ' num2str(ii)])
    sgtitle('AMOS Long March 3B Re-entry - Ordinary Chondrite') %CZ-3B R/B
end

% load('meteorMass.mat')
% posit5=72;
% heatCoef=Gauss_HeatCoef(files(posit5).alt,files(posit5).init_vel,files(posit5).init_size,files(posit5).material);

% for ii=1:4    
%     subplot(2,2,ii)
%     plot(allSat_fragm{ii}(:,3),allSat_fragm{ii}(:,2),'.b')
%     hold on
% %     plot(files(posit5).mag(:,1),files(posit5).alt,'Linewidth',1.2)
% %     plot(files(posit5).mag(:,2),files(posit5).alt,'--','Linewidth',1.2)
% 
%     velAvg=velAvgVEC(ii);
%     r=rVEC(ii);
%     H_melt=material_meltingheat(find(material_name==material(ii)));
%     colorMat=material_color(find(material_name==material(ii)));
%     
%     [rho_air,~,~,~,~,~,~] = atmos(altitude*1000);
%     heatCoef=Gauss_HeatCoef(altitude,velAvg/1000,r*100,material);
%     mass_rateAvg = heatCoef.*((pi*r^2) .* rho_air' .* (velAvg).^3) ./ (2*H_melt);
%     lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff;
%     lum_eneAvg(lum_eneAvg==0)=nan;
%     magAvg = 6.8 - 1.086 * log(lum_eneAvg);
%     plot(magAvg(:,1),altitude,colorMat,'Linewidth',1.2)
%     plot(magAvg(:,2),altitude,append(colorMat,'--'),'Linewidth',1.2)
% %     
% %     mass_rateAvg = files(posit5).heat_transf_coef.*(files(posit5).crosSec .* files(posit5).rho_air .* (files(posit5).vel).^3) ./ (2*files(posit5).Hmelt);
% %     lum_eneAvg = 0.5 .* files(posit5).vel.^2 .* mass_rateAvg .* lum_eff ;
% %     lum_eneAvg(lum_eneAvg==0)=nan;
% %     magAvg = 6.8 - 1.086 * log(lum_eneAvg);    
% %     plot(magAvg(:,1),files(posit5).alt,'r','Linewidth',1.2)
% %     plot(magAvg(:,2),files(posit5).alt,'r--','Linewidth',1.2)
%     
% %     mass_rateAvg = files(posit5).heat_transf_coef.*(files(posit5).crosSec .* files(posit5).rho_air .* (velAvg).^3) ./ (2*files(posit5).Hmelt);
% %     lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff ;
% %     lum_eneAvg(lum_eneAvg==0)=nan;
% %     magAvg = 6.8 - 1.086 * log(lum_eneAvg);    
% %     plot(magAvg(:,1),files(posit5).alt,'r','Linewidth',1.2)
% %     plot(magAvg(:,2),files(posit5).alt,'r--','Linewidth',1.2)
%     
% %     mass_rateAvg = files(posit5).coef_gauss_avg.*(crosSecAvg .* rho_airAvg .* (velAvg).^3) ./ (2*files(ii).Hmelt);
% %     lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff ;
% %     lum_eneAvg(lum_eneAvg==0)=nan;
% %     magAvg = 6.8 - 1.086 * log(lum_eneAvg);
% % 
% %     plot(magAvg(:,1),altitude,'r','Linewidth',1.2)
% %     plot(magAvg(:,2),altitude,'r--','Linewidth',1.2)
%     
%     ylim([0 100])
%     xlim([-10 10])
%     xlabel('Magnitude') 
%     ylabel('Altitude [km]')
%     set(gca, 'XDir','reverse')
%     grid on
%     title(['Frafgment ' num2str(ii)])
%     sgtitle('Long March 3B Re-entry') %CZ-3B R/B
% end
% 
