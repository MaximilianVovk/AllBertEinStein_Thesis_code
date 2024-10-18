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

altitude=1:100;

for ii=1:4    
    subplot(2,2,ii)
    plot(allSat_fragm{ii}(:,3),allSat_fragm{ii}(:,2),'.b')
    hold on
    plot(Gauss_HeatCoef(altitude,6.3,5,"iron"),altitude,'Linewidth',1.2)
    
    ylim([0 100])
    xlim([-10 10])
    xlabel('Magnitude') 
    ylabel('Altitude [km]')
    set(gca, 'XDir','reverse')
    grid on
    title(['Frafgment ' num2str(ii)])
end

