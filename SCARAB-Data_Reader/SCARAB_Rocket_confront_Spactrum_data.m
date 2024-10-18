clc
clear all
close all
%%
load('meteorRocketMass.mat')
load('SpectrumRocketData.mat')

altitude_consider=100;

for jj=1:length(EinSteindeorbit500kmSSO(:,1))
    rockSpect_sec(jj,1)=mean(rocket.sec(round(rocket.alt)==round(EinSteindeorbit500kmSSO(jj,2))));
    rockSpect_alt(jj,1)=mean(rocket.alt(round(rocket.alt)==round(EinSteindeorbit500kmSSO(jj,2))));
    rockSpect_vel(jj,1)=mean(rocket.vel(round(rocket.alt)==round(EinSteindeorbit500kmSSO(jj,2))))/1000;
    rockSpect_ang(jj,1)=mean(rocket.flight_ang(round(rocket.alt)==round(EinSteindeorbit500kmSSO(jj,2))));
end

RMS_sec=sqrt( mean( ( rockSpect_sec - EinSteindeorbit500kmSSO(:,1) ).^2 ) );
RMS_alt=sqrt( mean( ( rockSpect_alt - EinSteindeorbit500kmSSO(:,2) ).^2 ) );
RMS_vel=sqrt( mean( ( rockSpect_vel - EinSteindeorbit500kmSSO(:,3) ).^2 ) );
RMS_ang=sqrt( mean( ( rockSpect_ang -  EinSteindeorbit500kmSSO(:,4) ).^2 ) );

RMS_sec100km=sqrt( mean( ( rockSpect_sec(rockSpect_alt>altitude_consider) - EinSteindeorbit500kmSSO((rockSpect_alt>altitude_consider),1) ).^2 ) );
RMS_alt100km=sqrt( mean( ( rockSpect_alt(rockSpect_alt>altitude_consider) - EinSteindeorbit500kmSSO((rockSpect_alt>altitude_consider),2) ).^2 ) );
RMS_vel100km=sqrt( mean( ( rockSpect_vel(rockSpect_alt>altitude_consider) - EinSteindeorbit500kmSSO((rockSpect_alt>altitude_consider),3) ).^2 ) );
RMS_ang100km=sqrt( mean( ( rockSpect_ang(rockSpect_alt>altitude_consider) -  EinSteindeorbit500kmSSO((rockSpect_alt>altitude_consider),4) ).^2 ) );

rockSpect_sec = seconds(rockSpect_sec);
rockSpect_sec.Format = 'mm:ss'; 

secSpectrum=seconds(EinSteindeorbit500kmSSO(:,1));
secSpectrum.Format = 'mm:ss'; 

figure
subplot(2,1,1)
bar([RMS_alt RMS_vel RMS_ang],'g')
ylabel('standard deviation \sigma')
% ylim([0 0.5])
xticklabels({'\sigma Altitude','\sigma Velocity','\sigma Zenith Ang.'})
grid on
title('Standard deviation from Isar Aerospace Data')


subplot(2,1,2)
bar([RMS_alt100km RMS_vel100km RMS_ang100km],'g')
ylabel('standard deviation \sigma')
% ylim([0 0.5])
xticklabels({'\sigma Altitude','\sigma Velocity','\sigma Zenith Ang.'})
grid on
title('Standard deviation from Isar Aerospace Data up to 150 km')


figure
plot(EinSteindeorbit500kmSSO(:,3),EinSteindeorbit500kmSSO(:,2),'r','LineWidth',1.2)
hold on
xlabel('Velocity [km/s]')
ylabel('Altitude [km]')
plot(rockSpect_vel,rockSpect_alt,'b')
% plot(rocket.vel(rocket.alt>38)/1000,rocket.alt(rocket.alt>38),'r')
legend('Isar Aerospace Spectrum','SCARAB Rocket','Location','southwest')
title('Velocity of SCARAB and Isar Aerospace Simulations')
grid on

figure
plot(90+EinSteindeorbit500kmSSO(:,4),EinSteindeorbit500kmSSO(:,2),'r','LineWidth',1.2)
hold on
xlabel('Zenith Ang.[°]')
ylabel('Altitude [km]')
plot(90+rockSpect_ang,rockSpect_alt,'b')
% plot(90+rocket.flight_ang(rocket.alt>38),rocket.alt(rocket.alt>38),'r')
legend('Isar Aerospace Spectrum','SCARAB Rocket','Location','southwest')
title('Zenith Angle of SCARAB and Isar Aerospace Simulations')
grid on

figure
plot(secSpectrum,EinSteindeorbit500kmSSO(:,2),'r','LineWidth',1.2)
hold on
xlabel('Time [mm:ss]')
ylabel('Altitude [km]')
plot(rockSpect_sec,rockSpect_alt,'b')
% plot(rocket.vel(rocket.alt>38)/1000,rocket.alt(rocket.alt>38),'r')
legend('Isar Aerospace Spectrum','SCARAB Rocket','Location','southwest')
title('Altitude of SCARAB and Isar Aerospace Simulations')
grid on

figure
plot(secSpectrum,EinSteindeorbit500kmSSO(:,3),'r','LineWidth',1.2)
hold on
xlabel('Time [mm:ss]')
ylabel('Velocity [km/s]')
plot(rockSpect_sec,rockSpect_vel,'b')
% plot(rocket.vel(rocket.alt>38)/1000,rocket.alt(rocket.alt>38),'r')
legend('Isar Aerospace Spectrum','SCARAB Rocket','Location','southwest')
title('Velocity of SCARAB and Isar Aerospace Simulations')
grid on

figure
plot(secSpectrum,90+EinSteindeorbit500kmSSO(:,4),'r','LineWidth',1.2)
hold on
xlabel('Time [mm:ss]')
ylabel('Zenith Ang.[°]')
plot(rockSpect_sec,90+rockSpect_ang,'b')
% plot(90+rocket.flight_ang(rocket.alt>38),rocket.alt(rocket.alt>38),'r')
legend('Isar Aerospace Spectrum','SCARAB Rocket','Location','southwest')
title('Zenith Angle of SCARAB and Isar Aerospace Simulations')
grid on

%%
% FolderName = 'C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Reports\Images\Oct';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    figure(iFig)
%     saveas(gcf,fullfile(FolderName, [num2str(iFig),'_Theory_diffInclFlightAng.jpg']))
set(gca,'FontSize',20)
end