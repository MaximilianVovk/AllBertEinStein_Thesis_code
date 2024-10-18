clc
clear all
close all
%% CREATE FILE LongMarch3B_Fragments.mat
run('config.m')

load('sat1_6.3kms.mat')
load('sat2_6.7kms.mat')
load('sat3_7.1kms.mat')
load('sat4_7.5kms.mat')

allSat_fragm{1}=sat1;
allSat_fragm{2}=sat2;
allSat_fragm{3}=sat3;
allSat_fragm{4}=sat4;

LongMarch3B(1).name='Maunakea OBS - Fragm 1';
LongMarch3B(2).name='Maunakea OBS - Fragm 2';
LongMarch3B(3).name='Maunakea OBS - Fragm 3';
LongMarch3B(4).name='Maunakea OBS - Fragm 4';
LongMarch3B(5).name='Haleakala OBS - Fragm 1';
LongMarch3B(6).name='Haleakala OBS - Fragm 2';
LongMarch3B(7).name='Haleakala OBS - Fragm 3';
LongMarch3B(8).name='Haleakala OBS - Fragm 4';

numTimes=0;
OBS=0;
for ii=1:8
    if ii>=5
        OBS=1;
        numTimes=1;
    end
    
    LongMarch3B(ii).OBS = OBS;
    LongMarch3B(ii).fragm = ii-(numTimes*4);
    
    if ii==4 %fragment 4 invert the 2 observatories
        OBS=1;
    end
    if ii==8
        OBS=0;
    end
    
    LongMarch3B(ii).sec = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),1);
    sec = seconds(LongMarch3B(ii).sec);
    sec.Format = 'mm:ss';
    LongMarch3B(ii).time = sec;
    LongMarch3B(ii).alt = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),2);
    LongMarch3B(ii).gr_trk = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),4);
    LongMarch3B(ii).ang = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),7);
    LongMarch3B(ii).vel = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),5);
    LongMarch3B(ii).distOBS = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),8);
    LongMarch3B(ii).mag = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),6);
    LongMarch3B(ii).ABSmag = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),3);
    LongMarch3B(ii).azimuth = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),9);
    LongMarch3B(ii).elevationORaltitude = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),10);
    LongMarch3B(ii).lon = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),11);
    LongMarch3B(ii).lat = allSat_fragm{ii-(numTimes*4)}((allSat_fragm{ii-(numTimes*4)}(:,13)==OBS),12);
    
    headAng=rad2deg( atan2( (LongMarch3B(ii).alt(1:end-1)+r_planet).*deg2rad(LongMarch3B(ii).lat(1:end-1)-LongMarch3B(ii).lat(2:end)) , (LongMarch3B(ii).alt(1:end-1)+r_planet).*deg2rad(LongMarch3B(ii).lon(1:end-1)-LongMarch3B(ii).lon(2:end)) ) );
    headAng=[headAng ; headAng(end)];
    
    LongMarch3B(ii).head_ang = headAng;
    
end

save('LongMarch3B_Fragments.mat','LongMarch3B')