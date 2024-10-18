clc
clear all
close all

%%
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

[X,Y] = meshgrid([5 10 20],1:100);
Z=zeros(100,3);
load('meteorMass.mat')

[X,Y] = meshgrid([7 7.5 8 8.5 9 10],[1 2 3 5 10 15 45 90]);
angles=[1 2 3 5 10 15 45 90];


iron_vel=[];
iron_size=[];
iron_a1=[];
iron_b1=[];
iron_c1=[];

basaltmax_vel=[];
basaltmax_size=[];
basaltmax_a1=[];
basaltmax_b1=[];
basaltmax_c1=[];

basaltmin_vel=[];
basaltmin_size=[];
basaltmin_a1=[];
basaltmin_b1=[];
basaltmin_c1=[];

ordchon_vel=[];
ordchon_size=[];
ordchon_a1=[];
ordchon_b1=[];
ordchon_c1=[];

carchon_vel=[];
carchon_size=[];
carchon_a1=[];
carchon_b1=[];
carchon_c1=[];

granite_vel=[];
granite_size=[];
granite_a1=[];
granite_b1=[];
granite_c1=[];

sandstone_vel=[];
sandstone_size=[];
sandstone_a1=[];
sandstone_b1=[];
sandstone_c1=[];



Z=zeros(8,6);

matrix_holder=[];
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))         
%% plot     
% if files(ii).init_ang~=10 && files(ii).init_alt~=120 && files(ii).init_ang~=1.83
% if files(ii).init_alt~=120 && files(ii).init_ang~=1.83
% if files(ii).init_alt~=120 && files(ii).init_ang~=1.83 && files(ii).init_vel~=8
% if files(ii).init_ang~=1.83
% if files(ii).init_alt~=120
if files(ii).init_ang~=1.888

        switch files(ii).material
            case "iron"
                iron_vel=[iron_vel; files(ii).init_vel];
                iron_size=[iron_size; files(ii).init_size];
                iron_a1=[iron_a1; files(ii).gauss_StD_Mean.a1];
                iron_b1=[iron_b1; files(ii).gauss_StD_Mean.b1];
                iron_c1=[iron_c1; files(ii).gauss_StD_Mean.c1];
            case "basaltmax"
                basaltmax_vel=[basaltmax_vel; files(ii).init_vel];
                basaltmax_size=[basaltmax_size; files(ii).init_size];
                basaltmax_a1=[basaltmax_a1; files(ii).gauss_StD_Mean.a1];
                basaltmax_b1=[basaltmax_b1; files(ii).gauss_StD_Mean.b1];
                basaltmax_c1=[basaltmax_c1; files(ii).gauss_StD_Mean.c1];
            case "minbasalt"
                basaltmin_vel=[basaltmin_vel; files(ii).init_vel];
                basaltmin_size=[basaltmin_size; files(ii).init_size];
                basaltmin_a1=[basaltmin_a1; files(ii).gauss_StD_Mean.a1];
                basaltmin_b1=[basaltmin_b1; files(ii).gauss_StD_Mean.b1];
                basaltmin_c1=[basaltmin_c1; files(ii).gauss_StD_Mean.c1];
            case "ordchon"
                ordchon_vel=[ordchon_vel; files(ii).init_vel];
                ordchon_size=[ordchon_size; files(ii).init_size];
                ordchon_a1=[ordchon_a1; files(ii).gauss_StD_Mean.a1];
                ordchon_b1=[ordchon_b1; files(ii).gauss_StD_Mean.b1];
                ordchon_c1=[ordchon_c1; files(ii).gauss_StD_Mean.c1];
            case "carchon"
                carchon_vel=[carchon_vel; files(ii).init_vel];
                carchon_size=[carchon_size; files(ii).init_size];
                carchon_a1=[carchon_a1; files(ii).gauss_StD_Mean.a1];
                carchon_b1=[carchon_b1; files(ii).gauss_StD_Mean.b1];
                carchon_c1=[carchon_c1; files(ii).gauss_StD_Mean.c1];
            case "granite"
                granite_vel=[granite_vel; files(ii).init_vel];
                granite_size=[granite_size; files(ii).init_size];
                granite_a1=[granite_a1; files(ii).gauss_StD_Mean.a1];
                granite_b1=[granite_b1; files(ii).gauss_StD_Mean.b1];
                granite_c1=[granite_c1; files(ii).gauss_StD_Mean.c1];
            case "sandstone"
                sandstone_vel=[sandstone_vel; files(ii).init_vel];
                sandstone_size=[sandstone_size; files(ii).init_size];
                sandstone_a1=[sandstone_a1; files(ii).gauss_StD_Mean.a1];
                sandstone_b1=[sandstone_b1; files(ii).gauss_StD_Mean.b1];
                sandstone_c1=[sandstone_c1; files(ii).gauss_StD_Mean.c1];
        end
        
end
    else

    end
    
end

TypeFit='poly22';%'poly22'; %linearinterp

coefGauss.iron_a1 = fit([iron_vel, iron_size],iron_a1,TypeFit);
coefGauss.iron_b1 = fit([iron_vel, iron_size],iron_b1,TypeFit);
coefGauss.iron_c1 = fit([iron_vel, iron_size],iron_c1,TypeFit);

coefGauss.basaltmax_a1 = fit([basaltmax_vel, basaltmax_size],basaltmax_a1,TypeFit);
coefGauss.basaltmax_b1 = fit([basaltmax_vel, basaltmax_size],basaltmax_b1,TypeFit);
coefGauss.basaltmax_c1 = fit([basaltmax_vel, basaltmax_size],basaltmax_c1,TypeFit);

coefGauss.basaltmin_a1 = fit([basaltmin_vel, basaltmin_size],basaltmin_a1,TypeFit);
coefGauss.basaltmin_b1 = fit([basaltmin_vel, basaltmin_size],basaltmin_b1,TypeFit);
coefGauss.basaltmin_c1 = fit([basaltmin_vel, basaltmin_size],basaltmin_c1,TypeFit);

coefGauss.ordchon_a1 = fit([ordchon_vel, ordchon_size],ordchon_a1,TypeFit);
coefGauss.ordchon_b1 = fit([ordchon_vel, ordchon_size],ordchon_b1,TypeFit);
coefGauss.ordchon_c1 = fit([ordchon_vel, ordchon_size],ordchon_c1,TypeFit);

coefGauss.carchon_a1 = fit([carchon_vel, carchon_size],carchon_a1,TypeFit);
coefGauss.carchon_b1 = fit([carchon_vel, carchon_size],carchon_b1,TypeFit);
coefGauss.carchon_c1 = fit([carchon_vel, carchon_size],carchon_c1,TypeFit);

coefGauss.granite_a1 = fit([granite_vel, granite_size],granite_a1,TypeFit);
coefGauss.granite_b1 = fit([granite_vel, granite_size],granite_b1,TypeFit);
coefGauss.granite_c1 = fit([granite_vel, granite_size],granite_c1,TypeFit);

coefGauss.sandstone_a1 = fit([sandstone_vel, sandstone_size],sandstone_a1,TypeFit);
coefGauss.sandstone_b1 = fit([sandstone_vel, sandstone_size],sandstone_b1,TypeFit);
coefGauss.sandstone_c1 = fit([sandstone_vel, sandstone_size],sandstone_c1,TypeFit);

save('meteorGaussCoef.mat','coefGauss')

figure(1)
subplot(1,3,1)
plot( coefGauss.iron_a1, [iron_vel, iron_size],iron_a1 )
title('a_1 Gauss Coef.')
zlabel('a_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,2)
plot( coefGauss.iron_b1, [iron_vel, iron_size],iron_b1 )
title('b_1 Gauss Coef.')
zlabel('b_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,3)
plot( coefGauss.iron_c1, [iron_vel, iron_size],iron_c1 )
title('c_1 Gauss Coef.')
zlabel('c_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
sgtitle('Meteor Iron')

figure(2)
subplot(1,3,1)
plot( coefGauss.basaltmax_a1, [basaltmax_vel, basaltmax_size],basaltmax_a1 )
title('a_1 Gauss Coef.')
zlabel('a_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,2)
plot( coefGauss.basaltmax_b1, [basaltmax_vel, basaltmax_size],basaltmax_b1 )
title('b_1 Gauss Coef.')
zlabel('b_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,3)
plot( coefGauss.basaltmax_c1, [basaltmax_vel, basaltmax_size],basaltmax_c1 )
title('c_1 Gauss Coef.')
zlabel('c_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
sgtitle('Meteor Basalt Acidic')

figure(3)
subplot(1,3,1)
plot( coefGauss.basaltmin_a1 , [basaltmin_vel, basaltmin_size],basaltmin_a1 )
title('a_1 Gauss Coef.')
zlabel('a_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,2)
plot( coefGauss.basaltmin_b1 , [basaltmin_vel, basaltmin_size],basaltmin_b1 )
title('b_1 Gauss Coef.')
zlabel('b_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,3)
plot( coefGauss.basaltmin_c1 , [basaltmin_vel, basaltmin_size],basaltmin_c1 )
title('c_1 Gauss Coef.')
zlabel('c_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
sgtitle('Meteor Basalt Basic')

figure(4)
subplot(1,3,1)
plot( coefGauss.ordchon_a1, [ordchon_vel, ordchon_size],ordchon_a1 )
title('a_1 Gauss Coef.')
zlabel('a_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,2)
plot( coefGauss.ordchon_b1, [ordchon_vel, ordchon_size],ordchon_b1 )
title('b_1 Gauss Coef.')
zlabel('b_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,3)
plot( coefGauss.ordchon_c1, [ordchon_vel, ordchon_size],ordchon_c1 )
title('c_1 Gauss Coef.')
zlabel('c_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
sgtitle('Meteor Ordinary Chondrite')

figure(5)
subplot(1,3,1)
plot( coefGauss.carchon_a1, [carchon_vel, carchon_size],carchon_a1 )
title('a_1 Gauss Coef.')
zlabel('a_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,2)
plot( coefGauss.carchon_b1, [carchon_vel, carchon_size],carchon_b1 )
title('b_1 Gauss Coef.')
zlabel('b_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,3)
plot( coefGauss.carchon_c1, [carchon_vel, carchon_size],carchon_c1 )
title('c_1 Gauss Coef.')
zlabel('c_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
sgtitle('Meteor Carbonaceous Chondrite')

figure(6)
subplot(1,3,1)
plot( coefGauss.granite_a1, [granite_vel, granite_size],granite_a1 )
title('a_1 Gauss Coef.')
zlabel('a_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,2)
plot( coefGauss.granite_b1, [granite_vel, granite_size],granite_b1 )
title('b_1 Gauss Coef.')
zlabel('b_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,3)
plot( coefGauss.granite_c1, [granite_vel, granite_size],granite_c1 )
title('c_1 Gauss Coef.')
zlabel('c_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
sgtitle('Meteor Granite')

figure(7)
subplot(1,3,1)
plot( coefGauss.sandstone_a1, [sandstone_vel, sandstone_size],sandstone_a1 )
title('a_1 Gauss Coef.')
zlabel('a_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,2)
plot( coefGauss.sandstone_b1, [sandstone_vel, sandstone_size],sandstone_b1 )
title('b_1 Gauss Coef.')
zlabel('b_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
subplot(1,3,3)
plot( coefGauss.sandstone_c1, [sandstone_vel, sandstone_size],sandstone_c1 )
title('c_1 Gauss Coef.')
zlabel('c_1')
ylabel('Diameter [cm]') 
xlabel('Velocity [km/s]')
sgtitle('Meteor Sandstone')



% plot( ironFunction_a1, [iron_vel, iron_size],iron_a1 )
% 
% figure(ii)
% zlabel('Gauss Coef.') 
% xlabel('Velocity [km/s]')
% ylabel('Diameter [cm]')
