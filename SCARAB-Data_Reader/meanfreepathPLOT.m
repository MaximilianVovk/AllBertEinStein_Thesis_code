clc
clear all
close all
%%
wide=1000;

hh{1}=0:wide:11000;
hh{2}=11000:wide:20000;
hh{3}=20000:wide:47000;
hh{4}=47000:wide:51000;
hh{5}=51000:wide:61000;
hh{6}=61000:wide:71000;
hh{7}=71000:wide:81000;
hh{8}=81000:wide:120000;
hh{9}=120000:wide:160000;
hh{10}=160000:wide:200000;

figure
for ii=1:length(hh)
[rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(hh{ii});
R_specif=287;
mean_free_path=nu_air.*rho_air./P_air.*sqrt(pi.*R_specif.*T_air./2);

f{ii}=fit(h',mean_free_path','exp1');
semilogx(f{ii}(h),h/1000,'r','LineWidth',1.2);
grid on
hold on
semilogx(mean_free_path,h/1000,'b.');
xlabel('Mean Free Path [m]') 
ylabel('Altitude [km]')

clearvars -except hh f
end

% figure
% plot(f,hh,mean_free_path)
% grid on
% ylabel('Mean Free Path [m]') 
% xlabel('Altitude [m]')