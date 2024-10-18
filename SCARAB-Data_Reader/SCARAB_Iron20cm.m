clc
clear all
close all

%%
fileroot='C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\SCARAB\';
files= dir(fullfile([append(fileroot,'Meteors Iron 20cm')]));
%files=dir(fullfile('C:\Users\maxivy\Documents\DRAMA\TEST folder','*.txt'));


material_name=["iron" "basaltmax" "minbasalt" "ordchon" "carchon" "granite" "sandstone"];
material_density=[7870.0 2400.0 3100.0 3500.0 2800.0 2750.0 2000.0];
material_meltingheat=[272000.0 400000.0 506000.0 265000.0 265000.0 250000.0 680000.0];
material_color=["g" "b" "c" "r" "m" "k" "y"];
lum_eff = [0.1 0.001];

[X,Y] = meshgrid([7 7.5 8 8.5 9 10],[1 2 3 5 10 15 45 90]);
angles=[1 2 3 5 10 15 45 90];

x=[];
y=[];
z=[];

mass_loss_inc=[];
inc=[];

Z=zeros(8,6);

matrix_holder=[];
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))        
        %% metadata info
        materia_unit=regexp(files(ii).name,'([a-z]+)','match');
        files(ii).material=materia_unit{2};
        files(ii).rho=material_density(material_name==files(ii).material);
        files(ii).Hmelt=material_meltingheat(material_name==files(ii).material);
        
        metadata=str2double(regexp(files(ii).name,'([ \d\.]+)','match'));
        files(ii).init_size=metadata(1);
        files(ii).init_alt=metadata(2);
        files(ii).init_vel=metadata(3);
        files(ii).init_ang=metadata(4);
        files(ii).init_inc=metadata(5);
        
        %% mass info
        filename=append(files(ii).folder,'/',files(ii).name,'/0.1/mas.hst');
        fid = fopen( filename );
        lineRead = fgets(fid);
        
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)','tokens'))
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    if string(cell_holder{1}(1))~='#'
                        matrix_holder=[matrix_holder; (str2double(cell_holder(1:4)))];
                    end
                end
            lineRead = fgets(fid);
        end
        files(ii).sec= matrix_holder(:,1);
        sec = seconds(matrix_holder(:,1));
        sec.Format = 'mm:ss.SSS'; 
        files(ii).time= sec;
        files(ii).alt= matrix_holder(:,2);
        files(ii).gr_trk= matrix_holder(:,3);
        files(ii).mass= matrix_holder(:,4);
        files(ii).mass_rate=[0;abs((files(ii).mass(1:end-1)-files(ii).mass(2:end))./(files(ii).sec(1:end-1)-files(ii).sec(2:end)))];
        files(ii).vel=1000.*[files(ii).init_vel;abs(sqrt((files(ii).alt(1:end-1)-files(ii).alt(2:end)).^2+(files(ii).gr_trk(1:end-1)-files(ii).gr_trk(2:end)).^2)./(files(ii).sec(1:end-1)-files(ii).sec(2:end)))];
        files(ii).crosSec=pi*(files(ii).mass*3/4 / (pi*files(ii).rho)).^(2/3);
        files(ii).crosSec(files(ii).crosSec==0)=min(files(ii).crosSec(files(ii).crosSec>0));
        
        filename=append(files(ii).folder,'/',files(ii).name,'/0.1/atm.hst');
        fid = fopen( filename );
        lineRead = fgets(fid);
        matrix_holder=[];
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)','tokens'))
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    if string(cell_holder{1}(1))~='#'
                        matrix_holder=[matrix_holder; (str2double(cell_holder(1:6)))];
                    end
                end
            lineRead = fgets(fid);
        end
        for jj=1:length(files(ii).sec)
            files(ii).rho_air(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),5));
        end
        
        files(ii).heat_transf_coef=files(ii).mass_rate.*(2*files(ii).Hmelt) ./ (files(ii).crosSec .* files(ii).rho_air .* (files(ii).vel).^3);
        files(ii).heat_transf_coef(files(ii).heat_transf_coef>1)=1;      
        files(ii).lum_ene = 0.5 .* files(ii).vel.^2 .* files(ii).mass_rate.* lum_eff ;%.* 1e10 ./ (dist.^2); lum_eff
        files(ii).lum_ene(files(ii).lum_ene==0)=nan;
        files(ii).mag = 6.8 - 1.086 * log(files(ii).lum_ene);
        matrix_holder=[];
            
%% plot 
        if files(ii).init_vel==7.97
            files(ii).init_vel=8;
            files(ii).init_ang=2;
        end
        
        if files(ii).init_vel==7.6
            files(ii).init_vel=7.5;
        end
        
        if files(ii).init_ang==3
            files(ii).init_vel=8.3;
            files(ii).init_ang=4;
        end
            
        if files(ii).init_inc==97.5 %&& files(ii).init_ang~=3
            switch files(ii).init_vel
                case 7
                    Z(angles==files(ii).init_ang,1)=files(ii).mass(1)-files(ii).mass(end);
                case 7.5
                    Z(angles==files(ii).init_ang,2)=files(ii).mass(1)-files(ii).mass(end);
                case 8
                    Z(angles==files(ii).init_ang,3)=files(ii).mass(1)-files(ii).mass(end);
                case 8.5
                    Z(angles==files(ii).init_ang,4)=files(ii).mass(1)-files(ii).mass(end);
                case 9
                    Z(angles==files(ii).init_ang,5)=files(ii).mass(1)-files(ii).mass(end);
                case 10
                    Z(angles==files(ii).init_ang,6)=files(ii).mass(1)-files(ii).mass(end);
            end
            x=[x; files(ii).init_vel];
            y=[y; files(ii).init_ang];
            z=[z; files(ii).mass(1)-files(ii).mass(end)];
            
            
        end
        
        
            if files(ii).init_vel==7.5 && files(ii).init_ang==2
                mass_loss_inc=[mass_loss_inc files(ii).mass(1)-files(ii).mass(end)];
                inc=[inc files(ii).init_inc];
            end
        
    else

    end
    
end

x=[x; 8.5; 7.5 ; 9.5; 10];
y=90-[y; 1; 1; 1; 6];
z=[z; 0; 0; 0; 20];

% Z(3,1)=mean([Z(2,1) Z(4,1)]);
% Z(5,1)=mean([Z(4,1) Z(7,1)]);
% Z(6,1)=mean([Z(5,1) Z(7,1)]);
% 
% Z(5,3)=mean([Z(4,3) Z(7,3)]);
% Z(6,3)=mean([Z(5,3) Z(7,3)]);

figure
bar(inc,mass_loss_inc)
hold on
q = fit(inc',mass_loss_inc','cubicinterp');%'cubicinterp');
plot( q, inc,mass_loss_inc)
ylabel('Mass Loss [kg]')
xlabel('Orbit Inclination Angle [°]')
grid on
title('Orbit Inclination Dependence')

figure
f = fit([x, y],z,'cubicinterp');%'cubicinterp');
plot( f, [x, y],z )
% 
% figure
% surf(X,Y,Z)
% ii=1;
% % for ii=1:2
%     figure(ii)
    zlim([0 32.5])
    ylim([90-20 90])
    zlabel('Mass Loss [kg]') 
    xlabel('Velocity [km/s]')
    ylabel('Zenith Angle [°]')
% end
% save('meteorIronMass.mat','files')
