clc
clear all
close all

%%
fileroot='C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\SCARAB\';
% files= dir(fullfile([append(fileroot,'Meteors')]));
files_rocket= dir(fullfile([append(fileroot,'A_1\a-1-steel-500km97.5inc')]));
%files=dir(fullfile('C:\Users\maxivy\Documents\DRAMA\TEST folder','*.txt'));

lum_eff = [0.1 0.001];

rocket.material="Rocket";
rocket.rho=2870;
rocket.Hmelt=337700; %titanium
rocket.init_size=100;
rocket.init_alt=500000;
rocket.init_vel=7.5;
rocket.init_ang=0.00242;
rocket.init_inc=97.5;
rocket.crosSec=4.5*2;
rocket.sec= [];
rocket.time= [];
rocket.alt= [];
rocket.gr_trk= [];
rocket.mass= [];
rocket.mass_rate= [];
rocket.vel= [];
rocket.rho_air= [];
rocket.heat_transf_coef= [];    
rocket.lum_ene = [];
rocket.mag = [];
rocket.flight_ang= [];
rocket.head_ang= [];
rocket.lat= [];
rocket.lon=[];


matrix_holder=[];
for ii=1:length(files_rocket)
    if files_rocket(ii).isdir==1 && ~isempty(regexp(files_rocket(ii).name,'.1','match'))        

        files_rocket(ii).material=rocket.material;
        files_rocket(ii).rho=rocket.rho;
        files_rocket(ii).Hmelt=rocket.Hmelt;
        files_rocket(ii).init_size=rocket.init_size;
        files_rocket(ii).init_alt=rocket.init_alt;
        files_rocket(ii).init_vel=rocket.init_vel;
        files_rocket(ii).init_ang=rocket.init_ang;
        files_rocket(ii).init_inc=rocket.init_inc;
        
        %% mass info
        filename=append(files_rocket(ii).folder,'\',files_rocket(ii).name,'\mas.hst');
        fid = fopen( filename );
        lineRead = fgets(fid);
        matrix_holder=[];
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)','tokens'))
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    if string(cell_holder{1}(1))~='#'
                        matrix_holder=[matrix_holder; (str2double(cell_holder(1:4)))];
                    end
                end
            lineRead = fgets(fid);
        end
        files_rocket(ii).sec= matrix_holder(:,1);
        sec = seconds(matrix_holder(:,1));
        sec.Format = 'mm:ss.SSS'; 
        files_rocket(ii).time= sec;
        files_rocket(ii).alt= matrix_holder(:,2);
        files_rocket(ii).gr_trk= matrix_holder(:,3);
        files_rocket(ii).mass= matrix_holder(:,4);
        files_rocket(ii).mass_rate=[0;abs((files_rocket(ii).mass(1:end-1)-files_rocket(ii).mass(2:end))./(files_rocket(ii).sec(1:end-1)-files_rocket(ii).sec(2:end)))];
        files_rocket(ii).vel=1000.*[files_rocket(ii).init_vel;abs(sqrt((files_rocket(ii).alt(1:end-1)-files_rocket(ii).alt(2:end)).^2+(files_rocket(ii).gr_trk(1:end-1)-files_rocket(ii).gr_trk(2:end)).^2)./(files_rocket(ii).sec(1:end-1)-files_rocket(ii).sec(2:end)))];
        files_rocket(ii).crosSec=4*2;
        
        filename=append(files_rocket(ii).folder,'\',files_rocket(ii).name,'\atm.hst');
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
        jjarray=1;
        for jj=round(files_rocket(ii).sec(1)):round(files_rocket(ii).sec(end))
            files_rocket(ii).rho_air(jjarray,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj)),5));
            jjarray=jjarray+1;
        end
        
        files_rocket(ii).heat_transf_coef=files_rocket(ii).mass_rate.*(2*files_rocket(ii).Hmelt) ./ (files_rocket(ii).crosSec .* files_rocket(ii).rho_air .* (files_rocket(ii).vel).^3);
        files_rocket(ii).heat_transf_coef(files_rocket(ii).heat_transf_coef>1)=1;      
        files_rocket(ii).lum_ene = 0.5 .* files_rocket(ii).vel.^2 .* files_rocket(ii).mass_rate.* lum_eff ;%.* 1e10 ./ (dist.^2); lum_eff
        files_rocket(ii).lum_ene(files_rocket(ii).lum_ene==0)=nan;
        files_rocket(ii).mag = 6.8 - 1.086 * log(files_rocket(ii).lum_ene);
        matrix_holder=[];
            
        filename=append(files_rocket(ii).folder,'\',files_rocket(ii).name,'\trj.hst');
        fid = fopen( filename );
        lineRead = fgets(fid);
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\S+)','tokens'))
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    if string(cell_holder{1}(1))~='#'
                        matrix_holder=[matrix_holder; (str2double(cell_holder(1:8)))];
                    end
                end
            lineRead = fgets(fid);
        end
        jjarray=1;
        for jj=round(files_rocket(ii).sec(1)):round(files_rocket(ii).sec(end))
            files_rocket(ii).flight_ang(jjarray,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj)),6));
            files_rocket(ii).head_ang(jjarray,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj)),7));
            files_rocket(ii).lat(jjarray,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj)),3));
            files_rocket(ii).lon(jjarray,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj)),4));
            jjarray=jjarray+1;
        end
        

        %% SAVE in Rocket
        
        rocket.sec= [rocket.sec; files_rocket(ii).sec];
        rocket.time= [rocket.time; files_rocket(ii).time];
        rocket.alt= [rocket.alt; files_rocket(ii).alt];
        rocket.gr_trk= [rocket.gr_trk; files_rocket(ii).gr_trk];
        rocket.mass= [rocket.mass; files_rocket(ii).mass];
        rocket.mass_rate= [rocket.mass_rate; files_rocket(ii).mass_rate];
        rocket.vel= [rocket.vel; files_rocket(ii).vel];
        rocket.rho_air= [rocket.rho_air; files_rocket(ii).rho_air];
        rocket.heat_transf_coef= [rocket.heat_transf_coef; files_rocket(ii).heat_transf_coef];  
        rocket.lum_ene = [rocket.lum_ene; files_rocket(ii).lum_ene];
        rocket.mag = [rocket.mag; files_rocket(ii).mag];
        rocket.flight_ang= [rocket.flight_ang; files_rocket(ii).flight_ang];
        rocket.head_ang= [rocket.head_ang; files_rocket(ii).head_ang];
        rocket.lat= [rocket.lat; files_rocket(ii).lat];
        rocket.lon= [rocket.lon; files_rocket(ii).lon];
        
    else

    end
    
end

load('RocketFinalA1.mat')
sec = seconds(sara(:,1)+rocket.sec(end));
sec.Format = 'mm:ss.SSS'; 
rocket.sec= [rocket.sec; sara(:,1)+rocket.sec(end)];
rocket.time= [rocket.time; sec];
rocket.alt= [rocket.alt; sara(:,2)];
rocket.gr_trk= [rocket.gr_trk; sara(:,6)+rocket.gr_trk(end)];
rocket.mass= [rocket.mass; ones(1,length(sara(:,1)))'*rocket.mass(end)];
rocket.mass_rate= [rocket.mass_rate; ones(1,length(sara(:,1)))'*0];
rocket.vel= [rocket.vel; sara(:,5)];
[rho_airia,a_air,T_air,P_air,nu_air,h,sigma] = atmos(sara(:,2));
rocket.rho_air= [rocket.rho_air; rho_airia];
rocket.heat_transf_coef= [rocket.heat_transf_coef; ones(1,length(sara(:,1)))'*0];  
rocket.lum_ene = [rocket.lum_ene; zeros(2,length(sara(:,1)))'];
rocket.mag = [rocket.mag; zeros(2,length(sara(:,1)))'*nan];
rocket.flight_ang= [rocket.flight_ang; sara(:,12)];
rocket.head_ang= [rocket.head_ang; 90-sara(:,13)];
rocket.lat= [rocket.lat; sara(:,3)];
rocket.lon= [rocket.lon; sara(:,4)];

save('meteorRocketMass.mat','rocket')