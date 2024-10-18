clc
clear all
close all

%%
fileroot='C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\SCARAB\';
% files= dir(fullfile([append(fileroot,'Meteors')]));
% files_rocket= dir(fullfile([append(fileroot,'B_4\b-4-10cm-110km7.99kms-1.73deg')]));
files_rocket= dir(fullfile([append(fileroot,'B_4\b-4-10cm-130km7.97kms-1.83deg')]));
% files_rocket= dir(fullfile([append(fileroot,'B_4\b-4-5cm-130km7.97kms-1.83deg')]));
%files=dir(fullfile('C:\Users\maxivy\Documents\DRAMA\TEST folder','*.txt'));

lum_eff = [0.1 0.001]; 

cage.material="Aluminum";
cage.rho=2813;
cage.Hmelt=400000; %titanium
cage.init_size=20;
cage.init_alt=110000;
cage.init_vel=7.99;
cage.init_ang=1.73;
cage.init_inc=97.5;
cage.crosSec=[];
cage.sec= [];
cage.time= [];
cage.alt= [];
cage.gr_trk= [];
cage.mass= [];
cage.mass_rate= [];
cage.vel= [];
cage.rho_air= [];
cage.heat_transf_coef= [];    
cage.lum_ene = [];
cage.mag = [];
cage.flight_ang= [];
cage.head_ang= [];
cage.lat= [];
cage.lon=[];


matrix_holder=[];
for ii=1:length(files_rocket)
    if files_rocket(ii).isdir==1 && ~isempty(regexp(files_rocket(ii).name,'.2','match'))        

        files_rocket(ii).material=cage.material;
        files_rocket(ii).rho=cage.rho;
        files_rocket(ii).Hmelt=cage.Hmelt;
        files_rocket(ii).init_size=cage.init_size;
        files_rocket(ii).init_alt=cage.init_alt;
        files_rocket(ii).init_vel=cage.init_vel;
        files_rocket(ii).init_ang=cage.init_ang;
        files_rocket(ii).init_inc=cage.init_inc;
        
        %% mass info
        filename=append(files_rocket(ii).folder,'/',files_rocket(ii).name,'/mas.hst');
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
        files_rocket(ii).crosSec=pi*(files_rocket(ii).mass ./ (pi*0.005*files_rocket(ii).rho));
        files_rocket(ii).crosSec(files_rocket(ii).crosSec==0)=min(files_rocket(ii).crosSec(files_rocket(ii).crosSec>0));
        
        
        filename=append(files_rocket(ii).folder,'/',files_rocket(ii).name,'/atm.hst');
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
        cage.crosSec=[cage.crosSec; files_rocket(ii).crosSec];
        cage.sec= [cage.sec; files_rocket(ii).sec];
        cage.time= [cage.time; files_rocket(ii).time];
        cage.alt= [cage.alt; files_rocket(ii).alt];
        cage.gr_trk= [cage.gr_trk; files_rocket(ii).gr_trk];
        cage.mass= [cage.mass; files_rocket(ii).mass];
        cage.mass_rate= [cage.mass_rate; files_rocket(ii).mass_rate];
        cage.vel= [cage.vel; files_rocket(ii).vel];
        cage.rho_air= [cage.rho_air; files_rocket(ii).rho_air];
        cage.heat_transf_coef= [cage.heat_transf_coef; files_rocket(ii).heat_transf_coef];  
        cage.lum_ene = [cage.lum_ene; files_rocket(ii).lum_ene];
        cage.mag = [cage.mag; files_rocket(ii).mag];
        cage.flight_ang= [cage.flight_ang; files_rocket(ii).flight_ang];
        cage.head_ang= [cage.head_ang; files_rocket(ii).head_ang];
        cage.lat= [cage.lat; files_rocket(ii).lat];
        cage.lon= [cage.lon; files_rocket(ii).lon];
        
    else

    end
    
end

save('meteorcageMass.mat','cage')