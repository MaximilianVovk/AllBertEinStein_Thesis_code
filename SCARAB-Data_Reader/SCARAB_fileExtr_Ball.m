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

% files_rocket = dir(fullfile([append(fileroot,'B_4\b-4-10cm-130km7.97kms-1.83deg')]));

lum_eff = [0.1 0.001];

ball.material="Iron";
ball.rho=7870;
ball.Hmelt=272000;
ball.init_size=10;
ball.init_alt=110000;
ball.init_vel=7.99;
ball.init_ang=1.73;
ball.init_inc=97.5;
ball.crosSec=[];
ball.sec= [];
ball.time= [];
ball.alt= [];
ball.gr_trk= [];
ball.mass= [];
ball.mass_rate= [];
ball.vel= [];
ball.rho_air= [];
ball.heat_transf_coef= [];    
ball.lum_ene = [];
ball.mag = [];
ball.flight_ang= [];
ball.head_ang= [];
ball.lat= [];
ball.lon=[];


matrix_holder=[];
for ii=1:length(files_rocket)
    if files_rocket(ii).isdir==1 && ~isempty(regexp(files_rocket(ii).name,'.1','match')) && isempty(regexp(files_rocket(ii).name,'3.1','match'))       

        files_rocket(ii).material=ball.material;
        files_rocket(ii).rho=ball.rho;
        files_rocket(ii).Hmelt=ball.Hmelt;
        files_rocket(ii).init_size=ball.init_size;
        files_rocket(ii).init_alt=ball.init_alt;
        files_rocket(ii).init_vel=ball.init_vel;
        files_rocket(ii).init_ang=ball.init_ang;
        files_rocket(ii).init_inc=ball.init_inc;
        
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
        files_rocket(ii).crosSec=pi*(files_rocket(ii).mass*3/4 / (pi*files_rocket(ii).rho)).^(2/3);
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
        ball.crosSec=[ball.crosSec; files_rocket(ii).crosSec];
        ball.sec= [ball.sec; files_rocket(ii).sec];
        ball.time= [ball.time; files_rocket(ii).time];
        ball.alt= [ball.alt; files_rocket(ii).alt];
        ball.gr_trk= [ball.gr_trk; files_rocket(ii).gr_trk];
        ball.mass= [ball.mass; files_rocket(ii).mass];
        ball.mass_rate= [ball.mass_rate; files_rocket(ii).mass_rate];
        ball.vel= [ball.vel; files_rocket(ii).vel];
        ball.rho_air= [ball.rho_air; files_rocket(ii).rho_air];
        ball.heat_transf_coef= [ball.heat_transf_coef; files_rocket(ii).heat_transf_coef];  
        ball.lum_ene = [ball.lum_ene; files_rocket(ii).lum_ene];
        ball.mag = [ball.mag; files_rocket(ii).mag];
        ball.flight_ang= [ball.flight_ang; files_rocket(ii).flight_ang];
        ball.head_ang= [ball.head_ang; files_rocket(ii).head_ang];
        ball.lat= [ball.lat; files_rocket(ii).lat];
        ball.lon= [ball.lon; files_rocket(ii).lon];
        
    else

    end
    
end

save('meteorBallMass.mat','ball')