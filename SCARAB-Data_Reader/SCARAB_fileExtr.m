clc
clear all
close all

%%
fileroot='C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\SCARAB\';
% files= dir(fullfile([append(fileroot,'Meteors')]));
files= dir(fullfile([append(fileroot,'Meteors Iron 20cm')]));
%files=dir(fullfile('C:\Users\maxivy\Documents\DRAMA\TEST folder','*.txt'));


material_name=["iron" "basaltmax" "minbasalt" "ordchon" "carchon" "granite" "sandstone"];
material_density=[7870.0 2400.0 3100.0 3500.0 2800.0 2750.0 2000.0];
material_meltingheat=[272000.0 400000.0 506000.0 265000.0 265000.0 250000.0 680000.0];
material_color=["g" "b" "c" "r" "m" "k" "y"];
lum_eff = [0.1 0.001];

qq=1;
matrix_holder=[];
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match')) %&& ~isempty(regexp(files(ii).name,'1.83','match'))      
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
            
        filename=append(files(ii).folder,'/',files(ii).name,'/0.1/trj.hst');
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
        
        for jj=1:length(files(ii).sec)
            files(ii).flight_ang(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),6));
            files(ii).head_ang(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),7));
            files(ii).lat(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),3));
            files(ii).lon(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),4));
        end
        
        files(ii).flight_ang(1,1)=-files(ii).init_ang;
% % %         files(ii).head_ang(1,1)=;
        files(ii).lat(1,1)=0;
        files(ii).lon(1,1)=0;

        
%         files(ii).flight_ang= matrix_holder(:,6);
%         files(ii).head_ang= matrix_holder(:,7);   
%         files(ii).lat= matrix_holder(:,3);
%         files(ii).lon= matrix_holder(:,4); 

        filename=append(files(ii).folder,'/',files(ii).name,'/0.1/aero.hst');
        fid = fopen( filename );
        lineRead = fgets(fid);
        matrix_holder=[];
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\S+)','match'))
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    if string(cell_holder{1}(1))~='#'
                        matrix_holder=[matrix_holder; (str2double(cell_holder(1:7)))];
                    end
                end
            lineRead = fgets(fid);
        end
        for jj=1:length(files(ii).sec)
            files(ii).Ma(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),6));
            files(ii).Kn(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),7));
        end

        filename=append(files(ii).folder,'/',files(ii).name,'/0.1/att.hst');
        fid = fopen( filename );
        lineRead = fgets(fid);
        matrix_holder=[];
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\S+)','match'))
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    if string(cell_holder{1}(1))~='#'
                        matrix_holder=[matrix_holder; (str2double(cell_holder(1:9)))];
                    end
                end
            lineRead = fgets(fid);
        end
        for jj=1:length(files(ii).sec)
            files(ii).ang_phiX(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),4));
            files(ii).ang_thetaY(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),5));
            files(ii).ang_psiZ(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),6));

            files(ii).ang_velX(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),7));
            files(ii).ang_velY(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),8));
            files(ii).ang_velZ(jj,1)=mean(matrix_holder((find(round(matrix_holder(:,1))==jj-1)),9));
        end
        
        matrix_holder=[];
        
        
        
% files_mission(qq).material=files(ii).material;
% files_mission(qq).rho=files(ii).rho;
% files_mission(qq).Hmelt=files(ii).Hmelt;
% files_mission(qq).init_size=files(ii).init_size;
% files_mission(qq).init_alt=files(ii).init_alt;
% files_mission(qq).init_vel=files(ii).init_vel;
% files_mission(qq).init_ang=files(ii).init_ang;
% files_mission(qq).init_inc=files(ii).init_inc;
% files_mission(qq).crosSec=files(ii).crosSec;
% files_mission(qq).sec= files(ii).sec;
% files_mission(qq).time= files(ii).time;
% files_mission(qq).alt= files(ii).alt;
% files_mission(qq).gr_trk= files(ii).gr_trk;
% files_mission(qq).mass= files(ii).mass;
% files_mission(qq).mass_rate= files(ii).mass_rate;
% files_mission(qq).vel= files(ii).vel;
% files_mission(qq).rho_air= files(ii).rho_air;
% files_mission(qq).heat_transf_coef= files(ii).heat_transf_coef;   
% files_mission(qq).lum_ene = files(ii).lum_ene;
% files_mission(qq).mag = files(ii).mag;
% files_mission(qq).flight_ang= files(ii).flight_ang;
% files_mission(qq).head_ang= files(ii).head_ang;
% files_mission(qq).lat= files(ii).lat;
% files_mission(qq).lon= files(ii).lon;
%         
%         qq=qq+1;
    else

    end
    
end

save('SACARBIronKnMa.mat','files')
% save('meteorMissionMass.mat','files_mission')