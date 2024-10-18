clc
clear all
close all
%% DATA
run('config.m')
load('FRIPONstations.mat')

for ii=1:299
    station_name(ii)=FRIPONstations(ii).Station_Code;
end

% fileroot='C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Software\FRIPON Reader\202001-FIRST PART\';
% files= dir(fullfile([append(fileroot,'Meteors')]));
filesFRIPON=dir(fullfile('C:\Users\maxiv\Documents\TUM\4th Term Thesis\AllBertEinStein\Software\FRIPON Reader\202001-FIRST PART'));
 
matrix_holder=[];
for ii=1:length(filesFRIPON)
    if filesFRIPON(ii).isdir==1 && ~isempty(regexp(filesFRIPON(ii).name,'20','match'))        
        %% mass info
        filename=append(filesFRIPON(ii).folder,'/',filesFRIPON(ii).name,'/traj_SJ/',filesFRIPON(ii).name,'_trajectory.txt');
        fid = fopen( filename );
        lineRead = fgets(fid);
        
        iiStation=1;
        iiWhile=1;
        matrix_holder=[];
        while ischar(lineRead)
                if~isempty(regexp(lineRead,'(\d+)','match'))
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    if iiWhile==1
                        Stations=convertCharsToStrings(cell_holder{1,1});
                    end
                    matrix_holder=[matrix_holder; str2double(cell_holder(1:18))];
                    matrix_name_holder(iiWhile)=convertCharsToStrings(cell_holder{1,1});
                    if ~isequal(Stations(iiStation),convertCharsToStrings(cell_holder{1,1}))
                        iiStation=iiStation+1;
                        Stations(iiStation)=convertCharsToStrings(cell_holder{1,1});
                    end
                    iiWhile=iiWhile+1;
                end
            lineRead = fgets(fid);
        end
        
        for jj=1:length(Stations)
            filesFRIPON(ii).OBS(jj)= Stations(jj);
            filesFRIPON(ii).OBS_lon(jj)= FRIPONstations(Stations(jj)==station_name).lon;
            filesFRIPON(ii).OBS_lat(jj)= FRIPONstations(Stations(jj)==station_name).lat;
            filesFRIPON(ii).OBS_alt(jj)= FRIPONstations(Stations(jj)==station_name).alt/1000;
            
            THEsecond=matrix_holder((Stations(jj)==matrix_name_holder),4);
            minsec=min(matrix_holder(:,4));
            filesFRIPON(ii).sec{jj}= THEsecond-minsec;
            sec = seconds(filesFRIPON(ii).sec{jj});
            sec.Format = 'mm:ss.SSS'; 
            filesFRIPON(ii).time{jj}= sec;
            
            filesFRIPON(ii).lon{jj}= matrix_holder((Stations(jj)==matrix_name_holder),8);
            filesFRIPON(ii).lat{jj}= matrix_holder((Stations(jj)==matrix_name_holder),9);
            filesFRIPON(ii).alt{jj}= matrix_holder((Stations(jj)==matrix_name_holder),10);
            filesFRIPON(ii).OBS_dist{jj}= matrix_holder((Stations(jj)==matrix_name_holder),17);
            filesFRIPON(ii).mag{jj}= matrix_holder((Stations(jj)==matrix_name_holder),16);
            filesFRIPON(ii).ABSmag{jj}= matrix_holder((Stations(jj)==matrix_name_holder),18);

            x=matrix_holder((Stations(jj)==matrix_name_holder),5);
            y=matrix_holder((Stations(jj)==matrix_name_holder),6);
            z=matrix_holder((Stations(jj)==matrix_name_holder),7);
            filesFRIPON(ii).vel{jj}=[norm([x(1) y(1) z(1)]-[x(2) y(2) z(2)])/abs(filesFRIPON(ii).sec{jj}(2) - filesFRIPON(ii).sec{jj}(1)) ; vecnorm([x(1:end-1) y(1:end-1) z(1:end-1)]-[x(2:end) y(2:end) z(2:end)],2,2)./abs(filesFRIPON(ii).sec{jj}(2:end) - filesFRIPON(ii).sec{jj}(1:end-1))];
            filesFRIPON(ii).gr_trk{jj}= [0 ; vecnorm([(filesFRIPON(ii).alt{jj}(1:end-1)+r_planet).*deg2rad(filesFRIPON(ii).lat{jj}(1)-filesFRIPON(ii).lat{jj}(2:end))  (filesFRIPON(ii).alt{jj}(1:end-1)+r_planet).*deg2rad(filesFRIPON(ii).lon{jj}(1)-filesFRIPON(ii).lon{jj}(2:end))],2,2)/1000]; 
            dist=[vecnorm([x(2) y(2) z(2)]-[x(1) y(1) z(1)],2,2) ; vecnorm([x(2:end) y(2:end) z(2:end)]-[x(1:end-1) y(1:end-1) z(1:end-1)],2,2)];
            alt_diff=[filesFRIPON(ii).alt{jj}(1)-filesFRIPON(ii).alt{jj}(2); filesFRIPON(ii).alt{jj}(1:end-1)-filesFRIPON(ii).alt{jj}(2:end)];
            filesFRIPON(ii).head_ang{jj}= [rad2deg(atan2( (filesFRIPON(ii).alt{jj}(1)+r_planet).*deg2rad(filesFRIPON(ii).lat{jj}(1)-filesFRIPON(ii).lat{jj}(2)) , (filesFRIPON(ii).alt{jj}(1)+r_planet).*deg2rad(filesFRIPON(ii).lon{jj}(1)-filesFRIPON(ii).lon{jj}(2)) )) ; rad2deg(atan2( (filesFRIPON(ii).alt{jj}(1:end-1)+r_planet).*deg2rad(filesFRIPON(ii).lat{jj}(1:end-1)-filesFRIPON(ii).lat{jj}(2:end)) , (filesFRIPON(ii).alt{jj}(1:end-1)+r_planet).*deg2rad(filesFRIPON(ii).lon{jj}(1:end-1)-filesFRIPON(ii).lon{jj}(2:end)) ))];
            filesFRIPON(ii).incl_ang{jj}= rad2deg(asin( alt_diff ./ dist )); %[rad2deg(atan2( alt_diff(1) , filesFRIPON(ii).gr_trk{jj}(2) )) ; rad2deg(atan2( alt_diff(2:end) , filesFRIPON(ii).gr_trk{jj}(2:end) ))];
        end
        
        filename=append(filesFRIPON(ii).folder,'/',filesFRIPON(ii).name,'/dynamique_test/',filesFRIPON(ii).name,'_fripouille_c_report.txt');
        fid = fopen( filename );
        lineRead = fgets(fid);
        matrix_holder=[];
        FinalPercMass=[];
        while ischar(lineRead)
            
            if ~isempty(regexp('Final Mass : 7.65983645492e-146 %','Final Mass : ','match'))
                FinalPerc=(regexp('Final Mass : 7.65983645492e-146 %','[0-9.e-]','match'));
                FinalPercMass=str2double(append(FinalPerc{:}));
            end
            
            if~isempty(regexp(lineRead,'rho','match'))
                for iiDens=1:6
                    lineRead = fgets(fid);
                    cell_holder=regexp(lineRead,'(\S+)','match');
                    matrix_holder=[matrix_holder; (str2double(cell_holder(1:3)))];
                    
                end
                
                break
                
            end
            lineRead = fgets(fid);
        end
        
        filesFRIPON(ii).rho=matrix_holder(:,1);
        filesFRIPON(ii).mass=matrix_holder(:,3);
        filesFRIPON(ii).size=matrix_holder(:,2);
%         filesFRIPON(ii).mass_fin=matrix_holder(:,3).*FinalPercMass/100;

    else
    end
end

% % % save('FRIPONmeteors.mat','filesFRIPON')