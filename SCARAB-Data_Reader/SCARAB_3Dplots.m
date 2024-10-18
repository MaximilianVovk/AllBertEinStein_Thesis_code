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
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))               
%% subplot
        if ~isempty(files(ii).coef_avg)
            switch files(ii).material
                case "iron"
                    subplot(2,4,1)
                    title('Iron')
                case "basaltmax"
                    subplot(2,4,2)
                    title('Basalt Acidic')
                case "minbasalt"
                    subplot(2,4,3)
                    title('Basalt Basic')
                case "ordchon"
                    subplot(2,4,4)
                    title('Ordinary Chondrite')
                case "carchon"
                    subplot(2,4,5)
                    title('Carbonaceous Chondrite')
                case "granite"
                    subplot(2,4,6)
                    title('Granite')
                case "sandstone"
                    subplot(2,4,7)
                    title('Sandstone')
            end
            if files(ii).init_size==5
                Z(:,1)=files(ii).coef_gauss;
                surf(X,Y,Z)
                Z=zeros(100,3);
                zlabel('Heat transfer coef.') 
                xlabel('Size [cm]')
                ylabel('Altitude [km]')
            elseif files(ii).init_size==10
                Z(:,2)=files(ii).coef_gauss;
            elseif files(ii).init_size==20
                Z(:,3)=files(ii).coef_gauss;
            end
        end
    end
end