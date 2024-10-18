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

matrix_holder=[];
load('meteorMass.mat')
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))               
        %% subplot
        
        if files(ii).init_size==5
              figure(1)
                switch files(ii).material
                    case "iron"
                        subplot(2,4,1)
                        title('Iron')
                    case "basaltmax"
                        subplot(2,4,2)
                        title('Basalt Acidic')
                        holder_split=regexp(files(ii).name,'m-basaltmax5cm-','split');
                        bmax5cm=[bmax5cm convertCharsToStrings(holder_split{2})];
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
                plot(files(ii).heat_transf_coef(files(ii).heat_transf_coef<1),files(ii).alt(files(ii).heat_transf_coef<1))
                ylim([0 100])
                xlabel('Heat transfer coef.') 
                ylabel('Altitude [km]')
                grid on
                hold on
                xlim([0 0.4])
%                 if ~isempty(files(ii).coef_avg)
%                     plot(files(ii).coef_avg,1:100,'b.','Linewidth',1.2)
%                     plot(files(ii).coef_gauss_avg,1:100,'r-','Linewidth',1.2)       
%                 end
%                 if ~isempty(files(ii).coef_avg)
%                     plot(files(ii).coef_avg,1:100,material_color(material_name==files(ii).material),'Linewidth',1.2)
%                     plot(files(ii).coef_gauss_avg,1:100,append(material_color(material_name==files(ii).material),'--'),'Linewidth',1.2)       
%                 end
                
                sgtitle('Meteor of 5 cm')
                subplot(2,4,8)
                plot(linspace(0,1,100),nan)
                legend(bmax5cm,'location','west');
                axis off
                
        elseif files(ii).init_size==10
              figure(2)
                switch files(ii).material
                    case "iron"
                        subplot(2,4,1)
                        title('Iron')
                    case "basaltmax"
                        subplot(2,4,2)
                        title('Basalt Acidic')
                        holder_split=regexp(files(ii).name,'m-basaltmax10cm-','split');
                        bmax10cm=[bmax10cm convertCharsToStrings(holder_split{2})];
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
                plot(files(ii).heat_transf_coef(files(ii).heat_transf_coef<1),files(ii).alt(files(ii).heat_transf_coef<1))
                ylim([0 100])
                xlabel('Heat transfer coef.') 
                ylabel('Altitude [km]')
                grid on
                hold on
                xlim([0 0.25])
%                 if ~isempty(files(ii).coef_avg)
%                     plot(files(ii).coef_avg,1:100,'b.','Linewidth',1.2)
%                     plot(files(ii).coef_gauss_avg,1:100,'r-','Linewidth',1.2)       
%                 end
%                 if ~isempty(files(ii).coef_avg)
%                     plot(files(ii).coef_avg,1:100,material_color(material_name==files(ii).material),'Linewidth',1.2)
%                     plot(files(ii).coef_gauss_avg,1:100,append(material_color(material_name==files(ii).material),'--'),'Linewidth',1.2)                
%                 end
                
                sgtitle('Meteor of 10 cm')
                subplot(2,4,8)
                plot(linspace(0,1,100),nan)
                legend(bmax10cm,'location','west');
                axis off
                
        elseif files(ii).init_size==20
              figure(3)
                switch files(ii).material
                    case "iron"
                        subplot(2,4,1)
                        title('Iron')
                    case "basaltmax"
                        subplot(2,4,2)
                        title('Basalt Acidic')
                        holder_split=regexp(files(ii).name,'m-basaltmax20cm-','split');
                        bmax20cm=[bmax20cm convertCharsToStrings(holder_split{2})];
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
                plot(files(ii).heat_transf_coef(files(ii).heat_transf_coef<1),files(ii).alt(files(ii).heat_transf_coef<1))
                ylim([0 100])
%                 xlim([0 1])
                xlabel('Heat transfer coef.') 
                ylabel('Altitude [km]')
                grid on
                hold on
                xlim([0 0.2])
%                 if ~isempty(files(ii).coef_avg)
%                     plot(files(ii).coef_avg,1:100,'b.','Linewidth',1.2)
%                     plot(files(ii).coef_gauss_avg,1:100,'r-','Linewidth',1.2)       
%                 end
                
                sgtitle('Meteor of 20 cm')
                subplot(2,4,8)
                plot(linspace(0,1,100),nan)
                legend(bmax20cm,'location','west');
                axis off
        end
        
    else

    end
    
end

