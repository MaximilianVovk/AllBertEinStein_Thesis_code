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

bmax10cm=[];

matrix_holder=[];
load('meteorMass.mat')
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))        
        %% subplot
        
        if files(ii).init_size==5
            size5=figure(1);
            size5.WindowState = 'maximized';
        elseif files(ii).init_size==10
            size10=figure(2);
            size10.WindowState = 'maximized';
        elseif files(ii).init_size==20
            size20=figure(3);
            size20.WindowState = 'maximized';
        end
        
        switch files(ii).material
            case "iron"
                subplot(3,3,1)
                title('Iron')
            case "basaltmax"
                subplot(3,3,2)
                title('Basalt Acidic')
                if files(ii).init_size==10
                    holder_split=regexp(files(ii).name,'m-basaltmax10cm-','split');
                    bmax10cm=[bmax10cm convertCharsToStrings(holder_split{2})];
                end
            case "minbasalt"
                subplot(3,3,3)
                title('Basalt Basic')
            case "ordchon"
                subplot(3,3,4)
                title('Ordinary Chondrite')
            case "carchon"
                subplot(3,3,5)
                title('Carbonaceous Chondrite')
            case "granite"
                subplot(3,3,6)
                title('Granite')
            case "sandstone"
                subplot(3,3,7)
                title('Sandstone')
        end
        plot(files(ii).heat_transf_coef100km,1:100) %heat_transf_coef100km
        xlabel('Heat transfer coef.') 
        ylabel('Altitude / km')
        xlim([0 0.4])
        grid on
        hold on
        if files(ii).init_size==5
            sgtitle('Meteor of 5 cm')
        elseif files(ii).init_size==10
            sgtitle('Meteor of 10 cm')
        elseif files(ii).init_size==20
            sgtitle('Meteor of 20 cm')
        end
                
        subplot(3,3,9)
        plot(linspace(0,1,100),nan)
        legend(bmax10cm,'location','west');
        axis off
    else

    end
    
end

% size5=figure(1);
% subplot(2,4,1)
% plot(0,0)
% subplot(2,4,2)
% plot(0,0)
% subplot(2,4,3)
% plot(0,0)
% subplot(2,4,4)
% plot(0,0)
% subplot(2,4,5)
% plot(0,0)
% subplot(2,4,6)
% plot(0,0)
% subplot(2,4,7)
% plot(0,0)
% 
% size10=figure(2);
% subplot(2,4,1)
% plot(0,0)
% subplot(2,4,2)
% plot(0,0)
% subplot(2,4,3)
% plot(0,0)
% subplot(2,4,4)
% plot(0,0)
% subplot(2,4,5)
% plot(0,0)
% subplot(2,4,6)
% plot(0,0)
% subplot(2,4,7)
% plot(0,0)
% 
% size20=figure(3);
% subplot(2,4,1)
% plot(0,0)
% subplot(2,4,2)
% plot(0,0)
% subplot(2,4,3)
% plot(0,0)
% subplot(2,4,4)
% plot(0,0)
% subplot(2,4,5)
% plot(0,0)
% subplot(2,4,6)
% plot(0,0)
% subplot(2,4,7)
% plot(0,0)
% 
% 
% for ii=1:length(files)
%     if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))
%         if files(ii).init_size==5
%             size5=figure(1);
%             size5.WindowState = 'maximized';
%         elseif files(ii).init_size==10
%             size10=figure(2);
%             size10.WindowState = 'maximized';
%         elseif files(ii).init_size==20
%             size20=figure(3);
%             size20.WindowState = 'maximized';
%         end
%         
%         switch files(ii).material
%             case "iron"
%                 subplot(2,4,1)
%                 title('Iron Gauss')
%             case "basaltmax"
%                 subplot(2,4,2)
%                 title('Basalt Acidic Gauss')
%             case "minbasalt"
%                 subplot(2,4,3)
%                 title('Basalt Basic Gauss')
%             case "ordchon"
%                 subplot(2,4,4)
%                 title('Ordinary Chondrite Gauss')
%             case "carchon"
%                 subplot(2,4,5)
%                 title('Carbonaceous Chondrite Gauss')
%             case "granite"
%                 subplot(2,4,6)
%                 title('Granite Gauss')
%             case "sandstone"
%                 subplot(2,4,7)
%                 title('Sandstone Gauss')
%         end
%         plot(files(ii).coef_gauss,1:100,'--','Linewidth',1.2)
%         xlabel('Heat transfer coef.') 
%         ylabel('Altitude [km]')
%         grid on
%         hold on
%         if files(ii).init_size==5
%             sgtitle('Meteor of 5 cm')
%         elseif files(ii).init_size==10
%             sgtitle('Meteor of 10 cm')
%         elseif files(ii).init_size==20
%             sgtitle('Meteor of 20 cm')
%         end
%     else
% 
%     end
%     
% end