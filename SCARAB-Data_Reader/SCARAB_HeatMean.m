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

ironcoef20cm_pos=[];
bmaxcoef20cm_pos=[];
bmincoef20cm_pos=[];
cchocoef20cm_pos=[];
ochocoef20cm_pos=[];
grancoef20cm_pos=[];
sandcoef20cm_pos=[];

ironcoef10cm_pos=[];
bmaxcoef10cm_pos=[];
bmincoef10cm_pos=[];
cchocoef10cm_pos=[];
ochocoef10cm_pos=[];
grancoef10cm_pos=[];
sandcoef10cm_pos=[];

ironcoef5cm_pos=[];
bmaxcoef5cm_pos=[];
bmincoef5cm_pos=[];
cchocoef5cm_pos=[];
ochocoef5cm_pos=[];
grancoef5cm_pos=[];
sandcoef5cm_pos=[];


matrix_holder=[];
load('meteorMass.mat')
for ii=1:length(files)       
        %% sorting
        
        if files(ii).init_size==5
%               figure(1)
                switch files(ii).material
                    case "iron"
%                         subplot(2,4,1)
%                         title('Iron')
                        ironcoef5cm_pos(end+1,1)=ii;
                        ironcoef5cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "basaltmax"
%                         subplot(2,4,2)
%                         title('Basalt Acidic')
                        holder_split=regexp(files(ii).name,'m-basaltmax5cm-','split');
                        bmax5cm=[bmax5cm convertCharsToStrings(holder_split{2})];
                        bmaxcoef5cm_pos(end+1,1)=ii;
                        bmaxcoef5cm_pos(end,2)=length(files(ii).heat_transf_coef);
                        
                    case "minbasalt"
%                         subplot(2,4,3)
%                         title('Basalt Basic')
                        bmincoef5cm_pos(end+1,1)=ii;
                        bmincoef5cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "ordchon"
%                         subplot(2,4,4)
%                         title('Ordinary Chondrite')
                        ochocoef5cm_pos(end+1,1)=ii;
                        ochocoef5cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "carchon"
%                         subplot(2,4,5)
%                         title('Carbonaceous Chondrite')
                        cchocoef5cm_pos(end+1,1)=ii;
                        cchocoef5cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "granite"
%                         subplot(2,4,6)
%                         title('Granite')
                        grancoef5cm_pos(end+1,1)=ii;
                        grancoef5cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "sandstone"
%                         subplot(2,4,7)
%                         title('Sandstone')
                        sandcoef5cm_pos(end+1,1)=ii;
                        sandcoef5cm_pos(end,2)=length(files(ii).heat_transf_coef);
                end
%                 plot(files(ii).heat_transf_coef,files(ii).alt,material_color(material_name==files(ii).material))
%                 ylim([0 100])
%                 xlabel('Heat transfer coeficient') 
%                 ylabel('Altitude [km]')
%                 grid on
%                 hold on
%                 if files(ii).mass(end)==0
%                     plot(files(ii).heat_transf_coef(end-1),files(ii).alt(end-1),append("x",material_color(material_name==files(ii).material)),'Linewidth',2)
%                 end
% 
%               sgtitle('Meteor of 5 cm')
% %               subplot(2,4,8)
% %               plot(linspace(0,1,100),nan)
% %                 legend(bmax5cm,'location','westoutside');
% %                 legend(bmax5cm,'location','west');
% %                 axis off
                
        elseif files(ii).init_size==10
%               figure(2)
                switch files(ii).material
                    case "iron"
%                         subplot(2,4,1)
%                         title('Iron')
                        ironcoef10cm_pos(end+1,1)=ii;
                        ironcoef10cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "basaltmax"
%                         subplot(2,4,2)
%                         title('Basalt Acidic')
                        holder_split=regexp(files(ii).name,'m-basaltmax10cm-','split');
                        bmax10cm=[bmax10cm convertCharsToStrings(holder_split{2})];
                        bmaxcoef10cm_pos(end+1,1)=ii;
                        bmaxcoef10cm_pos(end,2)=length(files(ii).heat_transf_coef);
                        
                    case "minbasalt"
%                         subplot(2,4,3)
%                         title('Basalt Basic')
                        bmincoef10cm_pos(end+1,1)=ii;
                        bmincoef10cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "ordchon"
%                         subplot(2,4,4)
%                         title('Ordinary Chondrite')
                        ochocoef10cm_pos(end+1,1)=ii;
                        ochocoef10cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "carchon"
%                         subplot(2,4,5)
%                         title('Carbonaceous Chondrite')
                        cchocoef10cm_pos(end+1,1)=ii;
                        cchocoef10cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "granite"
%                         subplot(2,4,6)
%                         title('Granite')
                        grancoef10cm_pos(end+1,1)=ii;
                        grancoef10cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "sandstone"
%                         subplot(2,4,7)
%                         title('Sandstone')
                        sandcoef10cm_pos(end+1,1)=ii;
                        sandcoef10cm_pos(end,2)=length(files(ii).heat_transf_coef);
                end
%                 plot(files(ii).heat_transf_coef,files(ii).alt,material_color(material_name==files(ii).material))
%                 ylim([0 100])
%                 xlabel('Heat transfer coeficient') 
%                 ylabel('Altitude [km]')
%                 grid on
%                 hold on
%                 if files(ii).mass(end)==0
%                     plot(files(ii).heat_transf_coef(end-1),files(ii).alt(end-1),append("x",material_color(material_name==files(ii).material)),'Linewidth',2)
%                 end
%                 
%               sgtitle('Meteor of 10 cm')
% %               subplot(2,4,8)
% %               plot(linspace(0,1,100),nan)
% %                 legend(bmax10cm,'location','westoutside');
% %                 legend(bmax10cm,'location','west');
% %                 axis off
                
        elseif files(ii).init_size==20
%               figure(3)
                switch files(ii).material
                    case "iron"
%                         subplot(2,4,1)
%                         title('Iron')
                        ironcoef20cm_pos(end+1,1)=ii;
                        ironcoef20cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "basaltmax"
%                         subplot(2,4,2)
%                         title('Basalt Acidic')
                        holder_split=regexp(files(ii).name,'m-basaltmax20cm-','split');
                        bmax20cm=[bmax20cm convertCharsToStrings(holder_split{2})];
                        bmaxcoef20cm_pos(end+1,1)=ii;
                        bmaxcoef20cm_pos(end,2)=length(files(ii).heat_transf_coef);
                        
                    case "minbasalt"
%                         subplot(2,4,3)
%                         title('Basalt Basic')
                        bmincoef20cm_pos(end+1,1)=ii;
                        bmincoef20cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "ordchon"
%                         subplot(2,4,4)
%                         title('Ordinary Chondrite')
                        ochocoef20cm_pos(end+1,1)=ii;
                        ochocoef20cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "carchon"
%                         subplot(2,4,5)
%                         title('Carbonaceous Chondrite')
                        cchocoef20cm_pos(end+1,1)=ii;
                        cchocoef20cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "granite"
%                         subplot(2,4,6)
%                         title('Granite')
                        grancoef20cm_pos(end+1,1)=ii;
                        grancoef20cm_pos(end,2)=length(files(ii).heat_transf_coef);
                    case "sandstone"
%                         subplot(2,4,7)
%                         title('Sandstone')
                        sandcoef20cm_pos(end+1,1)=ii;
                        sandcoef20cm_pos(end,2)=length(files(ii).heat_transf_coef);
                end
%                 plot(files(ii).heat_transf_coef,files(ii).alt,material_color(material_name==files(ii).material))
%                 ylim([0 100])
%                 xlabel('Heat transfer coeficient') 
%                 ylabel('Altitude [km]')
%                 grid on
%                 hold on
%                 if files(ii).mass(end)==0
%                     plot(files(ii).heat_transf_coef(end-1),files(ii).alt(end-1),append("x",material_color(material_name==files(ii).material)),'Linewidth',2)
%                 end
%                 
%               sgtitle('Meteor of 20 cm')
% %               subplot(2,4,8)
% %               plot(linspace(0,1,100),nan)
% %                 legend(bmax20cm,'location','westoutside');
% %                 legend(bmax20cm,'location','west');
% %                 axis off
                
        end
        
%     else
% 
%     end
    
end

allsize_pos{1}=ironcoef20cm_pos;
allsize_pos{2}=bmaxcoef20cm_pos;
allsize_pos{3}=bmincoef20cm_pos;
allsize_pos{4}=cchocoef20cm_pos;
allsize_pos{5}=ochocoef20cm_pos;
allsize_pos{6}=grancoef20cm_pos;
allsize_pos{7}=sandcoef20cm_pos;

allsize_pos{8}=ironcoef10cm_pos;
allsize_pos{9}=bmaxcoef10cm_pos;
allsize_pos{10}=bmincoef10cm_pos;
allsize_pos{11}=cchocoef10cm_pos;
allsize_pos{12}=ochocoef10cm_pos;
allsize_pos{13}=grancoef10cm_pos;
allsize_pos{14}=sandcoef10cm_pos;

allsize_pos{15}=ironcoef5cm_pos;
allsize_pos{16}=bmaxcoef5cm_pos;
allsize_pos{17}=bmincoef5cm_pos;
allsize_pos{18}=cchocoef5cm_pos;
allsize_pos{19}=ochocoef5cm_pos;
allsize_pos{20}=grancoef5cm_pos;
allsize_pos{21}=sandcoef5cm_pos;

coefallsize=zeros(100,size(ironcoef20cm_pos,1));
for qq=1:size(allsize_pos,2)
    for ii=1:size(allsize_pos{qq},1)
        for jj=1:100
            coefallsize(jj,ii)=mean(files(allsize_pos{qq}(ii,1)).heat_transf_coef(find(round(files(allsize_pos{qq}(ii,1)).alt)==jj)));
            if isnan(coefallsize(jj,ii))
                coefallsize(jj,ii)=0;
                
            end
        end
    end
    
    files(allsize_pos{qq}(end,1)).coef_avg=(mean(coefallsize,2));%smooth

end

save('meteorMass.mat','files')
% plot(files(allsize_pos{qq}(end,1)).coef_avg,1:100,'Linewidth',1.2)
% hold on
% grid on
% plot(smooth(mean(coefallsize,2)),1:100,'Linewidth',1.2)
% plot(smooth(mean(coefallsize,2),10),1:100,'Linewidth',1.2)
% plot(smooth(mean(coefallsize,2),12),1:100,'Linewidth',1.2)
% legend('mean','smooth','smooth 10','smooth 20')
% run('SCARAB_Subplots.m')
% run('SCARAB_3Dplots.m')