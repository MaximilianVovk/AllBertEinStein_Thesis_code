clc
clear all
close all
%%
load('SACARBMissionKnMa.mat')
iron5cm=[];basaltmax5cm=[];minbasalt5cm=[];ordchon5cm=[];carchon5cm=[];granite5cm=[];sandstone5cm=[];
iron10cm=[];basaltmax10cm=[];minbasalt10cm=[];ordchon10cm=[];carchon10cm=[];granite10cm=[];sandstone10cm=[];
iron20cm=[];basaltmax20cm=[];minbasalt20cm=[];ordchon20cm=[];carchon20cm=[];granite20cm=[];sandstone20cm=[];
for ii=1:length(files)
    mean_heat = mean(files(ii).heat_transf_coef(files(ii).heat_transf_coef>0));
        if files(ii).init_size==5
            switch files(ii).material
            case "iron"
                iron5cm=[iron5cm  append(num2str(mean_heat), ' & ')];
            case "basaltmax"
                basaltmax5cm=[basaltmax5cm  append(num2str(mean_heat), ' & ')];
            case "minbasalt"
                minbasalt5cm=[minbasalt5cm  append(num2str(mean_heat), ' & ')];
            case "ordchon"
                ordchon5cm=[ordchon5cm  append(num2str(mean_heat), ' & ')];
            case "carchon"
                carchon5cm=[carchon5cm  append(num2str(mean_heat), ' & ')];
            case "granite"
                granite5cm=[granite5cm  append(num2str(mean_heat), ' & ')];
            case "sandstone"
                sandstone5cm=[sandstone5cm  append(num2str(mean_heat), ' & ')];
            end
        elseif files(ii).init_size==10
            switch files(ii).material
            case "iron"
                iron10cm=[iron10cm  append(num2str(mean_heat), ' & ')];
            case "basaltmax"
                basaltmax10cm=[basaltmax10cm  append(num2str(mean_heat), ' & ')];
            case "minbasalt"
                minbasalt10cm=[minbasalt10cm  append(num2str(mean_heat), ' & ')];
            case "ordchon"
                ordchon10cm=[ordchon10cm  append(num2str(mean_heat), ' & ')];
            case "carchon"
                carchon10cm=[carchon10cm  append(num2str(mean_heat), ' & ')];
            case "granite"
                granite10cm=[granite10cm  append(num2str(mean_heat), ' & ')];
            case "sandstone"
                sandstone10cm=[sandstone10cm  append(num2str(mean_heat), ' & ')];
            end
        elseif files(ii).init_size==20
            switch files(ii).material
            case "iron"
                iron20cm=[iron20cm  append(num2str(mean_heat), ' & ')];
            case "basaltmax"
                basaltmax20cm=[basaltmax20cm  append(num2str(mean_heat), ' & ')];
            case "minbasalt"
                minbasalt20cm=[minbasalt20cm  append(num2str(mean_heat), ' & ')];
            case "ordchon"
                ordchon20cm=[ordchon20cm  append(num2str(mean_heat), ' & ')];
            case "carchon"
                carchon20cm=[carchon20cm  append(num2str(mean_heat), ' & ')];
            case "granite"
                granite20cm=[granite20cm  append(num2str(mean_heat), ' & ')];
            case "sandstone"
                sandstone20cm=[sandstone20cm  append(num2str(mean_heat), ' & ')];
            end
        end
%     xticklabels({'Iron','Acidic Basalt','Basic Basalt','Ordinary Chondrite','Carbonaceous Chondrite','Granite','Sandstone'})
%     yticklabels({'7.6 km/s 2° LEO','7.97 km/s 1.83° SSO','8.5 km/s 5° SSO','8 km/s 5° SSO','9 km/s 10° SSO','9 km/s 5° SSO'})

end

% T5cm = table(iron5cm',basaltmax5cm',minbasalt5cm',ordchon5cm',carchon5cm',granite5cm',sandstone5cm')
% T10cm = table(iron10cm',basaltmax10cm',minbasalt10cm',ordchon10cm',carchon10cm',granite10cm',sandstone10cm')
% T20cm = table(iron20cm',basaltmax20cm',minbasalt20cm',ordchon20cm',carchon20cm',granite20cm',sandstone20cm')

sprintf('Iron & %s \\\\ \n Acidic Basalt & %s \\\\ \n Basic Basalt & %s \\\\ \n Carbonaceous Chondrite & %s \\\\ \n Ordinary Chondrite & %s \\\\ \n Granite & %s \\\\ \n Sandstone & %s \\\\ \n',iron5cm ,basaltmax5cm ,minbasalt5cm ,ordchon5cm ,carchon5cm ,granite5cm ,sandstone5cm)
sprintf('Iron & %s \\\\ \n Acidic Basalt & %s \\\\ \n Basic Basalt & %s \\\\ \n Carbonaceous Chondrite & %s \\\\ \n Ordinary Chondrite & %s \\\\ \n Granite & %s \\\\ \n Sandstone & %s \\\\ \n',iron10cm ,basaltmax10cm ,minbasalt10cm ,ordchon10cm ,carchon10cm ,granite10cm ,sandstone10cm)
sprintf('Iron & %s \\\\ \n Acidic Basalt & %s \\\\ \n Basic Basalt & %s \\\\ \n Carbonaceous Chondrite & %s \\\\ \n Ordinary Chondrite & %s \\\\ \n Granite & %s \\\\ \n Sandstone & %s \\\\ \n',iron20cm ,basaltmax20cm ,minbasalt20cm ,ordchon20cm ,carchon20cm ,granite20cm ,sandstone20cm)
% T5cm = disp([iron5cm ,basaltmax5cm ,minbasalt5cm ,ordchon5cm ,carchon5cm ,granite5cm ,sandstone5cm]
% T10cm = disp([iron10cm ,basaltmax10cm ,minbasalt10cm ,ordchon10cm ,carchon10cm ,granite10cm ,sandstone10cm]
% T20cm = disp([iron20cm ,basaltmax20cm ,minbasalt20cm ,ordchon20cm ,carchon20cm ,granite20cm ,sandstone20cm]
