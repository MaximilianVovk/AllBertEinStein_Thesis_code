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

for ii=1:3
    figure(ii)
    subplot(2,4,8)
    plot(1,1,'b.')
    hold on
    plot(1,1,'r-')
    legend('Mean Data','Gaussian Fit','location','west');
    axis off  
end

matrix_holder=[];
load('meteorMass.mat')
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))
        for jj=1:100
            heatcoef(jj)=mean(files(ii).heat_transf_coef(find(round(files(ii).alt)==jj)));
            if isnan(heatcoef(jj))
                heatcoef(jj)=0;
            end
        end
        files(ii).heat_transf_coef100km=heatcoef';
        alti=1:100;
        coefGauss = fit(alti',files(ii).heat_transf_coef100km,'gauss1'); 
        files(ii).gauss_StD_Mean=coefGauss;
        files(ii).coef_gauss=coefGauss(alti);
        
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
                subplot(2,4,1)
                title('Iron Gauss')
                if ~isempty(files(ii).coef_avg)
                    alti=1:100;
                    ironGauss5 = fit(alti',files(ii).coef_avg,'gauss1');
                    plot(ironGauss5,alti,files(ii).coef_avg)
                    legend('hide')
                    view(90,-90)
                    
                    files(ii).coef_gauss_avg=ironGauss5(alti);
                    if files(ii).init_size==10
                        ironGauss10=ironGauss5;
                    elseif files(ii).init_size==20
                        ironGauss20=ironGauss5;
                    end
                end
            case "basaltmax"
                subplot(2,4,2)
                title('Basalt Acidic Gauss')
                if ~isempty(files(ii).coef_avg)
                    alti=1:100;
                    basaltmaxGauss5 = fit(alti',files(ii).coef_avg,'gauss1');
                    plot(basaltmaxGauss5,alti,files(ii).coef_avg)
                    legend('hide')
                    view(90,-90)
                    
                    files(ii).coef_gauss_avg=basaltmaxGauss5(alti);
                    if files(ii).init_size==10
                        basaltmaxGauss10=basaltmaxGauss5;
                    elseif files(ii).init_size==20
                        basaltmaxGauss20=basaltmaxGauss5;
                    end
                end
            case "minbasalt"
                subplot(2,4,3)
                title('Basalt Basic Gauss')
                if ~isempty(files(ii).coef_avg)
                    alti=1:100;
                    basaltminGauss5 = fit(alti',files(ii).coef_avg,'gauss1');
                    plot(basaltminGauss5,alti,files(ii).coef_avg)
                    legend('hide')
                    view(90,-90)
                    
                    files(ii).coef_gauss_avg=basaltminGauss5(alti);
                    if files(ii).init_size==10
                        basaltminGauss10=basaltminGauss5;
                    elseif files(ii).init_size==20
                        basaltminGauss20=basaltminGauss5;
                    end
                end
            case "ordchon"
                subplot(2,4,4)
                title('Ordinary Chondrite Gauss')
                if ~isempty(files(ii).coef_avg)
                    alti=1:100;
                    ordchonGauss5 = fit(alti',files(ii).coef_avg,'gauss1');
                    plot(ordchonGauss5,alti,files(ii).coef_avg)
                    legend('hide')
                    view(90,-90)
                    
                    files(ii).coef_gauss_avg=ordchonGauss5(alti);
                    if files(ii).init_size==10
                        ordchonGauss10=ordchonGauss5;
                    elseif files(ii).init_size==20
                        ordchonGauss20=ordchonGauss5;
                    end
                end
            case "carchon"
                subplot(2,4,5)
                title('Carbonaceous Chondrite Gauss')
                if ~isempty(files(ii).coef_avg)
                    alti=1:100;
                    carchonGauss5 = fit(alti',files(ii).coef_avg,'gauss1');
                    plot(carchonGauss5,alti,files(ii).coef_avg)
                    legend('hide')
                    view(90,-90)
                    
                    files(ii).coef_gauss_avg=carchonGauss5(alti);
                    if files(ii).init_size==10
                        carchonGauss10=carchonGauss5;
                    elseif files(ii).init_size==20
                        carchonGauss20=carchonGauss5;
                    end
                end
            case "granite"
                subplot(2,4,6)
                title('Granite Gauss')
                if ~isempty(files(ii).coef_avg)
                    alti=1:100;
                    graniteGauss5 = fit(alti',files(ii).coef_avg,'gauss1');
                    plot(graniteGauss5,alti,files(ii).coef_avg)
                    legend('hide')
                    view(90,-90)
                    
                    files(ii).coef_gauss_avg=graniteGauss5(alti);
                    if files(ii).init_size==10
                        graniteGauss10=graniteGauss5;
                    elseif files(ii).init_size==20
                        graniteGauss20=graniteGauss5;
                    end
                end
            case "sandstone"
                subplot(2,4,7)
                title('Sandstone Gauss')
                if ~isempty(files(ii).coef_avg)
                    alti=1:100;
                    sandstoneGauss5 = fit(alti',files(ii).coef_avg,'gauss1');
                    plot(sandstoneGauss5,alti,files(ii).coef_avg)
                    legend('hide')
                    view(90,-90)
                    
                    files(ii).coef_gauss_avg=sandstoneGauss5(alti);
                    if files(ii).init_size==10
                        sandstoneGauss10=sandstoneGauss5;
                    elseif files(ii).init_size==20
                        sandstoneGauss20=sandstoneGauss5;
                    end
                end
        end
        xlim([0 100])
        ylabel('Heat transfer coef.') 
        xlabel('Altitude [km]')
        legend('hide')
        grid on
        hold on
        if files(ii).init_size==5
            sgtitle('Meteor of 5 cm')
        elseif files(ii).init_size==10
            sgtitle('Meteor of 10 cm')
        elseif files(ii).init_size==20
            sgtitle('Meteor of 20 cm')
        end
                
        
    else

    end
    
end

% save('meteorMass.mat','files')