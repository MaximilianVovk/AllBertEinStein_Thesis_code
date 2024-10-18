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

newcolors = [0 0.4470 0.7410
             0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330];
         
% for ii=1:3
%     figure(ii)
%     subplot(2,4,8)
%     plot(1,1,'b.')
%     hold on
%     plot(1,1,'r-')
%     legend('Mean Data','Gaussian Fit','location','west');
%     axis off  
% end

matrix_holder=[];
load('meteorMass.mat')
for ii=1:length(files)
    if files(ii).isdir==1 && ~isempty(regexp(files(ii).name,'m-','match'))               
        %% subplot
        
        if files(ii).init_size==5
              figure(1)
              colororder(newcolors)
                switch files(ii).material
                    case "iron"
                        subplot(2,4,1)
                        title(['Iron ' num2str(files(ii).mass(1)) ' kg'])
                    case "basaltmax"
                        subplot(2,4,2)
                        title(['Basalt Acidic ' num2str(files(ii).mass(1)) ' kg'])
                        holder_split=regexp(files(ii).name,'m-basaltmax5cm-','split');
                        bmax5cm=[bmax5cm convertCharsToStrings(holder_split{2})];
                    case "minbasalt"
                        subplot(2,4,3)
                        title(['Basalt Basic ' num2str(files(ii).mass(1)) ' kg'])
                    case "ordchon"
                        subplot(2,4,4)
                        title(['Ordinary Chondrite ' num2str(files(ii).mass(1)) ' kg'])
                    case "carchon"
                        subplot(2,4,5)
                        title(['Carbonaceous Chondrite ' num2str(files(ii).mass(1)) ' kg'])
                    case "granite"
                        subplot(2,4,6)
                        title(['Granite ' num2str(files(ii).mass(1)) ' kg'])
                    case "sandstone"
                        subplot(2,4,7)
                        title(['Sandstone ' num2str(files(ii).mass(1)) ' kg'])
                end
                plot([files(ii).mag(1:end-1,2); flip(files(ii).mag(1:end-1,1)); files(ii).mag(1)],[files(ii).alt(1:end-1); flip(files(ii).alt(1:end-1)); files(ii).alt(1)])
                ylim([0 100])
                xlim([-10 10])
                xlabel('Magnitude') 
                ylabel('Altitude [km]')
                set(gca, 'XDir','reverse')
                grid on
                hold on
                if ~isempty(files(ii).coef_avg)
                    crosSecAvg=ones(100,1)*min(files(ii).crosSec(files(ii).crosSec>=0));
                    rho_airAvg=ones(100,1)*min(files(ii).rho_air);
                    velAvg=ones(100,1)*min(files(ii).vel);
                    for jj=1:100
                        crosSecAvg(jj)=mean(files(ii).crosSec(find(round(files(ii).alt)==jj)));
                        if isnan(crosSecAvg(jj))
                            crosSecAvg(jj)=min(files(ii).crosSec(files(ii).crosSec>=0));
                        end
                        rho_airAvg(jj)=mean(files(ii).rho_air(find(round(files(ii).alt)==jj)));
                        if isnan(rho_airAvg(jj))
                            rho_airAvg(jj)=max(files(ii).rho_air);
                        end
                        velAvg(jj)=mean(files(ii).vel(find(round(files(ii).alt)==jj)));
                        if isnan(velAvg(jj))
                            velAvg(jj)=min(files(ii).vel);
                        end
                    end
                    mass_rateAvg = files(ii).coef_gauss_avg.*(crosSecAvg .* rho_airAvg .* (velAvg).^3) ./ (2*files(ii).Hmelt);
                    lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff ;
                    lum_eneAvg(lum_eneAvg==0)=nan;
                    magAvg = 6.8 - 1.086 * log(lum_eneAvg);
                    plot([magAvg(:,2); flip(magAvg(:,1)); magAvg(1)], [1:100 flip(1:100) 1], material_color(material_name==files(ii).material),'Linewidth',1.2)
                end
                
                sgtitle('Meteor of 5 cm')
                subplot(2,4,8)
                plot(linspace(0,1,100),nan)
                legend(bmax5cm,'location','west');
                axis off
                
        elseif files(ii).init_size==10
              figure(2)
              colororder(newcolors)
                switch files(ii).material
                    case "iron"
                        subplot(2,4,1)
                        title(['Iron ' num2str(files(ii).mass(1)) ' kg'])
                    case "basaltmax"
                        subplot(2,4,2)
                        title(['Basalt Acidic ' num2str(files(ii).mass(1)) ' kg'])
                        holder_split=regexp(files(ii).name,'m-basaltmax10cm-','split');
                        bmax10cm=[bmax10cm convertCharsToStrings(holder_split{2})];
                    case "minbasalt"
                        subplot(2,4,3)
                        title(['Basalt Basic ' num2str(files(ii).mass(1)) ' kg'])
                    case "ordchon"
                        subplot(2,4,4)
                        title(['Ordinary Chondrite ' num2str(files(ii).mass(1)) ' kg'])
                    case "carchon"
                        subplot(2,4,5)
                        title(['Carbonaceous Chondrite ' num2str(files(ii).mass(1)) ' kg'])
                    case "granite"
                        subplot(2,4,6)
                        title(['Granite ' num2str(files(ii).mass(1)) ' kg'])
                    case "sandstone"
                        subplot(2,4,7)
                        title(['Sandstone ' num2str(files(ii).mass(1)) ' kg'])
                end
                plot([files(ii).mag(1:end-1,2); flip(files(ii).mag(1:end-1,1)); files(ii).mag(1)],[files(ii).alt(1:end-1); flip(files(ii).alt(1:end-1)); files(ii).alt(1)])
                ylim([0 100])
                xlim([-10 10])
                xlabel('Magnitude') 
                ylabel('Altitude [km]')
                set(gca, 'XDir','reverse')
                grid on
                hold on
                if ~isempty(files(ii).coef_avg)
                    crosSecAvg=ones(100,1)*min(files(ii).crosSec(files(ii).crosSec>=0));
                    rho_airAvg=ones(100,1)*min(files(ii).rho_air);
                    velAvg=ones(100,1)*min(files(ii).vel);
                    for jj=1:100
                        crosSecAvg(jj)=mean(files(ii).crosSec(find(round(files(ii).alt)==jj)));
                        if isnan(crosSecAvg(jj))
                            crosSecAvg(jj)=min(files(ii).crosSec(files(ii).crosSec>=0));
                        end
                        rho_airAvg(jj)=mean(files(ii).rho_air(find(round(files(ii).alt)==jj)));
                        if isnan(rho_airAvg(jj))
                            rho_airAvg(jj)=max(files(ii).rho_air);
                        end
                        velAvg(jj)=mean(files(ii).vel(find(round(files(ii).alt)==jj)));
                        if isnan(velAvg(jj))
                            velAvg(jj)=min(files(ii).vel);
                        end
                    end
                    mass_rateAvg = files(ii).coef_gauss_avg.*(crosSecAvg .* rho_airAvg .* (velAvg).^3) ./ (2*files(ii).Hmelt);
                    lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff ;
                    lum_eneAvg(lum_eneAvg==0)=nan;
                    magAvg = 6.8 - 1.086 * log(lum_eneAvg);
                    plot([magAvg(:,2); flip(magAvg(:,1)); magAvg(1)], [1:100 flip(1:100) 1], material_color(material_name==files(ii).material),'Linewidth',1.2)
                end
                
                sgtitle('Meteor of 10 cm')
                subplot(2,4,8)
                plot(linspace(0,1,100),nan)
                legend(bmax10cm,'location','west');
                axis off
                
        elseif files(ii).init_size==20
              figure(3)
              colororder(newcolors)
                switch files(ii).material
                    case "iron"
                        subplot(2,4,1)
                        title(['Iron ' num2str(files(ii).mass(1)) ' kg'])
                    case "basaltmax"
                        subplot(2,4,2)
                        title(['Basalt Acidic ' num2str(files(ii).mass(1)) ' kg'])
                        holder_split=regexp(files(ii).name,'m-basaltmax20cm-','split');
                        bmax20cm=[bmax20cm convertCharsToStrings(holder_split{2})];
                    case "minbasalt"
                        subplot(2,4,3)
                        title(['Basalt Basic ' num2str(files(ii).mass(1)) ' kg'])
                    case "ordchon"
                        subplot(2,4,4)
                        title(['Ordinary Chondrite ' num2str(files(ii).mass(1)) ' kg'])
                    case "carchon"
                        subplot(2,4,5)
                        title(['Carbonaceous Chondrite ' num2str(files(ii).mass(1)) ' kg'])
                    case "granite"
                        subplot(2,4,6)
                        title(['Granite ' num2str(files(ii).mass(1)) ' kg'])
                    case "sandstone"
                        subplot(2,4,7)
                        title(['Sandstone ' num2str(files(ii).mass(1)) ' kg'])
                end
                plot([files(ii).mag(1:end-1,2); flip(files(ii).mag(1:end-1,1)); files(ii).mag(1)],[files(ii).alt(1:end-1); flip(files(ii).alt(1:end-1)); files(ii).alt(1)])
                ylim([0 100])
                xlim([-10 10])
                xlabel('Magnitude') 
                ylabel('Altitude [km]')
                set(gca, 'XDir','reverse')
                grid on
                hold on
                if ~isempty(files(ii).coef_avg)
                    crosSecAvg=ones(100,1)*min(files(ii).crosSec(files(ii).crosSec>=0));
                    rho_airAvg=ones(100,1)*min(files(ii).rho_air);
                    velAvg=ones(100,1)*min(files(ii).vel);
                    for jj=1:100
                        crosSecAvg(jj)=mean(files(ii).crosSec(find(round(files(ii).alt)==jj)));
                        if isnan(crosSecAvg(jj))
                            crosSecAvg(jj)=min(files(ii).crosSec(files(ii).crosSec>=0));
                        end
                        rho_airAvg(jj)=mean(files(ii).rho_air(find(round(files(ii).alt)==jj)));
                        if isnan(rho_airAvg(jj))
                            rho_airAvg(jj)=max(files(ii).rho_air);
                        end
                        velAvg(jj)=mean(files(ii).vel(find(round(files(ii).alt)==jj)));
                        if isnan(velAvg(jj))
                            velAvg(jj)=min(files(ii).vel);
                        end
                    end
                    mass_rateAvg = files(ii).coef_gauss_avg.*(crosSecAvg .* rho_airAvg .* (velAvg).^3) ./ (2*files(ii).Hmelt);
                    lum_eneAvg = 0.5 .* velAvg.^2 .* mass_rateAvg .* lum_eff ;
                    lum_eneAvg(lum_eneAvg==0)=nan;
                    magAvg = 6.8 - 1.086 * log(lum_eneAvg);
                    plot([magAvg(:,2); flip(magAvg(:,1)); magAvg(1)], [1:100 flip(1:100) 1], material_color(material_name==files(ii).material),'Linewidth',1.2)
                end
                
                sgtitle('Meteor of 20 cm')
                subplot(2,4,8)
                plot(linspace(0,1,100),nan)
                legend(bmax20cm,'location','west');
                axis off
        end
        
    else

    end
    
end

