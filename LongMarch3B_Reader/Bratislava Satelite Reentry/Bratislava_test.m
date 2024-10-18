clc
clear all
close all

%%

load('sat1_6.3kms.mat')
load('sat2_6.7kms.mat')
load('sat3_7.1kms.mat')
load('sat4_7.5kms.mat')

allSat_fragm{1}=sat1;
allSat_fragm{2}=sat2;
allSat_fragm{3}=sat3;
allSat_fragm{4}=sat4;

for ii=1:4
    maxTime=max(allSat_fragm{ii}(:,1));
    for jj=0:maxTime
        if ~isempty(max(allSat_fragm{ii}(find(round(allSat_fragm{ii}(:,1))==(jj)),3)))
%             valuesmag=allSat_fragm{ii}(find(round(allSat_fragm{ii}(:,2))==(jj)),3);
%             [maxiMag,pos]=max(valuesmag));
%             valuesmag(pos)=min(valuesmag);
%             [maxiMag,pos]=max(valuesmag);
%             maxMag(jj+1)=max(valuesmag(find(valuesmag<maxiMag-0.2)));

            posmag=find(round(allSat_fragm{ii}(:,1))==(jj));
            valuesmag=allSat_fragm{ii}(posmag,3);
            valuesalt=allSat_fragm{ii}(posmag,2);
            [Mag,pos]=max(valuesmag);
            [Mag,pos]=max(valuesmag(find(valuesmag<Mag)));
            maxMag(jj+1)=valuesmag(pos);
            maxAlt(jj+1)=valuesalt(pos);
            
%             allSat_fragm{ii}(jj,4)=max(valuesmag);

%         elseif max(allSat_fragm{ii}(find(round(allSat_fragm{ii}(:,1))==(jj)),3))==0
%             maxMag(jj+1)=nan;
%         else
%             maxMag(jj+1)=nan;
        end
    end
    subplot(2,2,ii)
    
satellite1 = fit(maxAlt',maxMag','gauss1');
plot(satellite1,maxAlt',maxMag')
% legend('hide')
view(90,-90)
    
%     plot(allSat_fragm{ii}(:,3),allSat_fragm{ii}(:,2),'.b')
%     hold on
%     plot(maxMag,maxAlt,'r')%,'Linewidth',1.2)
        clear maxMag maxAlt

% satellite1 = fit(allSat_fragm{ii}(:,2),allSat_fragm{ii}(:,3),'gauss1');
% plot(satellite1,allSat_fragm{ii}(:,2),allSat_fragm{ii}(:,3))
% % legend('hide')
% view(90,-90)
    
%     ylim([0 100])
%     xlim([-10 10])
    xlabel('Magnitude') 
    ylabel('Altitude [km]')
    set(gca, 'YDir','reverse')
    grid on
    title(['Frafgment ' num2str(ii)])
%     ylim([0 100])
end

% smoothing=1;
% 
% subplot(2,2,1)
% 
% % plot(smooth(sat1(:,3),smoothing),sat1(:,2))
% % hold on
% figure
% alt=1:100;
% satellite1 = fit(alt',allSat_fragm{ii}(:,4),'gauss1');
% plot(satellite1,alt',allSat_fragm{ii}(:,4))
% % legend('hide')
% view(90,-90)
% title('Frafgment 1')
%     
% subplot(2,2,2)
% % plot(smooth(sat2(:,3),smoothing),sat2(:,2))
% % hold on
% % plot(smooth(sat2(:,3)),sat2(:,2),'r','Linewidth',1.2);%smooth
% % satellite2 = fit(smooth(sat2(:,2),smoothing),sat2(:,3),'gauss1');
% % plot(satellite2,smooth(sat2(:,2),smoothing),sat2(:,3))
% % % legend('hide')
% % view(90,-90)
% title('Frafgment 2')
% 
% subplot(2,2,3)
% % plot(smooth(sat3(:,3),smoothing),sat3(:,2))
% % hold on
% % plot(smooth(sat3(:,3)),sat3(:,2),'r','Linewidth',1.2);%smooth
% % satellite3 = fit(smooth(sat3(:,2),smoothing),sat3(:,3),'gauss1');
% % plot(satellite3,smooth(sat3(:,2),smoothing),sat3(:,3))
% % % legend('hide')
% % view(90,-90)
% title('Frafgment 3')
% 
% subplot(2,2,4)
% % plot(smooth(sat4(:,3),smoothing),sat4(:,2))
% % hold on
% % plot(smooth(sat4(:,3)),sat4(:,2),'r','Linewidth',1.2);%smooth
% % satellite4 = fit(smooth(sat4(:,2),smoothing),sat4(:,3),'gauss1');
% % plot(satellite4,smooth(sat4(:,2),smoothing),sat4(:,3))
% % % legend('hide')
% % view(90,-90)
% title('Frafgment 4')
% 
% sgtitle('CZ-3B R/B Reentry')
% 
% for ii=1:4
%     subplot(2,2,ii)
%     xlabel('Abs.Mag [-]') 
%     ylabel('Altitude [km]')
%     grid on
% %     ylim([0 100])
% end
% 
% % 
% % subplot(2,4,8)
% % plot(linspace(0,1,100),nan)
% % legend(bmax10cm,'location','west');
% % axis off

