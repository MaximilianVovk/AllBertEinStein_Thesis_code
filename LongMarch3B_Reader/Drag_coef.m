function Cd = Drag_coef(Kn,A)
    switch nargin
        case 2
            Cd=(0.925+atan(Kn*35).^5*0.113).*(A/1.2).^(1/5);
        case 1
            Cd=0.925+atan(Kn*35).^5*0.113;
%         otherwise
%             Cd=0.925+atan(Kn*35).^5*0.113;
    end
% %% DATA
% Cd_data=[0.925 0.925 0.925 1 1.3 1.5 1.75 2 2 2 2 2 2 2 2 2 2];
% Kn_data=[0.001 0.002 0.01 0.035 0.1 0.15 0.35 1 1.5 2 3 5 6 7 8 9 10];
% 
% Cd=0.925+atan(Kn*35).^5*0.113;
end

% clc
% clear all
% close all
% 
% %%
% Cd=[0.925 0.925 0.925 1 1.3 1.5 1.75 2 2 2 2 2 2 2 2 2 2];
% Kn=[0.001 0.002 0.01 0.035 0.1 0.15 0.35 1 1.5 2 3 5 6 7 8 9 10];
% 
% DragFunc = fit(Kn',Cd','smoothingspline');
% % plot(DragFunc,Kn',Cd')
% % hold on
% semilogx(Kn,Cd)
% hold on
% Kn_more=0.01:0.001:1;
% semilogx([0.001 0.01 Kn_more 1 10],[0.925; 0.925; DragFunc(Kn_more);2; 2])
% Kn_more=0:0.001:1000;
% semilogx(Kn_more,0.925+atan(Kn_more*35).^5*0.113)
% % semilogx(Kn_more,0.925+atan(Kn_more*7.5)*0.69)
% R=2;
% L=1;
% R/L
% A=(pi * R^2)/(pi * R^2*L)^(2/3)
% % A=(pi * R^2)/(4/3 * pi * R^2*L)^(2/3);
% 
% semilogx(Kn_more,(0.925+atan(Kn_more*35).^5*0.113)*(A/1.2)^(1/5))
% grid on
% ylim([0.5 2.5])
% xlim([0 100])
% % atan(0.001:10)