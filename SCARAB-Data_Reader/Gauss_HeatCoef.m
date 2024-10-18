function [HeatCoef] = Gauss_HeatCoef(altitude,velocity,size,material)
functionval=fit([0;1;2],[1;2;3],'gauss1');
load('meteorGaussCoef.mat')

switch material
    case "iron"
        a1fun=coefGauss.iron_a1(velocity, size);
        b1fun=coefGauss.iron_b1(velocity, size);
        c1fun=coefGauss.iron_c1(velocity, size);
    case "basaltmax"
        a1fun=coefGauss.basaltmax_a1(velocity, size);
        b1fun=coefGauss.basaltmax_b1(velocity, size);
        c1fun=coefGauss.basaltmax_c1(velocity, size);
    case "minbasalt"
        a1fun=coefGauss.basaltmin_a1(velocity, size);
        b1fun=coefGauss.basaltmin_b1(velocity, size);
        c1fun=coefGauss.basaltmin_c1(velocity, size);
    case "ordchon"
        a1fun=coefGauss.ordchon_a1(velocity, size);
        b1fun=coefGauss.ordchon_b1(velocity, size);
        c1fun=coefGauss.ordchon_c1(velocity, size);
    case "carchon"
        a1fun=coefGauss.carchon_a1(velocity, size);
        b1fun=coefGauss.carchon_b1(velocity, size);
        c1fun=coefGauss.carchon_c1(velocity, size);
    case "granite"
        a1fun=coefGauss.granite_a1(velocity, size);
        b1fun=coefGauss.granite_b1(velocity, size);
        c1fun=coefGauss.granite_c1(velocity, size);
    case "sandstone"
        a1fun=coefGauss.sandstone_a1(velocity, size);
        b1fun=coefGauss.sandstone_b1(velocity, size);
        c1fun=coefGauss.sandstone_c1(velocity, size);
end

functionval.a1=a1fun;
functionval.b1=b1fun;
functionval.c1=c1fun;
HeatCoef=abs(functionval(altitude));
HeatCoef(HeatCoef<0.01)=0;
HeatCoef(HeatCoef>1)=1;
% HeatCoef(HeatCoef<0)=0;
end