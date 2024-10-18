%% COSTANTS
Stefan_Boltzmann_cost=5.67*10^(-8); % W⋅m−2⋅K−4
Boltzmann_cost= 1.380649*10^(-23); %m2 kg s-2 K-1
sea_lv_pres=101325; %101,325 Pa
R_gas=8.31446261815324;
R_specif=287;
H=(R_gas*250)/(0.029*9.81); %scale height
avogadro_N=6.0221408e+23;

%Earth
r_planet = 6371000;
G = 6.67e-11;
M = 5.975e24;
g0_planet = 9.81;

%Observer position
obs.x_init = 100000;
obs.y_init = 100000.0;
obs.z_init = 10000.0;

%integration variables
dt = 0.025;
t_max = 300;

t=1;

%% Asteroid data
ast.init_v = 8000;  %Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
ast.v(1) = ast.init_v;
ast.init_zenith = 45;    %zenith angle
ast.v_vertical(1) = -ast.init_v * cos(deg2rad(ast.init_zenith));
ast.v_horizontal(1) = ast.init_v * sin(deg2rad(ast.init_zenith));
ast.init_pos = r_planet + 120000;
ast.h(1) = ast.init_pos - r_planet;
ast.s(1) = 0;

[rho_air,a_air,T_air,P_air,nu_air,h,sigma] = atmos(ast.h(t));
atmo_density = [0.09, -0.38, -1.05, -1.74, -2.39, -2.98, -3.50, -4.07, -4.72, -5.45,-6.30, -7.00, -7.62, -7.97, -8.32, -8.67, -8.81];
ast.init_T = 300;
ast.T = ast.init_T;
ast.melting_T=1800; %avg stone 1800K avg iron 1800K BOOK
ast.boiling_T=2960; %avg stone 2960K avg iron 3508K BOOK
ast.T_sec=0;

ast.r(1) = 0.1; %in meters
ast.rho = 7800.0;   % density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
ast.A = (pi * ast.r(1)^2)/(4/3 * pi * ast.r(1)^3)^(2/3);% 1.2;    % shape factor (1.2 for sphere). %A=S/V^2/3
ast.drag_coeff = 0.5; % due to aircap actual Cd=2drag_coef
ast.heat_transfer = 0.02; %[-] kg/m^2/s vary with vel 2^-4 to 7^-3 and altitude heat transfer coefficient, between 0.05 and 0.6. and 0.5 Brown_Detlef_2004
ast.heat_of_ablation = 8000000;   %J/kg latentheat of melting q heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg. BOOK NEW evaporat 8e6 liquid 2e6
ast.heat_of_ablation_fus = 2000000; %289000 iron
ast.specif_heat = 895; %Stone 895 J/kg/K liquid 1100 J/kg.K iron 691 liquid 666 HEAT COND 689 J/kg·K and 1314 J/kg·K
ast.thermal_conductivity= 0.4;% of rocks falls usually in the range of 0.40–7.00
ast.emissivity=0.9; %Brown_Detlef_2004
ast.condensation_coef=0.5; %0.5 for stone and 1 for metals Brown_Detlef_2004
ast.vap_pres=0;
ast.surface_tension=1200e-3; %N/m
ast.radiation_abs_rate=3; %N/m/s
ast.sup_tension=1.23;

ast.ultim_tensile=2.25e7; %3.5*10^7;%1.7*10^7; avg for meteoroids
ast.ultim_compres=1.7*10^7;%1.5e8;
ast.dynamic_pres(1)=0.365*rho_air*ast.v(1)^2;
ast.stagn_pres(1)=rho_air*ast.v(1)^2;

ast.init_mass = 4/3 * pi * ast.r^3 * ast.rho;  % mass of projectile sphere in kg
ast.m(1)=ast.init_mass;
ast.lum_eff = 0.213 * ast.init_mass^(-0.17);      % 0.1 - 0.006
ast.ballistic_coeff = ast.init_mass / (ast.drag_coeff * pi * ast.r^2); % normal ballistic coeff wiki
ast.m_atom=23*1.6605402E-27; %mean atomic weight 23 mean molec weight of vapours 5 BOOK Pysics of meteor
ast.m_sec=0;