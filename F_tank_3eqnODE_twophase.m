function F = F_tank_3eqnODE_twophase(~,U,param)
% input U: rho, rho*e, T

rho  = U(1);
T    = U(3);

F    = zeros(3,1); 

p    = CO2.pVap(T);    


% get mass and heat flow
mdot = massflow(U,p,param);
Qdot = heatflow(U,param); % in kJ/s

% enthalpy
h      = (U(2) + p*1e3)/rho; % note 1e3 to go from MPa(=10^6 kg/(ms^2) to kJ/m3=10^3 kg /(ms^2)

% tank equations
F(1) = -mdot/param(5);
F(2) = (Qdot - mdot*h)/param(5);

% define perturbations needed for numerical derivatives
% here, eps is machine precision as defined by Matlab
T_pert   = T*sqrt(eps);
rho_pert = rho*sqrt(eps);

[f_g_rhoT, f_l_rhoT] = f_gl(T,rho);
[f_g_Tpert, f_l_Tpert] = f_gl(T+T_pert,rho);
[f_g_rhopert, f_l_rhopert] = f_gl(T,rho+rho_pert);

dfgdT  = (f_g_Tpert - f_g_rhoT)/T_pert;
dfldT  = (f_l_Tpert - f_l_rhoT)/T_pert;
dfgdU  = (f_g_rhopert - f_g_rhoT)/rho_pert;
dfldU  = (f_l_rhopert - f_l_rhoT)/rho_pert;

F(3) = (F(2) - (dfgdU + dfldU)*F(1))/(dfgdT + dfldT);



end

function [f_g,f_l] = f_gl(T,U1)
    rhol = CO2.rhoLiqSat(T);
    rhog = CO2.rhoVapSat(T);
    eg   = CO2.u_rhoT(rhog,T);
    el   = CO2.u_rhoT(rhol,T);
    alphag = (rhol - U1)/(rhol - rhog);
    alphal = (1-alphag);
    f_g  = rhog*eg*alphag;
    f_l  = rhol*el*alphal;
end
