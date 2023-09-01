function F = F_tank_3eqnODE_singlephase(~,U,param)
% input U: rho, rho*e, T

rho  = U(1);
T    = U(3);

F    = zeros(3,1); 

p    = CO2.p_rhoT(rho,T); % in MPa    


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

      
% get approximate derivatives of: de/drho and de/dT
e_rhoT = CO2.u_rhoT(rho,T);
dedrho = (CO2.u_rhoT(rho+rho_pert,T) - e_rhoT)/rho_pert;
dedT   = (CO2.u_rhoT(rho,T+T_pert) - e_rhoT)/T_pert;

F(3) = (F(2)/rho - (U(2)/(rho^2) + dedrho)*F(1))/dedT;
    


end

