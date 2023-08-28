function mdot = massflow(U,p,param) 
    % mass flow through valve
    mdot = param(4)*sqrt(U(1)*(p-param(2))*1e6).*((p-param(2))>0); % pressure in MPa, but mass flow should be in kg/s, so using 1e6 scaling
    
end