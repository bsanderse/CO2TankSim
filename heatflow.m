function Qdot = heatflow(U,param) 
    % heatflow through tank
    Qdot = param(3)*(param(1) - U(3)); % note: U(3) in K, param(3) in kW/K, result is in kW = kJ/s
    
end