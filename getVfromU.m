function V = getVfromU(U)
    % solve for V, given U(1)=rho, U(2) = rho*e and U(3)=T
    % assume U has size 3 x Nt
    % V1=U3=T, V2=U4=p, V3=U5=alpha_g, V4=U6=rhog, V5=U7=rho_l
    
    V = zeros(5,size(U,2));
    
    % T
    V(1,:) = U(3,:);
    % p
    V(2,:) = CO2.pVap(U(3,:));
    % rhol
    rhol   = CO2.rhoLiqSat(U(3,:));
    V(5,:) = rhol;
    % rhog
    rhog   = CO2.rhoVapSat(U(3,:));
    V(4,:) = rhog;
    
    % alpha_g
    V(3,:) = (rhol - U(1,:))./(rhol - rhog);
    
end
