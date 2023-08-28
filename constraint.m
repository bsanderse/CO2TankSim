function g = constraint(U)
    % U is size 3 x Nt
    nphase = getnphase(U);
    N = size(U,2);
    g = zeros(N,1);
    
    for k = 1:N
        if (nphase(k) == 1)
            g(k) = constraint_singlephase(U(:,k));
        elseif (nphase(k) == 2)
            V    = getVfromU(U(:,k));
            g(k) = constraint_twophase(V,U(:,k));            
        end
    end
        
    
end


function g_max = constraint_twophase(V,U)
    % U1 = rho, U2 = rho*e, V1=U3=T, V2=U4=p, V3=U5=alpha_g, V4=U6=rhog, V5=U7=rho_l
    
    g(1,1) = V(3)*V(4) + (1-V(3))*V(5) - U(1);
    g(2,1) = V(3)*V(4)*CO2.u_rhoT(V(4),V(1)) + (1-V(3))*V(5)*CO2.u_rhoT(V(5),V(1)) - U(2);
    g(3,1) = CO2.p_rhoT(V(4),V(1)) - CO2.p_rhoT(V(5),V(1));
    g(4,1) = CO2.G_rhoT(V(4),V(1)) - CO2.G_rhoT(V(5),V(1));
    g(5,1) = V(2) - CO2.p_rhoT(V(5),V(1));
    
    g_max = max(abs(g));
end


function g = constraint_singlephase(U)
    % U is size 3 x Nt
    g = abs(CO2.u_rhoT(U(1,:),U(3,:)) - U(2,:)./U(1,:)); 
    
end
