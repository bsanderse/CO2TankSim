function nphase = getnphase(U)
    % input U: rho, rho*e, T, possibly for different time instances:
    % U has size 3 x Nt
    % nphase = 1: single phase
    % nphase = 2: multi phase

    Nt = size(U,2);
    nphase = zeros(Nt,1);
    
    method = 2;
    
    switch method
        
        case 1
            ecurr = U(2,:)./U(1,:);
            
            % rho-e curve
            T_sample = linspace(CO2.Tt, CO2.Tc,100)';
            rhog_vap_sample = CO2.rhoVapSat(T_sample);
            eg_vap_sample   = CO2.u_rhoT(rhog_vap_sample,T_sample);
            rhol_vap_sample = CO2.rhoLiqSat(T_sample);
            el_vap_sample   = CO2.u_rhoT(rhol_vap_sample,T_sample);

            % interpolate e at current rho to find what corresponding e is
            for k=1:Nt
                if (U(1,k)>CO2.rhoc) % larger than critical value => use liquid curve
                    ecurr_vap = interp1(rhol_vap_sample,el_vap_sample,U(1,k));
                else
                    ecurr_vap = interp1(rhog_vap_sample,eg_vap_sample,U(1,k));
                end

                if (ecurr(k) < ecurr_vap)
                    nphase(k) = 2;
                else
                    nphase(k) = 1;
                end
            end
    
        case 2
            
            for k=1:Nt
                % single phase: 1 => 1
                % two phase: 0 => 2
                nphase(k) = 2 - CO2.isSinglePhase_rhoT(U(1,k),U(3,k));
            end

    end

end