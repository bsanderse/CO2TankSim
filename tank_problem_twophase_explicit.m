% tank problem

%%
addpath('/Users/sanderse/Dropbox/work/Programming/libs/SpanWagner/');
addpath('/Users/sanderse/Dropbox/work/Programming/RungeKutta/main/RKproperties');

close all
clearvars


%% select test case and time integration method
ODE_type = 1; % 1: 3 eqn ODE, own implementation; 2: 3 eqn ODE, Matlab solver
RK_method = 'RK44'; % in case ODE_type = 2: e.g. 'FE11','RK44',

test_case_ID = 2; % 1: Hammer; 2: Giljarhus, see testcases.m

%% conversions
unit_conversions;

%% load test case parameters
testcases;

% store parameters in vector
param(1) = T_amb;
param(2) = p_amb;
param(3) = etaA;
param(4) = Kv;
param(5) = V;

%% get consistent initial conditions for rho, rho*e and T that satisfy the constraint
rho_0 = CO2.rho_pT(p_0,T_0); % in kg/m3
e_0   = CO2.u_rhoT(rho_0,T_0); % in kJ/kg
U_0   = [rho_0;rho_0*e_0;T_0];


%% start time integration

dt_list = [1 2 4 8];
Nsim = length(dt_list);
U_list = zeros(3,Nsim);



switch ODE_type
 
    
    case 1 % own RK method
        
        for kk = 1:Nsim
            dt = dt_list(kk)
            
            disp('integrating ODE with own RK method');
            options.RK_method.name = RK_method;
            [A_RK,b_RK,c_RK] = getRKmethod(options.RK_method.name);
            RK_order         = check_orderconditions(A_RK,b_RK,c_RK,false);
            
            disp(['RK method: ' RK_method ', order ' num2str(RK_order)]);
            
            s = length(b_RK);% number of stages
            
            N = length(U_0);
            Nt = (t_end - t_0)/dt;
            
            U  = zeros(N,Nt+1);
            t  = zeros(Nt+1,1);
            
            nphase    = zeros(Nt+1,1);
            nphase(1) = getnphase(U_0);
            
            U(:,1)  = U_0;
            t(1)    = t_0;
            it_time = 1;
            
            while ( t(it_time) - t_end < -dt/100)
                
                
                Un   = U(:,it_time);
                tn   = t(it_time);
                F_RK = zeros(s,N);
                
                for j=1:s
                    
                    % intermediate stage value of stage j
                    Uj   = Un + dt*(A_RK(j,:)*F_RK).';
                    
                    % time level of this stage
                    tj   = tn + c_RK(j)*dt;
                    
                    % flux evaluation
                    F_RK(j,:) = F_tank_3eqnODE(tj,Uj,param);
                    
                end
                
                % update solution with the b-coefficients
                Unew = Un + dt*(b_RK.'*F_RK).';
                tnew = tn + dt;
                
                % update iteration counter
                it_time = it_time + 1;
                
                
                if (isnan(Unew))
                    error('NaN encountered');
                end
                
                
                % store solution in entire vector
                t(it_time)   = tnew;
                U(:,it_time) = Unew;
                nphase(it_time) = getnphase(Unew);
                
                if (it_time >1 && (nphase(it_time) ~= nphase(it_time-1)))
                    disp('phase switching');
                end
                
            end
            
            if (Nt+1>it_time) % truncate empty part of vector
                U = U(:,1:it_time);
                t = t(1:it_time);
                nphase = nphase(1:it_time);
            end
            
            disp(['simulation finished at t=' num2str(tnew)]);
            
            U_list(:,kk) = U(:,end);
        end
        
        
        
        
    case 2 % Matlab ODE
        
        disp('integrating ODE with matlab ode45 method');
        
        
        fcn = @(t,U) F_tank_3eqnODE(t,U,param);
        options = odeset('RelTol',1e-8,'AbsTol',1e-8);
        [t, U] = ode23s(fcn,[t_0 t_end], U_0, options);
        U = U.';
        
        
    otherwise
        error('wrong value for ODE_type');
        
        
        
end


%% postprocess
T_out = U(3,:);
% get pressure from rho and T
p_out = getpressure(U);

% evaluate accuracy with which constraint is satisfied
g_out = constraint(U);

V_out = getVfromU(U);
alphag_out = min(max(V_out(3,:),0),1); %.*(nphase == 2);




%% plot solution in p-T space
figure(1)

% p-T
plot(T_out,p_out/bar_to_MPA,'DisplayName','ODE simulation','Linewidth',2)
hold on
% plot vapour line
T_plot = [linspace(CO2.Tt-50, CO2.Tt, 50) linspace(CO2.Tt, 1100, 500)];
plot(T_plot, CO2.pVap(T_plot)/bar_to_MPA, '--','DisplayName', 'Saturation line','Linewidth',2);
grid on
xlabel('T [K]');
ylabel('p [bar]');
% legend('ODE simulation','saturation line');
legend('FontSize',14,'Location','NorthWest')
set(gca,'FontSize',16);



%% constraint accuracy
figure(2)
plot(t,g_out,'Linewidth',2);
xlabel('t [s]');
ylabel('constraint accuracy');
grid on
set(gca,'FontSize',16);

%% hold-up
figure(3)
plot(t,alphag_out,'Linewidth',2);
xlabel('t [s]');
ylabel('\alpha_{g} [-]');
grid on
set(gca,'FontSize',16);

%% rho vs. e
figure(4)
T_sample = linspace(CO2.Tt, CO2.Tc,100)';
rhog_vap_sample = CO2.rhoVapSat(T_sample);
eg_vap_sample = CO2.u_rhoT(rhog_vap_sample,T_sample);
rhol_vap_sample = CO2.rhoLiqSat(T_sample);
el_vap_sample = CO2.u_rhoT(rhol_vap_sample,T_sample);


semilogy(eg_vap_sample,rhog_vap_sample,'LineWidth',2,'DisplayName','gas');
hold on
semilogy(el_vap_sample,rhol_vap_sample,'LineWidth',2,'DisplayName','liquid');
semilogy(U(2,:)./U(1,:),U(1,:),'Linewidth',2,'DisplayName','simulation path');
grid on
set(gca,'FontSize',16);
xlabel('e [kJ/kg]');
ylabel('rho [kg/m3]');

%% make plot like in Giljarhus paper
figure(5)
yyaxis left
plot(t/3600,p_out/bar_to_MPA,'DisplayName','p [bar]','Linewidth',2)
hold on
plot(t/3600,T_out - C_to_K,'DisplayName','T [C]','Linewidth',2)

yyaxis right
plot(t/3600,alphag_out,'DisplayName','\alpha_g [-]','Linewidth',2);
legend
grid on
xlabel('t [hour]');
set(gca,'FontSize',16);

