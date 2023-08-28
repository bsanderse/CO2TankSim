% error plots
% Giljarhus test case, integrated until t=0.5 hours

% RK44 with dt=0.1
% U_ref = 1.0e+02 * ...
%     [0.118779029018086; ...
%     -8.730451803622472;...
%      2.781020589159342];

% ode45 with 1e-8 tol
U_ref = 1.0e+02 * ...
   [0.118778635811714; ...
  -8.778150112520649; ...
   2.781028605759860];

 
% dt_list 
dt_list = [0.5; 1; 2; 4; 8]; %; 16; 32];

% RK4
U_list = 1.0e+02 * ...
   [0.118780703012627   0.118778789084792   0.118782836543112   0.118780301531532   0.119033147777435; ...
  -8.528621956395915  -8.757719486983769  -8.276005183463706  -8.212041861048222  -1.870874598780055; ...
   2.780986545850507   2.781028090452308   2.780938918962780   2.780903688592661   2.776283464349440];

% RK2
% rho_list(:,2) = 1e2*[0.118755611935986; 0.118755611935986; 0.11875480645087;  0.118749211432196 ];
% rhoe_list(:,2) = 1e2*[-8.450452201674622; -8.450452201674622; -11.60003720634215;  -8.934812182265986];
% T_list(:,2) = 1e2*[2.781499260514579; 2.781499260514579;  2.78149983407106;  2.781499396862658];


%% errors

error = abs(U_list - U_ref);
% error_rho = abs(rho_list - rho_ref);
% error_T = abs(T_list - T_ref);


%% plots
figure(10)
loglog(dt_list,error(1,:),'s-','LineWidth',2);
grid on
xlabel('\Delta t');
ylabel('error in temperature at t=1800 s');
legend('RK4','RK2')
set(gca,'FontSize',16);

figure(11)
loglog(dt_list,error(3,:),'s-','LineWidth',2);
grid on
xlabel('\Delta t');
ylabel('error in density at t=1800 s');
legend('RK4','RK2')
set(gca,'FontSize',16);