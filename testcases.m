% function params = testcase(test_ID)

%% test cases
switch test_case_ID
    case 1
        % params Hammer
        disp('test case Hammer');
        T_0 = 300;
        p_0 = 100*bar_to_MPA;

        t_0 = 0;
        t_end = 200; % in seconds

        T_amb = 293.15;
        p_amb = 1*bar_to_MPA;
        etaA  = 1*J_to_kJ; % from W/K to kW/K

        Kv = 5e-7; % in m2
        V = pi/100;

        dt = 1;
    case 2
        % params Giljarhus
        disp('test case Giljarhus');
        T_0 = 300; %25 + C_to_K;
        p_0 = 80*bar_to_MPA;

        t_0 = 0;
        t_end = 200; %0.5*3600; % in seconds

        T_amb = 5 + C_to_K;
        p_amb = 6*bar_to_MPA;
        etaA  = 10*J_to_kJ; % from W/K to kW/K

        Kv = 8e-7;
        V = pi/100; % = (pi/4)*(0.2^2)
        
        dt = 1;
    otherwise
        
        error('wrong test case ID');
end
