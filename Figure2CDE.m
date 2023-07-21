%% Code generating Figures 2C, 2D, 2E

% Cleaning 
clear all, close all, clc;

% Maximum input-current values (C,D,E respectively) 
IMaxVals = [1.1 1.8 1.8];

% Initial conditions (C,D,E respectively) 
IC_C = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
IC_D = [0.050, 0.00, 0.062, 0.046, 0.035, 0.35];
IC_E = [0.005, 0.00, 0.052, 0.03, 0.08, 0.35];
ICVals = [IC_C; IC_D; IC_E];



for i = 1:length(IMaxVals)      % For each input current values...

    
    %% External inputs (as a vector)  
    I1fun = @(t) (t<=500).*IMaxVals(i);
    I = @(t) [I1fun(t); 0*t; 0*t; 0*t; 0*t; 0*t];       %Applied current (only) in V1
    
    %% Connectivity matrix
    W = [... %V1E  PPCE     PFCE         V1I     PPCI     PFCI
             1.0,   11.22,  1.29,       -2.3,    0,     0;...    % V1E
             4.57,  1.0,    10.57,       0,     -1.8,   0;...    % PPCE
             0.72,  9.87,   1.0,         0,      0,    -1.9;...  % PFCE

             2.0,    0,      0,          0.5,    0,     0;...    % V1I
             0,      2.0,    0,          0,      0.5,   0;...    % PPCI
             0,      0,      2.0,        0,      0,     0.5;     % PFCI
    ];
    
    %% Gain functions excitatory & inhibitory parts
    mu1E= 3;  nu1E =2.0;
    mu2E= 2;  nu2E =4.0;
    mu3E= 2;  nu3E =2.0;
    mu1I= 2;  nu1I =0.3;
    mu2I= 2;  nu2I =0.3;
    mu3I= 2;  nu3I =0.3;
    
    % FIRING RATE FUNCTION
    fr = @(v,mu,nu) 1/(1+exp(-mu*(v-nu)));       
    S = @(v) [ fr(v(1),mu1E,nu1E);  fr(v(2),mu2E,nu2E); fr(v(3),mu3E,nu3E); fr(v(4),mu1E,nu1E); fr(v(5),mu2E,nu2E); fr(v(6),mu3E,nu3E)]; 

    %% Time constants
    tauE = [30; 66.6;  38]; 
    tauI = [10; 10;  10];
    tau = [tauE; tauI];

    %% Decay constants
    betaE = [0.8; 0.3; 0.8]; betaI = [0.07; 0.1; 0.07];
    beta = [betaE; betaI];

    %% Right-hand side
    F = @(t,u)(-beta.*u + S(W*u + I(t)))./tau;

    %% Time step
    tSpan = [0:1:1500];
    [t,U] = ode23s(F,tSpan, ICVals(i,:));

    %% Plot trajectory
    figure, hold on;
    for subl = 1:6
       subplot(6,1,subl)
       plot(t,U(:,subl)); drawnow;
    end
    subplot(6,1,1)
        
    titles_list = ["Figure 2C"; "Figure 2D"; "Figure 2E"];
    title([titles_list(i)])
    subtitle(['IMax = ' num2str(IMaxVals(i))]);
    hold off;

end







