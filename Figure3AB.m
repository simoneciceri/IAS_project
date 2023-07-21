%% Code generating Figures 3A, 3B

% Cleaning 
clear all, close all, clc;

% Maximum input-current values (C,D,E respectively) 
IMaxVals = [1.1 1.8 2 3];
nI = 4; 

% Number of random initial conditions
nIC =   20;    

% Container storing number of spikes of each trajectory
v1Int = zeros(nI,nIC);           

% Plotting options
plotTrajectories = true;
plotHistograms = true;


for i = 1:length(IMaxVals)     % For each input current values...

    %% Random initial conditions
    vrand = 0.05;             % Uniform distribution width
    ICVals = vrand*[rand(1,nIC); rand(1,nIC); rand(1,nIC); rand(1,nIC); rand(1,nIC); rand(1,nIC)];
    
    for ik = 1:nIC    % For each (random) initial condition

        %% External inputs (as a vector)  
        I1fun = @(t) (t>30 & t<=500).*IMaxVals(i);
        I = @(t) [I1fun(t); 0*t; 0*t; 0*t; 0*t; 0*t];       %Applied current (only) in V1

        %% Connectivity matrices
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

        % Firing rate function
        fr = @(v,mu,nu) 1/(1+exp(-mu*(v-nu)));       
        S = @(v) [ fr(v(1),mu1E,nu1E);  fr(v(2),mu2E,nu2E); fr(v(3),mu3E,nu3E); fr(v(4),mu1E,nu1E); fr(v(5),mu2E,nu2E); fr(v(6),mu3E,nu3E)]; 

        %% Time constants
        tauE = [30; 66.6;  38]; 
        tauI = [10; 10; 10];
        tau = [tauE; tauI];

        %% Decay constants
        betaE = [0.8; 0.3; 0.8]; 
        betaI = [0.07; 0.1; 0.07];
        beta = [betaE; betaI];

        %% Right-hand side
        F = @(t,u)(-beta.*u + S(W*u + I(t)))./tau;

            
        %% Time step
        u0 = [0. 0. 0. 0. 0. 0.];
        tSpan = [-50 0];
        [t,U] = ode23s(F,tSpan,u0);
    
        tSpan = [0:1:1500];
        [t,U] = ode23s(F,tSpan,U(end,:)' + ICVals(:,ik));   % add noise

        %% Integral of V1E between t0 and 150
        t0 = 250;                    
        dt = tSpan(2)-tSpan(1); 
        idt = find(t>=t0);
        v1Int(i,ik) = sum(U(idt,1))*dt / 1000;

        %% Plot trajectory
        if plotTrajectories       
           figure(1);
           subplot(4,1,5-i); hold on;
           plot(t,U(:,1)); drawnow;
        end

    end
    %% Plot histograms
    if plotHistograms
        figure(2);
        subplot(4,1,5-i); hold on;
        h = histogram(v1Int(i,:),50, 'BinLimits',[0,0.4]); drawnow;
    end
     
end


if plotTrajectories 
    figure(2);
    linkaxes([subplot(4,1,1), subplot(4,1,2), subplot(4,1,3), subplot(4,1,4)],'xy')
    hold off;
end
if plotHistograms 
    figure(2);
    linkaxes([subplot(4,1,1), subplot(4,1,2), subplot(4,1,3), subplot(4,1,4)],'xy')
    hold off;
end




