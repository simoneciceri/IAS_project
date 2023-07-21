%% Code generating Figure 3C
% Note: it might take some time. Reducing nI and nIC makes it faster



clear all, close all, clc;

% Maximum input values
nI =    20; 
IMaxVals = linspace(1,3,nI);

% Number of random initial conditions
nIC =   50;      

% Container storing number of spikes of each trajectory
v1Int = zeros(nI,nIC);    

% Plotting options
plotTrajectories = false;
plotHistograms = false;


for i = 1:length(IMaxVals)      % For each input current values...

    %% Random initial conditions
    vrand = 0.05;             % Uniform distribution width
    ICVals = vrand*[rand(1,nIC); rand(1,nIC); rand(1,nIC); rand(1,nIC); rand(1,nIC); rand(1,nIC)];
    
    if plotTrajectories 
        figure(2), hold on; title(['IMax = ' num2str(IMaxVals(i))]);
    end

    for ik = 1:length(ICVals)   % For each (random) initial condition

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

        
        if plotTrajectories 
            figure(2);
            plot(t,U(:,1)); drawnow;
        end
    end
    
    if plotTrajectories
        hold off;
    end

  
    if plotHistograms
        figure(3); 
        title(['IMax = ' num2str(IMaxVals(i))]);
        h = histogram(v1Int(i,:),20); drawnow;
        pause
    end

end

%% Surface plot
n_bins = 50;
y = linspace( 0 , max(v1Int, [], 'all')+0.05, n_bins);

for i = 1:length(IMaxVals)
    figure(3);
    h = histogram(v1Int(i,:), y);
    z(:, i) = h.Values;
end
y(end) = [];
[X,Y] = meshgrid(IMaxVals, y);

% Probability of counting a certain amount of spikes
spike_prob = zeros(nI,1);
for i = 1:nI
    spike_prob(i) = sum( v1Int(i,:)>0.2 & v1Int(i,:)<0.35 ) / nIC;
end

% Plot the probability
figure(1), hold on;
surf(X,Y,z); 
colormap(flipud(gray));
colorbar;
drawnow;
hold off;

if not(plotHistograms)
    close(figure(3));
end
