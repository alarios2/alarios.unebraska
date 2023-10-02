% Zombie model with quarantines and latent infection.
% Model is from 2009 paper of Munz, Hudea, Imad, and Smith 
% published in Infectious Disease Modeling Research Papers.
close all; clear all

S_0 = 500; % Initial Suseptible (Human) Population
I_0 = 0;   % Initial Infected Population 
Z_0 = 0;   % Initial Zombie Population
R_0 = 0;   % Initial Removed Population (Dead)
Q_0 = 0;   % Initial Quantined Population

T = 4000; % Final time

% The paper does not give values for the following parameters, 
% so some of these are just guesses.
Pi    = 0.0001; % Human birth rate (maybe proportional to delta)
alpha = 0.005;  % Humans killing zombies coefficient
beta  = 0.0095; % Humans becoming zombies by being bitten coefficient
delta = 0.0001; % Natural human death rate coefficient
zeta  = 0.0001; % Dead Humans turn into zombies coefficient

rho   = 0.005;  % Rate of infected turning into zombies
kappa = 0.001;  % Rate of infected entering quarantine
sigma = 0.001;  % Rate of zombies  entering quarantine
gamma = 0.0001; % Rate of escape from quarantine, which are killed



% x = [S I Z R Q]
SIZRQ = @(t,x)[Pi - beta*x(1).*x(3) - delta*x(1);...
                    beta*x(1).*x(3) - (rho+delta+kappa)*x(2);...
                    rho*x(2) + zeta*x(4) - alpha*x(1).*x(3) - sigma*x(3);...
                    delta*(x(1) + x(2)) + alpha*x(1).*x(3) - zeta*x(4) + gamma*x(5);...
                    kappa*x(2) + sigma*x(3) - gamma*x(5)];
    
[t,x] = ode45(SIZRQ,[0,T],[S_0,I_0,Z_0,R_0,Q_0]);

plot(t,x(:,1)); hold on;
plot(t,x(:,2));
plot(t,x(:,3));
plot(t,x(:,4));
plot(t,x(:,5));
legend('Susceptible','Infected','Zombie','Removed','Quarantined','location','best')
xlabel('t');
hold on;




