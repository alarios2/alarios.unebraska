% Zombie model in basic S-Z-R form.
% Model is from 2009 paper of Munz, Hudea, Imad, and Smith 
% published in Infectious Disease Modeling Research Papers.
close all;

S_0 = 500; % Initial Suseptible (Human) Population
Z_0 = 0; % Initial Zombie Population
R_0 = 0; % Initial Removed Population (Dead)

T = 15; % Final time

% The paper does not give all the parameters, so results differ silightly.
Pi    = 0.0001; % Human birth rate (maybe proportional to delta)
alpha = 0.005;  % Humans killing zombies coefficient
beta  = 0.0095; % Humans becoming zombies by being bitten coefficient
delta = 0.0001; % Natural human death rate coefficient
zeta  = 0.0001; % Dead Humans turn into zombies coefficient

% x = [S Z R]
SZR = @(t,x)[Pi - beta*x(1).*x(2) - delta*x(1);...
             beta*x(1).*x(2) + zeta*x(3) - alpha*x(1).*x(2);...
             delta*x(1) + alpha*x(1).*x(2) - zeta*x(3)];
    
[t,x] = ode45(SZR,[0,T],[S_0,Z_0,R_0]);

% subplot(1,3,1);
plot(t,x(:,1)); hold on;
plot(t,x(:,2));
plot(t,x(:,3));
legend('Susceptible','Zombie','Removed','location','best')
xlabel('t');