%% 3D: Lorenz attractor
% beta = 8/3;sigma = 10;rho = 28;
% lorenz = @(t,x)[sigma*(x(2)-x(1));x(1).*(rho-x(3))-x(2);x(1).*x(2)-beta*x(3)];
% [t,x] = ode45(lorenz,[0,300],[0,1,1]);
% plot3(x(:,1),x(:,2),x(:,3));

close all;

S_0 = 500; % Initial Suseptible (Human) Population
Z_0 = 0; % Initial Zombie Population
R_0 = 0; % Initial Removed Population (Dead)

T = 50; % Final time

Pi = 0.0001;          % Human birth rate (maybe proportional to delta)
alpha = 0.0075;  % Humans killing zombies coefficient
beta = 0.0055;   % Humans becoming zombies by being bitten coefficient
delta = 0.0001;  % Natural human death rate coefficient
zeta = 0.09;     % Dead Humans turn into zombies coefficient




% x = [S Z R]
SZR = @(t,x)[Pi - beta*x(1).*x(2) - delta*x(1);...
             beta*x(1).*x(2) + zeta*x(3) - alpha*x(1).*x(2);...
             delta*x(1) + alpha*x(1).*x(2) - zeta*x(3)];
    
[t,x] = ode45(SZR,[0,T],[S_0,Z_0,R_0]);

% subplot(1,3,1);
plot(t,x(:,1)); hold on;
plot(t,x(:,2)); hold on;
plot(t,x(:,3));
legend('Susceptible','Zombie','Removed','location','best')
xlabel('t');
hold on;

% % subplot(1,3,2);
% 
% title('Zombie');
% xlabel('t');
% hold on;
% 
% % subplot(1,3,3);
% plot(t,x(:,3));
% title('Removed');
% xlabel('t');





