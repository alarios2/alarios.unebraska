% Predator Prey "Zombie" model
close all; clear all;

S_0 = 1; % Initial Suseptible (Human) Population (maybe in thousands)
Z_0 = 0.1; % Initial Zombie Population (maybe in thousands)

% % x(1) represents rabbits, x(2) represents foxes.
alpha = 4; 
beta = 2; 
delta = 3; 
gamma = 3;
SZ = @(t,x) [alpha*x(1) - beta*x(1)*x(2);
            -gamma*x(2) + delta*x(1)*x(2)];

[t,x] = ode45(SZ,[0,5],[S_0,Z_0]);

subplot(1,2,1)
plot(x(:,1),x(:,2));
axis([0,6,0,10]);
axis('square');
xlabel('Susceptible');
ylabel('Zombies');

subplot(1,2,2)
plot(t,x(:,1)); hold on;
plot(t,x(:,2)); 
legend('Susceptible','Zombie','location','best')
xlabel('t');