%Test VanderPol Osicaltor
clear; close all; clc;

addpath('funcs');
addpath('models');


%% Koopman Spectral Linearization
%%parameter
N_values=[3,5];
r=[0.5;0.5];
T=5;
op=1;
p=0.1;
x0=[0.1;0.1];
dynamics = @VanDelPolDynamics;

% x1-t standard solution
[t,x] = ode45(@vdp,[0,T],x0);
figure;
plot(t,x(:,1),'-','LineWidth',2,'Color','k')
hold on;


% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
 
% Compute solution
[t,y]=ode45(@(t,y)K*y,[0,T],init1);
plot(t,y(:,(N_values(i)^2+1)/2),'--','LineWidth',2);
hold on
legend('Real','N=3','N=5','N=7','N=9')
end


% x2-t standard solution
[t,x] = ode45(@vdp,[0,T],x0);
figure;
plot(t,x(:,2),'-','LineWidth',2,'Color','k')
hold on;


% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
 
% Compute solution
[t,y]=ode45(@(t,y)K*y,[0,T],init2);
plot(t,y(:,(N_values(i)^2+1)/2),'--','LineWidth',2);
hold on
legend('Real','N=3','N=5','N=7','N=9')
end


function dxdt = vdp(t,x)
dxdt  = [x(2);-x(1)+0.1*(1-x(1)^2)*x(2)];
end