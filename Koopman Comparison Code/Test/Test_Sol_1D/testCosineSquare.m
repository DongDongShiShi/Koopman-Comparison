%% Test Cosinequare
clear; close all; clc;

addpath('algorithm');
addpath('model');

%% Parameters
T=10;
x0=0.9;
p=struct();
p.k=1;
dynamics=@CosineSquareDynamics;
N_values=3:2:11;
tspan=linspace(0,T,1000);
r=0.3;
K_errors=zeros(length(N_values),1);
C_errors=zeros(length(N_values),1);



%% Koopman Spectral Linearization 
[t,xs] = ode45(@(t,x)dynamics(x,p),tspan,x0);
figure;
plot(t,xs,'LineWidth',2,'Color','k');
hold on;

for i=1:length(N_values)


% Koopman Spectral Linearization Matirx
[K,y0] = KSLinearization(dynamics,p,r,x0,N_values(i));

% Koopman Spectral Figure
[t,y]=ode45(@(t,y)K*y,tspan,y0);
plot(t,y(:,(N_values(i)+1)/2),'--','LineWidth',2);
hold on
K_errors(i)=mean(abs(y(:,(N_values(i)+1)/2)-xs));
end
legend('Real Values', 'N=3', 'N=5', 'N=7', 'N=9');
xlabel('t','interpreter', 'latex');
ylabel('x','interpreter', 'latex');
title('Koopman Spectral Linearization','FontSize',16);




%% Figure of Error wrt Truncation Order r=0.3, N=3,5,7,9,11
figure;
semilogy(N_values,K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
xlim([2,12])
ylim([1e-4,1e-1])
xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')


%% Test Radius N=11, r=0.3-0.8
N=9;
tspan=linspace(0,T,1000);
r_values=linspace(0.3,0.8,10);
K_errors=zeros(length(r_values),1);

for i = 1:length(r_values)
[K,y0] = KSLinearization(dynamics,p,r_values(i),x0,N);
[t,y]=ode45(@(t,y)K*y,tspan,y0);
K_errors(i)=mean(abs(y(:,(N+1)/2)-xs));
end

figure;
semilogy(r_values,K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
xlim([0.2,1])
ylim([1e-5,1e-2])
xlabel("radius r","FontSize",16,'interpreter', 'latex')
ylabel("MAE","FontSize",16,'interpreter', 'latex')


