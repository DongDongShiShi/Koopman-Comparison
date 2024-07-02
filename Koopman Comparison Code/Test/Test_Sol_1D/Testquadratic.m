%% Test Quadratic
clear; close all; clc;

addpath('algorithm');
addpath('model');

%% Parameters
T=10;
x0=0.08;
p=struct();
p.k=1;
dynamics=@QuadraticDynamics;
N_values=[3,5,7,9,11];
tspan=linspace(0,T,1000);
r=0.03;
K_errors=zeros(length(N_values),1);
C_errors=zeros(length(N_values),1);


%% Carleman Linearization
%Standard Solution
%[t,xs] = ode45(@(t,x)dynamics(x,p),tspan,x0);
%figure;
%plot(t,xs,'LineWidth',2,'Color','k');
%hold on;

%for i=1:length(N_values)

%Construct Carleman Matrix and Initial Condition
%powers=1:N_values(i);
%y0=x0.^powers;
%v=1:N_values(i)-1;
%C=diag(v,1);

%Carleman Figure
%[t,y]=ode45(@(t,y)C*y,tspan,y0);
%plot(t,y(:,1),'--','LineWidth',2)
%hold on
%C_errors(i)=mean((y(:,1)-xs).^2);
%end
          
%legend('Real Values', 'N=3', 'N=5', 'N=7', 'N=9','N=11');
%xlabel('t');
%ylabel('x');
%title('Carleman Linearization');
%hold off



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
legend('Real Values', 'N=3', 'N=5', 'N=7', 'N=9','N=11');
xlabel('t','interpreter', 'latex');
ylabel('x','interpreter', 'latex');
title('Koopman Spectral Linearization','FontSize',16);




%% Figure of Error wrt Truncation Order r=0.03, N=3,5,7,9,11
figure;
semilogy(N_values,K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
xlim([0,13])
ylim([1e-6,1e-1])
xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')


%% Test Radius N=11, r=0.01-0.05
N=11;
tspan=linspace(0,T,1000);
r_values=linspace(0.01,0.05,20);
K_errors=zeros(length(r_values),1);

for i = 1:length(r_values)
[K,y0] = KSLinearization(dynamics,p,r_values(i),x0,N);
[t,y]=ode45(@(t,y)K*y,tspan,y0);
K_errors(i)=mean(abs(y(:,(N+1)/2)-xs));
end

figure;
semilogy(r_values,K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
xlim([0,0.06])
ylim([1e-6,1e-2])
xlabel("radius r","FontSize",16,'interpreter', 'latex')
ylabel("MAE","FontSize",16,'interpreter', 'latex')




