%Test Simple Pendulum
clear; close all; clc;

addpath('2D/algorithm');
addpath('2D/model');


%% Koopman Spectral Linearization
%%parameter
N_values=3:2:9;
r=[1;1];
T=10;
tspan=linspace(0,T,1000);
op=1;
p=struct();
p.g=9.8;
p.L=9.8;
x0=[0.1;0.1];
dynamics = @SimplePendulumDynamics;
K_errors=zeros(length(N_values),2);

% x1-t standard solution
figure
[t,xs] = ode45(@(t,x)SP(t,x,p),tspan,x0);
plot(t,xs(:,1),'LineWidth',2,'Color','k')
xlabel('t')
ylabel('x1')
hold on;

%lb = [x0(1)-r(1), x0(2)-r(2)];
%ub = [x0(1)+r(1), x0(2)+r(2)];

% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
 
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init1);
plot(t,y(:,(N_values(i)^2+1)/2),'--','LineWidth',2);
K_errors(i,1)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,1)));
hold on
end

legend('Real','N=3','N=5','N=7','N=9','N=11')
title("x1-t of Simple Pendulum Dynamics","FontSize",16, 'interpreter', 'latex')

%% x2-t standard solution
figure
[t,xs] = ode45(@(t,x)SP(t,x,p),tspan,x0);
plot(t,xs(:,2),'LineWidth',2,'Color','k')
xlabel('t')
ylabel('x2')
hold on;

lb = [x0(1)-r(1), x0(2)-r(2)];
ub = [x0(1)+r(1), x0(2)+r(2)];

% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
 
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init2);
plot(t,y(:,(N_values(i)^2+1)/2),'--',"LineWidth",2);
hold on
K_errors(i,2)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,2)));
end
legend('Real','N=3','N=5','N=7','N=9','N=11', 'interpreter', 'latex')
title("x2-t of Simple Pendulum Dynamics","FontSize",16, 'interpreter', 'latex')

%% Error-Truncation order N=5,7,9,11,13
figure;
semilogy(N_values,K_errors(:,1),'o-','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(N_values,K_errors(:,2),'^-','Color','k','LineWidth',2,'MarkerSize',10)
title("MAE error wrt truncation order","FontSize",16, 'interpreter', 'latex')
legend('x1 error','x2 error', 'interpreter', 'latex')

xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
xlim([2,10])
ylim([1e-5,1e-1])

%% Test Error wrt radius r1-error of x1 and x2,N=13;
r1=linspace(0.05,1.5,10);
r2=1*ones(1,length(r1));
r_values=[r1;r2];
r_errors=zeros(length(r_values),2);
N=9;
for i = 1:length(r_values)
    [K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N, x0, r_values(:,i),op);
    [t,y1]=ode45(@(t,y)K*y,tspan,init1);
    [t,y2]=ode45(@(t,y)K*y,tspan,init2);
    r_errors(i,1)=mean(abs(y1(:,(N^2+1)/2)-xs(:,1)));
    r_errors(i,2)=mean(abs(y2(:,(N^2+1)/2)-xs(:,2)));
end
figure;
semilogy(r_values(1,:),r_errors(:,1),'-o','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(1,:),r_errors(:,2),'-^','Color','k','LineWidth',2,'MarkerSize',10)
legend('error of x1','error of x2', 'interpreter', 'latex')
title("MAE-radius r1", 'interpreter', 'latex')
xlabel("radius r","FontSize",16,'interpreter', 'latex')
ylabel("MAE","FontSize",16,'interpreter', 'latex')
xlim([0,1.6])
ylim([1e-5,1e-1])

%% Test Error wrt radius r2 N=9
r2=linspace(0.5,1.5,10);
r1=1*ones(1,length(r2));
r_values=[r1;r2];
r_errors=zeros(length(r_values),2);
N=9;
for i = 1:length(r_values)
    [K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N, x0, r_values(:,i),op);
    [t,y1]=ode45(@(t,y)K*y,tspan,init1);
    [t,y2]=ode45(@(t,y)K*y,tspan,init2);
    r_errors(i,1)=mean(abs(y1(:,(N^2+1)/2)-xs(:,1)));
    r_errors(i,2)=mean(abs(y2(:,(N^2+1)/2)-xs(:,2)));
end
figure;
semilogy(r_values(2,:),r_errors(:,1),'-o','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(2,:),r_errors(:,2),'-^','Color','k','LineWidth',2,'MarkerSize',10)
legend('error of x1','error of x2', 'interpreter', 'latex')
title("MAE-radius r2", 'interpreter', 'latex')
xlabel("radius r","FontSize",16,'interpreter', 'latex')
ylabel("MAE","FontSize",16,'interpreter', 'latex')
xlim([0.4,1.6])
ylim([1e-5,1e-1])

function dxdt = SP(t,x,p)
dxdt  = [x(2);- (p.g / p.L) * sin(x(1))];
end
