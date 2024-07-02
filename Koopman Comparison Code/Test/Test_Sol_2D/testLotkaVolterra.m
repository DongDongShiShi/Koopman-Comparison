%Test VanderPol Osicaltor
clear; close all; clc;

addpath('2D/algorithm');
addpath('2D/model');


%% Koopman Spectral Linearization
%%parameter
N_values=5:2:13;
r=[3;3];
T=10;
tspan=linspace(0,T,1000);
op=1;
p=struct();
p.alpha=1.1;
p.beta=0.4;
p.delta=0.1;
p.gamma=0.4;
x0=[5;5];
dynamics = @LotkaVolterraDynamics;
K_errors=zeros(length(N_values),2);

% x1-t standard solution
figure
[t,xs] = ode45(@(t,x)LK(t,x,p),tspan,x0);
plot(t,xs(:,1),'LineWidth',2,'Color','k')
xlabel('t')
ylabel('x1')
hold on;

% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
 
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init1);
plot(t,y(:,(N_values(i)^2+1)/2),'--','LineWidth',2);
K_errors(i,1)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,1)));
hold on

end
legend('Real','N=5','N=7','N=9','N=11','N=13')
title("x1-t of LotkaVolterra Dynamics","FontSize",16, 'interpreter', 'latex')

%% x2-t standard solution
figure
plot(t,xs(:,2),'LineWidth',2,'Color','k')
xlabel('t')
ylabel('x2')
hold on;

% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
 
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init2);
plot(t,y(:,(N_values(i)^2+1)/2),'--',"LineWidth",2);
hold on
K_errors(i,2)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,2)));
end
legend('Real','N=5','N=7','N=9','N=11','N=13', 'interpreter', 'latex')
title("x2-t of LotkaVolterra Dynamics","FontSize",16, 'interpreter', 'latex')

%% Error-Truncation order N=5,7,9,11,13
figure;
semilogy(N_values,K_errors(:,1),'o-','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(N_values,K_errors(:,2),'^-','Color','k','LineWidth',2,'MarkerSize',10)
title("MAE error wrt truncation order","FontSize",16, 'interpreter', 'latex')
legend('x1 error','x2 error', 'interpreter', 'latex')
xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
xlim([4,14])
ylim([1e-3,1e0])
%% Test Error wrt radius r1-error of x1 and x2,N=13;
r1=linspace(2.5,5,10);
r2=3*ones(1,length(r1));
r_values=[r1;r2];
r_errors=zeros(length(r_values),2);
N=13;
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
xlim([2,6])
ylim([1e-3,1e0])
%% Test Error wrt radius r2 N=13
r2=linspace(2.2,3.2,10);
r1=3*ones(1,length(r2));
r_values=[r1;r2];
r_errors=zeros(length(r_values),2);
N=13;
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
xlim([2,3.5])
ylim([1e-3,1e1])
function dxdt = LK(t,x,p)
dxdt  = [p.alpha*x(1)-p.beta*x(1)*x(2);p.delta*x(1)*x(2)-p.gamma*x(2)];
end


