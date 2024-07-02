%Test Kraichnan-Orszag
clear; close all; clc;

addpath('algorithm');
addpath('models');


%% Koopman Spectral Linearization
%%parameter
N_values=[3,5,7,9];
r=[0.1;0.1;0.1];
T=5;
op=1;
p=struct();
p.k1=1;
p.k2=1;
p.k3=1;
x0=[0.1;-0.2;0.3];
tspan=linspace(0,T,100);
dynamics = @KraichnanOrszagDynamics;
K_errors=zeros(length(N_values),3);


%%standard solution
[t,xs] = ode45(@(t,x)KO(t,x,p),tspan,x0);
% x1-t standard solution
figure
plot(t,xs(:,1),'LineWidth',2,'Color','k')
xlabel('t')
ylabel('x1')
hold on;


 
% Compute Koopman Spectral Linearization Matrix

 
for i=1:length(N_values)
[K,init1,init2,init3] = KoopmanLinearization_3D(dynamics, p, T, N_values(i), x0, r,op);
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init1);
plot(t,y(:,(N_values(i)^3+1)/2),'--','LineWidth',2);
K_errors(i,1)=mean(abs(y(:,(N_values(i)^3+1)/2)-xs(:,1)));
hold on

end
legend('Real','N=3','N=5','N=7','N=9','N=11')
title("x1-t of KraichnanOrszag Dynamics")

% x2-t standard solution
figure
plot(t,xs(:,2),'LineWidth',2,'Color','k')
xlabel('t')
ylabel('x2')
hold on;

% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2,init3] = KoopmanLinearization_3D(dynamics, p, T, N_values(i), x0, r,op);
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init2);
plot(t,y(:,(N_values(i)^3+1)/2),'--','LineWidth',2);
K_errors(i,2)=mean(abs(y(:,(N_values(i)^3+1)/2)-xs(:,2)));
hold on
end
legend('Real','N=3','N=5','N=7','N=9')
title("x2-t of KraichnanOrszag Dynamics")

% x3-t standard solution
figure
plot(t,xs(:,3),'LineWidth',2,'Color','k')
xlabel('t')
ylabel('x3')
hold on;

lb = [x0(1)-r(1), x0(2)-r(2), x0(3)-r(3)];
ub = [x0(1)+r(1), x0(2)+r(2), x0(3)+r(3)];

% Compute Koopman Spectral Linearization Matrix
for i=1:length(N_values)
[K,init1,init2,init3] = KoopmanLinearization_3D(dynamics, p, T, N_values(i), x0, r,op);
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init3);
plot(t,y(:,(N_values(i)^3+1)/2),'--','LineWidth',2);
K_errors(i,3)=mean(abs(y(:,(N_values(i)^3+1)/2)-xs(:,3)));
hold on
end
legend('Real','N=3','N=5','N=7','N=9','N=11')
title("x3-t of KraichnanOrszag Dynamics")

%% Error-N with N=3,5,7,9 r=[0.1,0.1,0.1] T=5
figure;
semilogy(N_values,K_errors(:,1),'o-','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(N_values,K_errors(:,2),'^-','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(N_values,K_errors(:,3),'s-','Color','k','LineWidth',2,'MarkerSize',10)
title("MAE error wrt truncation order","FontSize",16, 'interpreter', 'latex')
legend('x1 error','x2 error','x3 error ','interpreter', 'latex')
xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
xlim([2,10]);ylim([1e-8,1e-2])



%% Error-r1 N=9 r1=0.05-0.15
r1=linspace(0.05,0.15,10);
r2=0.1*ones(1,length(r1));
r3=0.1*ones(1,length(r1));
r_values=[r1;r2;r3];
r_errors=zeros(length(r_values),3);
N=9;
for i = 1:length(r_values)
    [K,init1,init2,init3] = KoopmanLinearization_3D(dynamics, p, T, N, x0, r_values(:,i),op);
    [t,y1]=ode45(@(t,y)K*y,tspan,init1);
    [t,y2]=ode45(@(t,y)K*y,tspan,init2);
    [t,y3]=ode45(@(t,y)K*y,tspan,init3);
    r_errors(i,1)=mean(abs(y1(:,(N^3+1)/2)-xs(:,1)));
    r_errors(i,2)=mean(abs(y2(:,(N^3+1)/2)-xs(:,2)));
    r_errors(i,3)=mean(abs(y3(:,(N^3+1)/2)-xs(:,3)));
end
figure;
semilogy(r_values(1,:),r_errors(:,1),'-o','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(1,:),r_errors(:,2),'-^','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(1,:),r_errors(:,3),'-s','Color','k','LineWidth',2,'MarkerSize',10)
legend('error of x1','error of x2','error of x3', 'interpreter', 'latex')
title("MAE-radius r1", 'interpreter', 'latex')
xlabel("radius r","FontSize",16,'interpreter', 'latex')
ylabel("MAE","FontSize",16,'interpreter', 'latex')
xlim([0.04,0.16])
ylim([1e-8,1e-2])

%% Error-r2 N=9 r1=0.05-0.15
r2=linspace(0.02,0.22,10);
r1=0.1*ones(1,length(r2));
r3=0.1*ones(1,length(r2));
r_values=[r1;r2;r3];
r_errors=zeros(length(r_values),3);
N=9;
for i = 1:length(r_values)
    [K,init1,init2,init3] = KoopmanLinearization_3D(dynamics, p, T, N, x0, r_values(:,i),op);
    [t,y1]=ode45(@(t,y)K*y,tspan,init1);
    [t,y2]=ode45(@(t,y)K*y,tspan,init2);
    [t,y3]=ode45(@(t,y)K*y,tspan,init3);
    r_errors(i,1)=mean(abs(y1(:,(N^3+1)/2)-xs(:,1)));
    r_errors(i,2)=mean(abs(y2(:,(N^3+1)/2)-xs(:,2)));
    r_errors(i,3)=mean(abs(y3(:,(N^3+1)/2)-xs(:,3)));
end
figure;
semilogy(r_values(2,:),r_errors(:,1),'-o','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(2,:),r_errors(:,2),'-^','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(2,:),r_errors(:,3),'-s','Color','k','LineWidth',2,'MarkerSize',10)
legend('error of x1','error of x2','error of x3', 'interpreter', 'latex')
title("MAE-radius r2", 'interpreter', 'latex')
xlabel("radius r","FontSize",16,'interpreter', 'latex')
ylabel("MAE","FontSize",16,'interpreter', 'latex')
xlim([0,0.25])
ylim([1e-8,1e-2])
%% Error-r2 N=9 r1=0.05-0.15
r3=linspace(0.02,0.22,10);
r1=0.1*ones(1,length(r3));
r2=0.1*ones(1,length(r3));
r_values=[r1;r2;r3];
r_errors=zeros(length(r_values),3);
N=9;
for i = 1:length(r_values)
    [K,init1,init2,init3] = KoopmanLinearization_3D(dynamics, p, T, N, x0, r_values(:,i),op);
    [t,y1]=ode45(@(t,y)K*y,tspan,init1);
    [t,y2]=ode45(@(t,y)K*y,tspan,init2);
    [t,y3]=ode45(@(t,y)K*y,tspan,init3);
    r_errors(i,1)=mean(abs(y1(:,(N^3+1)/2)-xs(:,1)));
    r_errors(i,2)=mean(abs(y2(:,(N^3+1)/2)-xs(:,2)));
    r_errors(i,3)=mean(abs(y3(:,(N^3+1)/2)-xs(:,3)));
end
figure;
semilogy(r_values(3,:),r_errors(:,1),'-o','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(3,:),r_errors(:,2),'-^','Color','k','LineWidth',2,'MarkerSize',10)
hold on
semilogy(r_values(3,:),r_errors(:,3),'-s','Color','k','LineWidth',2,'MarkerSize',10)
legend('error of x1','error of x2','error of x3', 'interpreter', 'latex')
title("MAE-radius r3", 'interpreter', 'latex')
xlabel("radius r","FontSize",16,'interpreter', 'latex')
ylabel("MAE","FontSize",16,'interpreter', 'latex')
xlim([0,0.22])
ylim([1e-8,1e-2])
function dxdt = KO(~,x,p)
dxdt  = [p.k1*x(2)*x(3);p.k2*x(1)*x(3);p.k3*-2*x(1)*x(2)];
end