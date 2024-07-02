%% Error-N Compared with Carleman Linearization
clear;
clc;
%% Quadratic 
%%Parameters
T=10;
x0=0.08;
p=struct();
p.k=1;
dynamics=@QuadraticDynamics;
N_values=[3,5,7,9,11];
tspan=linspace(0,T,1000);
r=0.03;
K_errors=zeros(length(N_values),1);C_errors=zeros(length(N_values),1);
K_time=zeros(length(N_values),1);C_time=zeros(length(N_values),1);

%%Carleman Linearization
%Standard Solution
[t,xs] = ode45(@(t,x)dynamics(x,p),tspan,x0);

for i=1:length(N_values)
%Construct Carleman Matrix and Initial Condition
tic;
powers=1:N_values(i);
y0=x0.^powers;
C=build_matrix(N_values(i),0,1);

%Carleman error
[t,y]=ode45(@(t,y)C*y,tspan,y0);
C_errors(i)=mean(abs(y(:,1)-xs));
C_time(i)=toc;
end

for i=1:length(N_values)
tic;
%%Koopman Spectral Linearization Matirx
[K,y0] = KSLinearization(dynamics,p,r,x0,N_values(i));

% Koopman Spectral error
[t,y]=ode45(@(t,y)K*y,tspan,y0);
K_errors(i)=mean(abs(y(:,(N_values(i)+1)/2)-xs));
K_time(i) = toc;
end

%%Figure of Error wrt Truncation Order r=0.03, N=3,5,7,9,11
figure;
semilogy(N_values,K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
hold on
semilogy(N_values,C_errors,'-s','LineWidth',2,'Color','k','MarkerSize',10)
xlim([0,13])
ylim([1e-6,1e-1])
xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
legend("Error of Koopman","Error of Carleman")


%% LotkaVolterra
%%parameter
N_values=5:2:11;
r=[5;5];
T=1;
tspan=linspace(0,T,1000);
op=1;
p=struct();
p.alpha=1.1;
p.beta=0.4;
p.delta=0.1;
p.gamma=0.4;
x0=[5;5];
dynamics = @LotkaVolterraDynamics;
K_errors=zeros(length(N_values),2);C_errors=zeros(length(N_values),2);
K_time=zeros(length(N_values),1);C_time=zeros(length(N_values),1);
%%standard solution
[t,xs] = ode45(@(t,x)LK(t,x,p),tspan,x0);

%Construct Carleman Matrix and Initial Condition
for i=1:length(N_values)
tic;
powers=1:N_values(i);
y0=tensor_power_sequence(x0,N_values(i));
F1=[1.1,0;0,-0.4];
F2=[0,-0.4,0,0;0,0.1,0,0];
C=build_matrix(N_values(i),F1,F2);

%Carleman error
[t,y]=ode45(@(t,y)C*y,tspan,y0);
C_errors(i,1)=mean(abs(y(:,1)-xs(:,1)));
C_errors(i,2)=mean(abs(y(:,2)-xs(:,2)));
C_time(i)=toc;
end
%Koopman Spectral
for i=1:length(N_values)
tic;
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init1);
K_errors(i,1)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,1)));
[t,y]=ode45(@(t,y)K*y,tspan,init2);
K_errors(i,2)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,2)));
K_time(i)=toc;
end
total_C_errors=sqrt(mean(C_errors.^2, 2));
total_K_errors=sqrt(mean(K_errors.^2, 2));
%%Figure of error-N
figure;
semilogy(N_values,total_K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
hold on
semilogy(N_values,total_C_errors,'-s','LineWidth',2,'Color','k','MarkerSize',10)

xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
legend("Error of Koopman","Error of Carleman")
xlim([3,13])
ylim([1e-8,1e1])



%% KrainchnanOrszag
%%parameter
N_values=3:2:9;
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
K_errors=zeros(length(N_values),3);C_errors=zeros(length(N_values),3);
K_time=zeros(length(N_values),1);C_time=zeros(length(N_values),1);
%%standard solution
[t,xs] = ode45(@(t,x)KO(t,x,p),tspan,x0);

%Construct Carleman Matrix and Initial Condition
for i=1:length(N_values)
tic;
powers=1:N_values(i);
y0=tensor_power_sequence(x0,N_values(i));
F1=zeros(3);
F2=[0,0,0,0,0,0,0,1,0;0,0,1,0,0,0,0,0,0;0,-2,0,0,0,0,0,0,0];
C=build_matrix(N_values(i),F1,F2);

%Carleman error
[t,y]=ode45(@(t,y)C*y,tspan,y0);
C_errors(i,1)=mean(abs(y(:,1)-xs(:,1)));
C_errors(i,2)=mean(abs(y(:,2)-xs(:,2)));
C_errors(i,3)=mean(abs(y(:,3)-xs(:,3)));
C_time(i)=toc;
end

%Koopman Error
for i=1:length(N_values)
tic;
[K,init1,init2,init3] = KoopmanLinearization_3D(dynamics, p, T, N_values(i), x0, r,op);
[t,y]=ode45(@(t,y)K*y,tspan,init1);
K_errors(i,1)=mean(abs(y(:,(N_values(i)^3+1)/2)-xs(:,1)));
[t,y]=ode45(@(t,y)K*y,tspan,init2);
K_errors(i,2)=mean(abs(y(:,(N_values(i)^3+1)/2)-xs(:,2)));
[t,y]=ode45(@(t,y)K*y,tspan,init3);
K_errors(i,3)=mean(abs(y(:,(N_values(i)^3+1)/2)-xs(:,3)));
K_time(i)=toc;
end
total_C_errors=sqrt(mean(C_errors.^2, 2));
total_K_errors=sqrt(mean(K_errors.^2, 2));


%figure of error-N

figure;
semilogy(N_values,total_K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
hold on
semilogy(N_values,total_C_errors,'-s','LineWidth',2,'Color','k','MarkerSize',10)

xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
legend("Error of Koopman","Error of Carleman")
xlim([2,11])
ylim([1e-8,1e1])



%% CosineSquare
%Parameters
T=10;
x0=0.9;
p=struct();
p.k=1;
dynamics=@CosineSquareDynamics;
N_values=3:2:11;
tspan=linspace(0,T,1000);
r=0.3;
K_errors=zeros(length(N_values),1);C_errors=zeros(length(N_values),1);
K_time=zeros(length(N_values),1);C_time=zeros(length(N_values),1);

% standard solution
[t,xs] = ode45(@(t,x)dynamics(x,p),tspan,x0);


% Carleman Error
cof=TaylorCof(@(x)cos(x)^2, 0.9, 12);

for i=1:length(N_values)
tic;
%Construct Carleman Matrix and Initial Condition
powers=1:N_values(i);
y0=x0.^powers;
C=build_matrix(N_values(i),cof{2:end});
b=zeros(length(y0),1);
b(1)=cof{1};
%Carleman error
[t,y]=ode45(@(t,y)C*y+b,tspan,y0);
C_errors(i)=mean(abs(y(:,1)-xs));
C_time(i)=toc;
end

%Koopman Error
% Koopman Spectral Linearization Matirx
for i=1:length(N_values)
tic;
[K,y0] = KSLinearization(dynamics,p,r,x0,N_values(i));
[t,y]=ode45(@(t,y)K*y,tspan,y0);
K_errors(i)=mean(abs(y(:,(N_values(i)+1)/2)-xs));
K_time(i)=toc;
end

%figure of error-N
figure;
semilogy(N_values,K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
hold on
semilogy(N_values,C_errors,'-s','LineWidth',2,'Color','k','MarkerSize',10)
xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
legend("Error of Koopman","Error of Carleman")
xlim([2,13])
ylim([1e-5,1e2])



%% Simple Pendulum
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
K_errors=zeros(length(N_values),2);C_errors=zeros(length(N_values),2);
K_time=zeros(length(N_values),1);C_time=zeros(length(N_values),1);

% standard solution
[t,xs] = ode45(@(t,x)SP(t,x,p),tspan,x0);

% Carleman Error
for i = 1:length(N_values)

% Construct Carleman Matrix and Initial Condition  
cof = TaylorCof(@(x)sin(x), x0(2), 9);
matrix_array = cell(9);
matrix_array{1} = [0,1;-1,0];

for j = 2:9
    matrix_array{j} = sparse(2,2^j);
    matrix_array{j}(2,1) = cof{j};
end
tic;
C  = build_matrix(N_values(i),matrix_array{:});
y0 = tensor_power_sequence(x0,N_values(i));
[t,y]=ode45(@(t,y)C*y,tspan,y0);
C_errors(i,1)=mean(abs(y(:,1)-xs(:,1)));
C_errors(i,2)=mean(abs(y(:,2)-xs(:,2)));
C_time(i)=toc;
end

%Koopman Error
for i=1:length(N_values)
tic;
[K,init1,init2] = KoopmanLinearization_2D(dynamics, p, T, N_values(i), x0, r,op);
% Compute solution
[t,y]=ode45(@(t,y)K*y,tspan,init1);
K_errors(i,1)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,1)));
[t,y]=ode45(@(t,y)K*y,tspan,init2);
K_errors(i,2)=mean(abs(y(:,(N_values(i)^2+1)/2)-xs(:,2)));
K_time(i)=toc;
end
total_C_errors=sqrt(mean(C_errors.^2, 2));
total_K_errors=sqrt(mean(K_errors.^2, 2));

%%figure of error-N
figure;
semilogy(N_values,total_K_errors,'-o','LineWidth',2,'Color','k','MarkerSize',10)
hold on
semilogy(N_values,total_C_errors,'-s','LineWidth',2,'Color','k','MarkerSize',10)

xlabel('Truncation order N','FontSize',16,'interpreter', 'latex')
ylabel('MAE','FontSize',16,'interpreter', 'latex')
legend("Error of Koopman","Error of Carleman")
xlim([1,10])
ylim([1e-5,1e-1])




function dxdt = KO(~,x,p)
dxdt  = [p.k1*x(2)*x(3);p.k2*x(1)*x(3);p.k3*-2*x(1)*x(2)];
end
function dxdt = LK(~,x,p)
dxdt  = [p.alpha*x(1)-p.beta*x(1)*x(2);p.delta*x(1)*x(2)-p.gamma*x(2)];
end
function dxdt = SP(t,x,p)
dxdt  = [x(2);- (p.g / p.L) * sin(x(1))];
end
function result = tensor_power_sequence(x0, N)
    result = x0;
    for k = 2:N
        result = [result; kron(result(end-length(x0)^(k-1)+1:end),x0)];
    end
end