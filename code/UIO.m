clear all;close all;clc

A = [0.8 1 0;0 0.5 -1;0 0 0.7];
B = [0 0.0058 0 -0.0058;-0.0058 0 0.0058 0;-0.0002 0.0002 -0.0002 0.0002];
C = [1 -1 0;0 1 1;0 -1 1];
F1 = C;
D = 0;
E = [0.5;0.4;0.3];
f1=0.2;
f2=0.5;
iter=300;
num=1:iter;
t=0.01*num;
%%
Bi = cell(1,5);
Bi{1} = B; 
for i=1:4
    mid = B;
    mid(:,i)=0;
    Bi{i+1}=mid;
end
%% UIO_Healthy
% E=[Bi{1}(:,1) Bi{1}(:,2) E];
% B2=[0 0 0 -0.0058;0 0 0.0058 0;0 0 -0.0002 0.0002];
rank(C*E);
rank(E);
%H = E / (C * E);
H = E*inv((C*E)'*(C*E))*(C*E)';
T = eye(3)-H*C;
% TT=T*B2;
TT=T*Bi{1};
A1 = T*A;
rank(obsv(A1,C));

K1 = (place(A1',C',[0.5,0.8,0.2]))';
L = (place(A',C',[0.5,0.8,0.2]))';
F = A1-K1*C;  %F=N
K = K1 + F*H; %K=K1+K2

% E = [0.5;0.4;0.3];

%% UIO_F1
E_1=[Bi{1}(:,1) E];
rank(C*E_1);
rank(E_1);
H_1 = E_1 / (C * E_1);
%H = E*inv((C*E)'*(C*E))*(C*E)'
T_1 = eye(3)-H_1*C;
TT_1=T_1*Bi{2};
A_1 = T_1*A;
rank(obsv(A_1,C));

K_1 = (place(A_1',C',[0.5,0.8,0.2]))';
F_1 = A_1-K_1*C;  %F=N
KF_1 = K_1 + F_1*H_1; %K=K1+K2
%% UIO_F2
E_2=[Bi{1}(:,2) E];

rank(C*E_2);
rank(E_2);
H_2 = E_2 / (C * E_2);
%H = E*inv((C*E)'*(C*E))*(C*E)'
T_2 = eye(3)-H_2*C;
TT_2=T_2*Bi{3};
A_2 = T_2*A;
rank(obsv(A_2,C));

K_2 = (place(A_2',C',[0.5,0.8,0.2]))';
F_2 = A_2-K_2*C;  %F=N
KF_2 = K_2 + F_2*H_2; %K=K1+K2
%%
u=181.2171*ones(4,iter); %输入
W=normrnd(0,0.3,1,iter);%扰动
%原系统
x = zeros(3,iter);
x(:,1)=[10;1;-10];
y = zeros(3,iter);
%UIO_Healthy
z = zeros(3,iter);
x_hat = zeros(3,iter);
x_hat(:,1)=[10;1;-10];
y_hat = zeros(3,iter);
%UIO_F1
z_f1 = zeros(3,iter);
x_f1 = zeros(3,iter);
x_f1(:,1)=[10;1;-10];
y_f1 = zeros(3,iter);
%UIO_F2
z_f2 = zeros(3,iter);
x_f2 = zeros(3,iter);
x_f2(:,1)=[10;1;-10];
y_f2 = zeros(3,iter);
%%
for k=1:iter
    % fault occurs
    if (k>=100 && k<=200)
        B(:,1)=f1*B(:,1);
    end
    %System
    y(:,k) = C * x(:,k);
    x(:,k+1) = A * x(:,k) + B * u(:,k) + E * W(:,k);
    
    %Unknown Input Observer 0
    x_hat(:,k) = z(:,k) + H * y(:,k);
    y_hat(:,k) = C * x_hat(:,k);
    z(:,k+1) = F * z(:,k) + TT * u(:,k) + K * y(:,k);
    
    %Unknown Input Observer 1
    x_f1(:,k) = z_f1(:,k) + H_1 * y(:,k);
    y_f1(:,k) = C * x_f1(:,k);
    z_f1(:,k+1) = F_1 * z_f1(:,k) + TT_1 * u(:,k) + KF_1 * y(:,k);
    
    %Unknown Input Observer 2
    x_f2(:,k) = z_f2(:,k) + H_2 * y(:,k);
    y_f2(:,k) = C * x_f2(:,k);
    z_f2(:,k+1) = F_2 * z_f2(:,k) + TT_1 * u(:,k) + KF_2 * y(:,k);
    
    %Residual Error
    e_hat(:,k) = y(:,k) - y_hat(:,k);
    e_f1(:,k) = y(:,k) - y_f1(:,k);
    e_f2(:,k) = y(:,k) - y_f2(:,k);
end
%% Plot the results
t = linspace(0,0.01*iter,iter);
figure(1)
plot(t,2+x(1,1:iter),'b','LineWidth',2),hold on
plot(t,x(2,1:iter),'g','LineWidth',2),hold on
plot(t,-2+x(3,1:iter),'r','LineWidth',2),hold on

plot(t,2+x_hat(1,:),'xb'),hold on
plot(t,x_hat(2,:),'xg'),hold on
plot(t,-2+x_hat(3,:),'xr'),hold off
xlabel('Time'),ylabel('State Estimation')

h = legend('$ x_1 $','$ x_2 $','$ x_3 $','$ \hat{x}_1 $','$ \hat{x}_2 $','$ \hat{x}_3 $');
set(h,'Interpreter','latex');
title('UIO healthy')
%%
figure(2)
plot(t,2+x(1,1:iter),'b','LineWidth',2),hold on
plot(t,x(2,1:iter),'g','LineWidth',2),hold on
plot(t,-2+x(3,1:iter),'r','LineWidth',2),hold on

plot(t,2+x_f2(1,:),'xb'),hold on
plot(t,x_f2(2,:),'xg'),hold on
plot(t,-2+x_f2(3,:),'xr'),hold off
xlabel('Time'),ylabel('State Estimation')

h = legend('$ x_1 $','$ x_2 $','$ x_3 $','$ \hat{x}_1 $','$ \hat{x}_2 $','$ \hat{x}_3 $');
set(h,'Interpreter','latex');
title('UIO fault2')
%%
figure(3)
plot(t,2+x(1,1:iter),'b','LineWidth',2),hold on
plot(t,x(2,1:iter),'g','LineWidth',2),hold on
plot(t,-2+x(3,1:iter),'r','LineWidth',2),hold on

plot(t,2+x_f1(1,:),'xb'),hold on
plot(t,x_f1(2,:),'xg'),hold on
plot(t,-2+x_f1(3,:),'xr'),hold off

xlabel('Time'),ylabel('State Estimation')

h = legend('$ x_1 $','$ x_2 $','$ x_3 $','$ \hat{x}_1 $','$ \hat{x}_2 $','$ \hat{x}_3 $');
set(h,'Interpreter','latex');
title('UIO fault1')
%%
figure(4)
plot(t,e_hat(1,:)),hold on
plot(t,e_hat(2,:)),hold on
plot(t,e_hat(3,:)),hold off
title('UIO healthy')
ylabel('State Estimation Error')
%%
figure(5)
plot(t,e_f1(1,:)),hold on
plot(t,e_f1(2,:)),hold on
plot(t,e_f1(3,:)),hold off
title('UIO fault1')
ylabel('State Estimation Error')
%%
figure(6)
plot(t,e_f2(1,:)),hold on
plot(t,e_f2(2,:)),hold on
plot(t,e_f2(3,:)),hold off
title('UIO fault2')
ylabel('State Estimation Error')



