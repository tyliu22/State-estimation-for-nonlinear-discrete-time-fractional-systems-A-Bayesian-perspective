%*************************************************************************%
%   分数阶粒子滤波仿真复现
%   论文：     fractional order PF  性能分析 RMSE & ERROR
%   目的：分数阶粒子滤波算法 RMSE测试
%         对系统噪声均值进行估计
%         函数实验:    D^{0.7} x_k = 3*sin(2*x_{k-1}) - exp(x_{k-1}) + w_k
%                              y_k = x_k + v_k
%   结果：较好的对状态进行估计
%
%   备注：分数阶粒子滤波的算法
%           RMSE测试       ERROR测试
%*************************************************************************%



clc;
clear all;

LineWidth = 1.5;
SimuTimes = 50;    % 仿真时长


  
%--------------ERROR--------------%

load('x_SampleParticle.mat') % x_SampleParticle
load('X_RealState.mat')      % X_RealState
load('x_EstiState.mat')      % x_EstiState

LineWidth = 1.5;
SimuTimes = 50;    % 仿真时长
NumParticle = 100;  % 粒子个数

%% estimation accuracy
t = 1 : SimuTimes;
figure;
plot(t, X_RealState,'b',t, x_EstiState, 'r--','linewidth',LineWidth);
hold on

for i = 1:NumParticle-1
    plot(t, x_SampleParticle(i,:),'c:');
    hold on
end
plot(t, x_SampleParticle(NumParticle,:),'c');
hold on
plot(t, X_RealState,'b',t, x_EstiState, 'r--','linewidth',LineWidth);
legend('real state','estimated state','particles');
%L2 = plot(t, x_SampleParticle(100,k),'y',t, X_RealState,'b',t, x_EstiState, 'r--','linewidth',LineWidth);
%Esitimated_state = legend('real state','estimated state','Location','best');
%set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 50 -10 2]) % 设置坐标轴在指定的区间
axis normal
set(gca,'FontSize',10); 
xlabel('iteration times','FontSize',7); 
ylabel('estimated state','FontSize',7);
% 设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)

%--------------ERROR--------------%
Error_State = abs(X_RealState-x_EstiState);
figure;
plot(t, Error_State,'b','linewidth',LineWidth);
Esitimated_state = legend('error','Location','best');
set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % 设置坐标轴在指定的区间
axis normal
set(gca,'FontSize',10); 
xlabel('$k$','FontSize',7,'Interpreter','latex') 
ylabel('error','FontSize',7);
% 设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)

%--------------RELATIVE ERROR--------------%
%-----------------相对误差-----------------%
Error_State = abs(X_RealState-x_EstiState)./abs(X_RealState);
figure;
plot(t, Error_State,'b','linewidth',LineWidth);
Esitimated_state = legend('relative error','Location','best');
set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % 设置坐标轴在指定的区间
axis normal
set(gca,'FontSize',10); 
xlabel('$k$','FontSize',7,'Interpreter','latex') 
ylabel('Relative Error','FontSize',7);
% 设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)



%--------------RELATIVE ERROR--------------%
%-----------------相对误差-----------------%
Error_State = abs(X_RealState-x_EstiState)./(abs(X_RealState)+abs(x_EstiState));
figure;
plot(t, Error_State,'b','linewidth',LineWidth);
Esitimated_state = legend('$\sigma$','Location','best');
set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % 设置坐标轴在指定的区间
axis normal
set(gca,'FontSize',10); 
xlabel('$k$','FontSize',7,'Interpreter','latex') 
ylabel('$\sigma$','FontSize',7);
% 设置坐标轴刻度字体名称，大小
set(gca,'FontName','Helvetica','FontSize',8)

