%*************************************************************************%
%   �����������˲����渴��
%   ���ģ�     fractional order PF  ���ܷ��� RMSE & ERROR
%   Ŀ�ģ������������˲��㷨 RMSE����
%         ��ϵͳ������ֵ���й���
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) - exp(x_{k-1}) + w_k
%                              y_k = x_k + v_k
%   ������ϺõĶ�״̬���й���
%
%   ��ע�������������˲����㷨
%           RMSE����       ERROR����
%*************************************************************************%



clc;
clear all;

LineWidth = 1.5;
SimuTimes = 50;    % ����ʱ��


  
%--------------ERROR--------------%

load('x_SampleParticle.mat') % x_SampleParticle
load('X_RealState.mat')      % X_RealState
load('x_EstiState.mat')      % x_EstiState

LineWidth = 1.5;
SimuTimes = 50;    % ����ʱ��
NumParticle = 100;  % ���Ӹ���

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
%axis([0 50 -10 2]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
xlabel('iteration times','FontSize',7); 
ylabel('estimated state','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)

%--------------ERROR--------------%
Error_State = abs(X_RealState-x_EstiState);
figure;
plot(t, Error_State,'b','linewidth',LineWidth);
Esitimated_state = legend('error','Location','best');
set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
xlabel('$k$','FontSize',7,'Interpreter','latex') 
ylabel('error','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)

%--------------RELATIVE ERROR--------------%
%-----------------������-----------------%
Error_State = abs(X_RealState-x_EstiState)./abs(X_RealState);
figure;
plot(t, Error_State,'b','linewidth',LineWidth);
Esitimated_state = legend('relative error','Location','best');
set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
xlabel('$k$','FontSize',7,'Interpreter','latex') 
ylabel('Relative Error','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)



%--------------RELATIVE ERROR--------------%
%-----------------������-----------------%
Error_State = abs(X_RealState-x_EstiState)./(abs(X_RealState)+abs(x_EstiState));
figure;
plot(t, Error_State,'b','linewidth',LineWidth);
Esitimated_state = legend('$\sigma$','Location','best');
set(Esitimated_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
xlabel('$k$','FontSize',7,'Interpreter','latex') 
ylabel('$\sigma$','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)

