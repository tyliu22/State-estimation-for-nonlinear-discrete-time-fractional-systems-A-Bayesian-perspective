clc;
clear all;

LineWidth = 1.5;
SimuTimes = 200;    % ����ʱ��
FontSize = 8;

P_alpha_W = 1.6; % alpha�ȶ��ֲ� ����ָ��
P_beta_W  = 0.8; % alpha�ȶ��ֲ� �߶Ȳ���
P_gamma_W = 0.1; % alpha�ȶ��ֲ� ƫб����
P_delta_W = 0;   % alpha�ȶ��ֲ� λ�ò���

P_alpha_V = 1.8; % alpha�ȶ��ֲ� ����ָ��
P_beta_V  = 0.9; % alpha�ȶ��ֲ� �߶Ȳ���
P_gamma_V = 0.2; % alpha�ȶ��ֲ� ƫб����
P_delta_V = 0;   % alpha�ȶ��ֲ� λ�ò���

W_noise = (stblrnd(P_alpha_W,P_beta_W,P_gamma_W,P_delta_W,SimuTimes,1))';  % ϵͳ����
V_noise = (stblrnd(P_alpha_V,P_beta_V,P_gamma_V,P_delta_V,SimuTimes,1))';  % ��������

t = 1 : SimuTimes;
figure;
subplot(211)
plot(t,W_noise,'b-','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 SimuTimes -0.5 1])
ylabel('$\omega$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h1 = legend('$\omega$','Location','best');
set(h1,'Interpreter','latex')

subplot(212)
plot(t,V_noise,'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 SimuTimes -1 1.2])
ylabel('$\nu$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h2 = legend('$\nu$','Location','best');
set(h2,'Interpreter','latex')


%% alpha stable distribution 

% alpha ����ָ�����������ֲ���β��

figure;
x = -5:.01:5;
beta = 0;
gam = 1;
delta = 0;
plot( x , stblpdf(x, .5,  beta, gam, delta, 'quick'), 'b-', ...
      x , stblpdf(x, 1,   beta, gam, delta, 'quick'), 'r--',...
      x , stblpdf(x, 1.5, beta, gam, delta, 'quick'), 'g:',...
      x , stblpdf(x, 2,   beta, gam, delta, 'quick'), 'c-.','linewidth',LineWidth)
axis([-5 5 0 .7]);
title('\alpha-stable densities, \beta = 0, \gamma = 1, \delta = 0');
alpha_state = legend('$\alpha$ = 0.5',...
                         '$\alpha$ = 1.0',...
                         '$\alpha$ = 1.5',...
                         '$\alpha$ = 2', 'Interpreter','latex');
set(alpha_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
% xlabel('$k$','FontSize',7,'Interpreter','latex') 
% ylabel('error','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)


% figure;
% x = -5:.01:5;
% beta = 1;
% gam = 1;
% delta = 0;
% plot( x , stblpdf(x,.5,beta,gam,delta,'quick'),'b-', ...
%   x , stblpdf(x,1,beta,gam,delta,'quick'),'r--',...
%   x , stblpdf(x,1.5,beta,gam,delta,'quick'),'g:',...
%   x , stblpdf(x,2,beta,gam,delta,'quick'),'c-.','linewidth',LineWidth)
% axis([-5 5 0 .6]);
% title('Skewed \alpha-stable densities, \beta = 0.5,\gamma = 1, \delta = 0');
% Skewed_state = legend('$\alpha$ = 0.5',...
%     '$\alpha$ = 1.0',...
%     '$\alpha$ = 1.5',...
%     '$\alpha$ = 2', 'Interpreter','latex');
% set(Skewed_state,'Interpreter','latex')
% set(gcf,'Position',[200 200 400 300]); 
% %axis([0 200 -6 6]) % ������������ָ��������
% axis normal
% set(gca,'FontSize',10); 
% % xlabel('$k$','FontSize',7,'Interpreter','latex') 
% % ylabel('error','FontSize',7);
% % ����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',8)

%%

% beta ƫбָ�����������ֲ�ƫб

figure;
x = -5:.01:5;
alpha = 1;
% beta = 0;
gam = 1;
delta = 0;
plot( x , stblpdf(x, alpha, -1, gam, delta, 'quick'),'b-', ...
      x , stblpdf(x, alpha, 0, gam, delta, 'quick'),'r--',...
      x , stblpdf(x, alpha, 1, gam, delta, 'quick'),'g:','linewidth',LineWidth)
axis([-5 5 0 .5]);
title('\alpha-stable densities, \alpha = 1, \gamma = 1, \delta = 0');
beta_state = legend('$\beta$ = -1',...
    '$\beta$ = 0',...
    '$\beta$ = 1', 'Interpreter','latex');
set(beta_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
% xlabel('$k$','FontSize',7,'Interpreter','latex') 
% ylabel('error','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)


% gamma ��ģָ��

figure;
x = -5:.01:5;
alpha = 1;
beta = 0;
% gam = 1;
delta = 0;
plot( x , stblpdf(x, alpha, beta, 0.5, delta, 'quick'),'b-', ...
      x , stblpdf(x, alpha, beta, 1, delta, 'quick'),'r--',...
      x , stblpdf(x, alpha, beta, 2, delta, 'quick'),'g:',...
      x , stblpdf(x, alpha, beta, 4, delta, 'quick'),'c-.','linewidth',LineWidth)
axis([-5 5 0 .7]);
title('\alpha-stable densities, \alpha = 1, \beta = 0, \delta = 0');
gamma_state = legend('$\gamma$ = 0.5',...
                     '$\gamma$ = 1',...
                     '$\gamma$ = 2',...
                     '$\gamma$ = 4', 'Interpreter','latex');
set(gamma_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
% xlabel('$k$','FontSize',7,'Interpreter','latex') 
% ylabel('error','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)


% delta ��ģָ��

figure;
x = -5:.01:5;
alpha = 1;
beta = 0;
gam = 1;
% delta = 0;
plot( x , stblpdf(x, alpha, beta, gam, -2, 'quick'),'b-', ...
      x , stblpdf(x, alpha, beta, gam, 0, 'quick'),'r--',...
      x , stblpdf(x, alpha, beta, gam, 2, 'quick'),'g:','linewidth',LineWidth)
axis([-5 5 0 .5]);
title('\alpha-stable densities, \alpha = 1, \beta = 0, \gamma = 1');
gamma_state = legend('$\delta$ = 0.5',...
                     '$\delta$ = 1',...
                     '$\delta$ = 2', 'Interpreter','latex');
set(gamma_state,'Interpreter','latex')
set(gcf,'Position',[200 200 400 300]); 
%axis([0 200 -6 6]) % ������������ָ��������
axis normal
set(gca,'FontSize',10); 
% xlabel('$k$','FontSize',7,'Interpreter','latex') 
% ylabel('error','FontSize',7);
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',8)

