%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �����׿������˲������渴��
%   ���ģ� 
%   Ŀ�ģ�FCDKF, FEKF, FPF, FUKF ���ƾ������ܱȽ�
%         ����ʵ��:    D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                              y_k = x_k + v_k
%   �����
%
%   ��ע��
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save��load��ʹ�÷�����
% save(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'),'data');
% load(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'));

% SystemParameter(1) = N;     % ���沽��
% SystemParameter(2) = q;     % ϵͳ������ֵ
% SystemParameter(3) = r;     % ����������ֵ
% SystemParameter(4) = Q;     % ϵͳ��������
% SystemParameter(5) = R;     % ������������
% SystemParameter(6) = alpha; % ϵͳ�״�

% FPF_EstimatedState   FPF_X_esti 
% FUKF_EstimatedState  FUKF_X_esti 
% FEKF_EstimatedState  FEKF_X_esti
% FCDKF_EstimatedState FCDKF_X_esti

clc
clear

% ���ܷ����ļ���·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_Gau_FPF\MIMOEstimatedState\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));

%���沽��
N = SystemParameter{1,1}{1};     % ���沽��
q = SystemParameter{1,1}{2};     % ϵͳ������ֵ
r = SystemParameter{1,1}{3};     % ����������ֵ
Q = SystemParameter{1,1}{4};     % ϵͳ��������
R = SystemParameter{1,1}{5};     % ������������
alpha = SystemParameter{1,1}{6}; % ϵͳ�״�


FPF_X_esti_500 = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState_500','.mat'))));
FPF_X_esti_250 = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState_250','.mat'))));
FPF_X_esti_100 = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState_100','.mat'))));

RealState = cell2mat(struct2cell(load(strcat(path,'RealState','.mat'))));


% FEKF_error = FEKF_X_esti - RealState;


k = 1:1:50;
LineWidth = 1.5;
FontSize = 7;

%% error
figure;
subplot(311)
plot(k,RealState(1,:),'b-',k,FPF_X_esti_250(1,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 25 -0.2 0.5])
ylabel('$x_1$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h1 = legend('$x_1$','$\hat{x}_1$','Location','best');
set(h1,'Interpreter','latex')

subplot(312)
plot(k,RealState(2,:),'b-',k,FPF_X_esti_250(2,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 25 -1 1.2])
ylabel('$x_2$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h2 = legend('$x_2$','$\hat{x}_2$','Location','best');
set(h2,'Interpreter','latex')

subplot(313)
plot(k,RealState(3,:),'b-',k,FPF_X_esti_250(3,:),'r--','linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 25 -0.5 0.8])
ylabel('$x_3$','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
h3 = legend('$x_3$','$\hat{x}_3$','Location','best');
set(h3,'Interpreter','latex')



%% RMSE
% SE
Error_FPF_X_esti_100 = (FPF_X_esti_100(:,:) - RealState(:,:)).^2; % .^2
Error_FPF_X_esti_250 = (FPF_X_esti_250(:,:) - RealState(:,:)).^2;
Error_FPF_X_esti_500 = (FPF_X_esti_500(:,:) - RealState(:,:)).^2;


%% error
figure;
subplot(311)
plot(k,Error_FPF_X_esti_100(1,:),'b-',k,Error_FPF_X_esti_250(1,:),'r--',k,Error_FPF_X_esti_500(1,:),':c', ...
    'linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 25 0 0.3])
ylabel('SE($x_1$)','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
he1 = legend('$N=10$','$N=50$','$N=100$','Location','best');
set(he1,'Interpreter','latex')

subplot(312)
plot(k,Error_FPF_X_esti_100(2,:),'b-',k,Error_FPF_X_esti_250(2,:),'r--',k,Error_FPF_X_esti_500(2,:),':c', ...
    'linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 25 0 4])
ylabel('SE($x_2$)','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
he2 = legend('$N=10$','$N=50$','$N=100$','Location','best');
set(he2,'Interpreter','latex')

subplot(313)
plot(k,Error_FPF_X_esti_100(3,:),'b-',k,Error_FPF_X_esti_250(3,:),'r--',k,Error_FPF_X_esti_500(3,:),':c', ...
    'linewidth',LineWidth);
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
axis([0 25 0 0.5])
ylabel('SE($x_3$)','FontSize',FontSize,'Interpreter','latex')
xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% ����������̶��������ƣ���С
set(gcf,'Position',[200 200 400 300]); 
set(gca,'FontName','Helvetica','FontSize',FontSize)
he3 = legend('$N=10$','$N=50$','$N=100$','Location','best');
set(he3,'Interpreter','latex')

%% 
% 'r' ��ɫ 'm' �ۺ�
% 'g' ��ɫ 'c' ��ɫ
% 'b' ��ɫ 'w' ��ɫ
% 'y' ��ɫ 'k' ��ɫ

% '-' ʵ�� '--' ����
% ':' ���� '-.' �㻮��

% '.' �õ�Ż��Ƹ����ݵ� '^' �������ǻ��Ƹ����ݵ�
% '+' ��'+'�Ż��Ƹ����ݵ� 'v' �������ǻ��Ƹ����ݵ�
% '*' ��'*'�Ż��Ƹ����ݵ� '>' �������ǻ��Ƹ����ݵ�
% ' .' ��'.'�Ż��Ƹ����ݵ� '<' �������ǻ��Ƹ����ݵ�
% 's'��squar �������λ��Ƹ����ݵ�'p' ������ǻ��Ƹ����ݵ�
% 'd'��diamond�����λ��Ƹ����ݵ� 'h' �������ǻ��Ƹ����ݵ�

