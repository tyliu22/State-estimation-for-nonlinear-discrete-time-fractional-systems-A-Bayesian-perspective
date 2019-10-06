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
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\FourAlgCompare\SISO_ErrorAccuracy\EstimatedState\';
SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));

N = SystemParameter(1);     % ���沽��
q = SystemParameter(2);     % ϵͳ������ֵ
r = SystemParameter(3);     % ����������ֵ
Q = SystemParameter(4);     % ϵͳ��������
R = SystemParameter(5);     % ������������
alpha = SystemParameter(6); % ϵͳ�״�

FPF_X_esti = cell2mat(struct2cell(load(strcat(path,'FPF_EstimatedState','.mat'))));
FUKF_X_esti = cell2mat(struct2cell(load(strcat(path,'FUKF_EstimatedState','.mat'))));
FEKF_X_esti = cell2mat(struct2cell(load(strcat(path,'FEKF_EstimatedState','.mat'))));
FCDKF_X_esti = cell2mat(struct2cell(load(strcat(path,'FCDKF_EstimatedState','.mat'))));

RealState = cell2mat(struct2cell(load(strcat(path,'RealState','.mat'))));





k = 1:1:50;
LineWidth = 1.5;
FontSize = 7;

%% error
% figure;
% subplot(221)
% plot(k,RealState(1,:),'b-',k,FEKF_X_esti(1,:),'r--','linewidth',LineWidth);
% % axis([xmin xmax ymin ymax])  ������������ָ��������
% axis normal
% axis([0 25 -15 5])
% ylabel('$x$','FontSize',FontSize,'Interpreter','latex')
% xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% % ����������̶��������ƣ���С
% set(gcf,'Position',[200 200 400 300]); 
% set(gca,'FontName','Helvetica','FontSize',FontSize)
% legend('real','FEKF','Location','best');
% 
% subplot(222)
% plot(k,RealState(1,:),'b-',k,FCDKF_X_esti(1,:),'r--','linewidth',LineWidth);
% % axis([xmin xmax ymin ymax])  ������������ָ��������
% axis normal
% axis([0 25 -15 5])
% ylabel('$x$','FontSize',FontSize,'Interpreter','latex')
% xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% % ����������̶��������ƣ���С
% set(gcf,'Position',[200 200 400 300]); 
% set(gca,'FontName','Helvetica','FontSize',FontSize)
% legend('real','FCDKF','Location','best');
% 
% subplot(223)
% plot(k,RealState(1,:),'b-',k,FUKF_X_esti(1,:),'r--','linewidth',LineWidth);
% % axis([xmin xmax ymin ymax])  ������������ָ��������
% axis normal
% axis([0 25 -15 5])
% ylabel('$x$','FontSize',FontSize,'Interpreter','latex')
% xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% % ����������̶��������ƣ���С
% set(gcf,'Position',[200 200 400 300]); 
% set(gca,'FontName','Helvetica','FontSize',FontSize)
% legend('real','FUKF','Location','best');
% 
% subplot(224)
% plot(k,RealState(1,:),'b-',k,FPF_X_esti(1,:),'r--','linewidth',LineWidth);
% % axis([xmin xmax ymin ymax])  ������������ָ��������
% axis normal
% axis([0 25 -15 5])
% ylabel('$x$','FontSize',FontSize,'Interpreter','latex')
% xlabel('$k$','FontSize',FontSize,'Interpreter','latex')
% % ����������̶��������ƣ���С
% set(gcf,'Position',[200 200 400 300]); 
% set(gca,'FontName','Helvetica','FontSize',FontSize)
% legend('real','FPF','Location','best');

%% RMSE
% SE
Error_FPF_X_esti   = (FPF_X_esti   - RealState).^2;
Error_FUKF_X_esti  = (FUKF_X_esti  - RealState).^2;
Error_FEKF_X_esti  = (FEKF_X_esti  - RealState).^2;
Error_FCDKF_X_esti = (FCDKF_X_esti - RealState).^2;

RMSE_FPF_X_esti   = zeros(1,N);
RMSE_FUKF_X_esti  = zeros(1,N);
RMSE_FEKF_X_esti  = zeros(1,N);
RMSE_FCDKF_X_esti = zeros(1,N);

RMSE_FPF_X_esti(1,1)   = Error_FPF_X_esti(1,1);
RMSE_FUKF_X_esti(1,1)  = Error_FUKF_X_esti(1,1);
RMSE_FEKF_X_esti(1,1)  = Error_FEKF_X_esti(1,1);
RMSE_FCDKF_X_esti(1,1) = Error_FCDKF_X_esti(1,1);

 for i = 2:1:N
     RMSE_FPF_X_esti(1,i)   = RMSE_FPF_X_esti(1,i-1)   + Error_FPF_X_esti(1,i);
     RMSE_FUKF_X_esti(1,i)  = RMSE_FUKF_X_esti(1,i-1)  + Error_FUKF_X_esti(1,i);
     RMSE_FEKF_X_esti(1,i)  = RMSE_FEKF_X_esti(1,i-1)  + Error_FEKF_X_esti(1,i);
     RMSE_FCDKF_X_esti(1,i) = RMSE_FCDKF_X_esti(1,i-1) + Error_FCDKF_X_esti(1,i);
 end
 
  for i = 1:1:N
     RMSE_FPF_X_esti(1,i)   = sqrt( RMSE_FPF_X_esti(1,i) / i );
     RMSE_FUKF_X_esti(1,i)  = sqrt( RMSE_FUKF_X_esti(1,i) / i );
     RMSE_FEKF_X_esti(1,i)  = sqrt( RMSE_FEKF_X_esti(1,i) / i );
     RMSE_FCDKF_X_esti(1,i) = sqrt( RMSE_FCDKF_X_esti(1,i) / i );
  end
 
figure;
semilogy(k,RMSE_FEKF_X_esti(1,:),'b-',k,RMSE_FCDKF_X_esti(1,:),'r--', ...
     k,RMSE_FUKF_X_esti(1,:),'c:', k,RMSE_FPF_X_esti(1,:),'g-.', 'linewidth',LineWidth);
set(gcf,'Position',[200 200 400 300]); 
% axis([xmin xmax ymin ymax])  ������������ָ��������
axis normal
%axis([0 N 0 1])
ylabel('RMSE','FontSize',FontSize)
xlabel('iteration times','FontSize',FontSize)
% ����������̶��������ƣ���С
set(gca,'FontName','Helvetica','FontSize',FontSize)
legend('FEKF','FCDKF','FUKF','FPF','Location','best');
 
 
 
 





%% Square Error
% figure;
% plot(k,Error_FEKF_X_esti(1,:),'b-',k,Error_FCDKF_X_esti(1,:),'r--', ...
%      k,Error_FUKF_X_esti(1,:),':c', k,Error_FPF_X_esti(1,:),'g-.', 'linewidth',LineWidth);
% set(gcf,'Position',[200 200 400 300]); 
% % axis([xmin xmax ymin ymax])  ������������ָ��������
% axis normal
% axis([0 N 0 4])
% ylabel('SE','FontSize',FontSize)
% xlabel('iteration times','FontSize',FontSize)
% % ����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',FontSize)
% legend('FEKF','FCDKF','FUKF','FPF','Location','best');


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

