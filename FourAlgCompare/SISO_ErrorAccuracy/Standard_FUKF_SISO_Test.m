%*************************************************************************%
%                  �������޼��������˲������渴��                         %
%_________________________________________________________________________%
%   ���� : fractional order FUKF
%   Ŀ�� : �������޼��������˲������渴��
%   ����ʵ�� :
%               D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                       y_k = x_k + v_k
%
%   ��� : 
%
%   ��ע : 
%
%*************************************************************************%

% save��load��ʹ�÷�����
% save(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'),'data');
% load(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'));

% SystemParameter = zeros(1,6);
% SystemParameter(1) = N;     % ���沽��
% SystemParameter(2) = q;     % ϵͳ������ֵ
% SystemParameter(3) = r;     % ����������ֵ
% SystemParameter(4) = Q;     % ϵͳ��������
% SystemParameter(5) = R;     % ������������
% SystemParameter(6) = alpha; % ϵͳ�״�

clc
clear

% ָ����ǰ·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\SystemDataSet\';
SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));

%���沽��
N = SystemParameter(1);     % ���沽��
q = SystemParameter(2);     % ϵͳ������ֵ
r = SystemParameter(3);     % ����������ֵ
Q = SystemParameter(4);     % ϵͳ��������
R = SystemParameter(5);     % ������������
alpha = SystemParameter(6); % ϵͳ�״�

Z_meas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

 % GL �����¶̼���ԭ��ĳ���
 L = N+1;

% ���� alpha �״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       % ΢�ֽ״�Ϊ alpha ʱ GL �����µ�ϵ��
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

I = eye(1,1);                 % ���ɵ�λ��

X_esti = zeros(1,N);          % ״̬���Ź���ֵ

% ϵͳ�������������
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;

%*************************************************************************%
%-------------------------�޼��������˲������ܲ���------------------------%
%*************************************************************************%

P_xesti = cell(1,N);            % ����������

% ��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0     = 10;              % ��ʼԤ�ⷽ����
P_xesti{1,1} = P_pred_0;        % ��ʼ���Ʒ�����

state_dim = 1;
L_sample  = 2 * state_dim +1;

SigmaPoints = zeros(1, L_sample);
SigmaWeight = zeros(1, L_sample);
%GammaPoints = zeros(1, L_sample);
ChiPoints   = zeros(1, L_sample);

for k=2:1:N
    % �ԳƲ���
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_esti(:,k-1), P_xesti{1,k-1}); 
    
    for i = 1 : 1 : L_sample
        ChiPoints(:,i) = f(SigmaPoints(:,i)) + q;
    end

    % Predicted state
        X_pre = 0;
    for i = 1 : 1 : L_sample
        X_pre = X_pre +  SigmaWeight(:,i)*ChiPoints(:,i);
    end
        % ��������
        rema = 0;
        if k>L
            for i = 2:1:L+1
               rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
            end
        else
            for i = 2:1:k
                rema = rema + bino_fir(1,i)*X_esti(1,k+1-i);
            end
        end
    X_pre = X_pre - rema;

    % Predicted state error covariance 
    P_xpre = 0*I;
    for i = 1 : 1 : L_sample 
        P_xpre = P_xpre + ...
                 SigmaWeight(:,i)*(ChiPoints(:,i)-X_pre)*(ChiPoints(:,i)-X_pre)' + ...
                 SigmaWeight(:,i)*(SigmaPoints(:,i)-X_esti(1,k-1))*(ChiPoints(:,i)-X_pre)'*bino_fir(1,2)' ...
                 + bino_fir(1,2)*SigmaWeight(:,i)*(SigmaPoints(:,i)-X_esti(1,k-1))*(ChiPoints(:,i)-X_pre)';
    end
        % ��������
        rema_P = 0;
        if k>L+1
            for i = 2:1:L+2
                rema_P = rema_P + bino_fir(1,i)*P_xesti{1,k+1-i}*bino_fir(1,i)';
            end
        else
            for i = 2:1:k
                rema_P = rema_P + bino_fir(1,i)*P_xesti{1,k+1-i}*bino_fir(1,i)';
            end
        end
    P_xpre = P_xpre + rema_P + Q;

    % measurement update
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_pre, P_xpre); 
    
    for i = 1 : 1 : L_sample
        ChiPoints(:,i) = h(SigmaPoints(:,i)) + r;
    end
    
    % Predicted measurement
    Z_pre = 0;
    for i = 1 : 1 : L_sample
        Z_pre = Z_pre + SigmaWeight(:,i)*ChiPoints(:,i);
    end
    
    % Predicted measurement error covariance 
    P_zpre = 0*I;
    for i = 1 : 1 : L_sample 
        P_zpre = P_zpre + SigmaWeight(:,i)*(ChiPoints(:,i)-Z_pre)* ...
                 (ChiPoints(:,i)-Z_pre)';
    end
    P_zpre = P_zpre + R;
   
    % cross-variance 
    P_xzpre = 0;
    for i = 1 : 1 : L_sample 
        P_xzpre = P_xzpre + SigmaWeight(:,i)*(SigmaPoints(:,i)-X_pre)* ...
                 (ChiPoints(:,i)-Z_pre)';
    end
    
    % Kalman gain
    Kk = P_xzpre/P_zpre;
    
    % estimated state
    X_esti(:,k) = X_pre + Kk*( Z_meas(1,k) - Z_pre );

    % estimation error covariance
    P_xesti{1,k} = P_xpre - Kk*P_zpre*Kk';

end


% % ������������ͼ
% k = 1:1:N;
% LineWidth = 1.5;
% 
% % square error
% figure;
% plot(k,X_state_real(1,:),'r',k,X_state_esti(1,:),'b--','linewidth',LineWidth);
% % set(gcf,'Position',[200 200 400 300]); 
% % axis([xmin xmax ymin ymax])������������ָ��������
%  axis normal
%  axis([ 0 N -6 6 ])
% ylabel('x','FontSize',8)
% xlabel('time(sec)','FontSize',8)
% % ����������̶��������ƣ���С
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
% legend('Real state 1','Estimation state 1','Location','best');


FUKF_X_esti = X_esti;

% ���ܷ����ļ���·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\EstimatedState\';
save(strcat(path,'FUKF_EstimatedState','.mat'),'FUKF_X_esti');


