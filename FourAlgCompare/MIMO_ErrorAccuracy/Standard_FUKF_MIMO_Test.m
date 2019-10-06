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
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOSystemDataSet\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));
% SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));
Z_meas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

%���沽��
N = SystemParameter{1,1}{1};     % ���沽��
q = SystemParameter{1,1}{2};     % ϵͳ������ֵ
r = SystemParameter{1,1}{3};     % ����������ֵ
Q = SystemParameter{1,1}{4};     % ϵͳ��������
R = SystemParameter{1,1}{5};     % ������������
alpha = SystemParameter{1,1}{6}; % ϵͳ�״�

%GL�����¶̼���ԭ��ĳ���
L = N+1;

%����alpha�״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);          %΢�ֽ״�Ϊ0.03ʱGL�����µ�ϵ��
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha(1,1)+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %΢�ֽ״�Ϊ1.2ʱGL�����µ�ϵ��
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha(1,2)+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,N);          %΢�ֽ״�Ϊ0.1ʱGL�����µ�ϵ��
bino_thi(1,1) = 1;
for i = 2:1:N
    bino_thi(1,i) = (1-(alpha(1,3)+1)/(i-1))*bino_thi(1,i-1);  
end

%����GL΢�ֶ�����ϵ������
gamma = cell(1,N);
temp_matrx = zeros(3,3);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

%ϵͳ��������
I = eye(3,3);                %���ɵ�λ��

% ϵͳ�������������
f=@(x)[0.2*x(1)*x(3) + 0.3*cos(x(2)) + 0.04; ...
       0.2*cos(x(1)) - 0.5*sin(x(3)) + 0.04; ...
       0.7*sin(x(2)) * cos(x(3)) + 0.04 ];
h=@(x)0.3*x(1)*x(2) - 0.2*sin(x(3));

%*************************************************************************%
%-------------------------�޼��������˲������ܲ���------------------------%
%*************************************************************************%
X_esti = zeros(3,N);          % ״̬���Ź���ֵ
X_esti(:,1) = [0.1 0.1 0.1]'; % ��ʼ״̬����
P_xesti = cell(1,N);          % ����������

% ��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0 = [10,0,0;0,10,0;0,0,10];   % ��ʼԤ�ⷽ����
P_xesti{1,1} = P_pred_0;             % ��ʼ���Ʒ�����

state_dim = 3;
L_sample  = 2 * state_dim +1;

SigmaPoints = zeros(state_dim, L_sample);
SigmaWeight = zeros(1, L_sample);
GammaPoints = zeros(state_dim, L_sample);
ChiPoints   = zeros(1, L_sample);

for k=2:1:N
    % �ԳƲ���
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_esti(:,k-1), P_xesti{1,k-1}); 

    for i = 1 : 1 : L_sample
        GammaPoints(:,i) = f(SigmaPoints(:,i)) + q;
    end

    % Predicted state
    X_pre = [0; 0; 0];
    for i = 1 : 1 : L_sample
        X_pre = X_pre +  SigmaWeight(:,i)*GammaPoints(:,i);
    end
        % ��������
        rema = 0;
        if k>L
            for i = 2:1:L+1
               rema = rema + gamma{1,i}*X_esti(:,k+1-i);
            end
        else
            for i = 2:1:k
                rema = rema + gamma{1,i}*X_esti(:,k+1-i);
            end
        end
    X_pre = X_pre - rema;

    % Predicted state error covariance 
    P_xpre = 0*I;
    for i = 1 : 1 : L_sample 
        P_xpre = P_xpre + ...
                 SigmaWeight(:,i)*(GammaPoints(:,i)-X_pre)*(GammaPoints(:,i)-X_pre)' + ...
                 SigmaWeight(:,i)*(SigmaPoints(:,i)-X_esti(:,k-1))*(GammaPoints(:,i)-X_pre)'*gamma{1,2}' ...
                 + ( SigmaWeight(:,i)*(SigmaPoints(:,i)-X_esti(:,k-1))*(GammaPoints(:,i)-X_pre)'*gamma{1,2}' )';
                 %+ gamma{1,2}*( SigmaWeight(:,i)*(SigmaPoints(:,i)-X_esti(:,k-1))*(GammaPoints(:,i)-X_pre)')';
    end
    % ��������
    rema_P = 0;

    for i = 2:1:k
        rema_P = rema_P + gamma{1,i}*P_xesti{1,k+1-i}*gamma{1,i}';
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
    P_zpre = 0;
    for i = 1 : 1 : L_sample 
        P_zpre = P_zpre + SigmaWeight(:,i)*(ChiPoints(:,i)-Z_pre)* ...
                 (ChiPoints(:,i)-Z_pre)';
    end
    P_zpre = P_zpre + R;
   
    % cross-variance 
    P_xzpre = [0; 0; 0];
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
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOEstimatedState\';
save(strcat(path,'FUKF_EstimatedState','.mat'),'FUKF_X_esti');


