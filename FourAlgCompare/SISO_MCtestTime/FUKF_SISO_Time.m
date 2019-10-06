%*************************************************************************%
%                  �������޼��������˲������渴��                         %
%_________________________________________________________________________%
%   ���� : 
%   Ŀ�� : fractional order UKF
%   ����ʵ�� :
%               D^{0.7} x_k = 3*sin(2*x_{k-1}) - exp(x_{k-1}) + w_k
%                       y_k = x_k + v_k
%
%   ��� : 
%
%   ��ע : ���棬����ʱ�估����
%
%*************************************************************************%

 clc
 clear
 
 tempp = 50;

 for kkk = 1:tempp

%���沽��
N = 50;

q = 1;                % ϵͳ������ֵ
r = 1;                % ����������ֵ
Q = 0.81;             % ϵͳ�����������
R = 0.25;             % ���������������

 % GL �����¶̼���ԭ��ĳ���
 L = N+1;

% ���� alpha �״ζ�Ӧ��GL����ϵ�� binomial coefficient 
bino_fir = zeros(1,N);       % ΢�ֽ״�Ϊ0.7ʱ GL �����µ�ϵ��
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

I = eye(1,1);                % ���ɵ�λ��

% err_state_FEKF  = zeros(kk_N,N);

X_state_real = zeros(1,N);         % ��ʵ״̬
Z_state_meas = zeros(1,N);         % ʵ�ʹ۲�ֵ

% ����
W_noise = sqrt(Q)*randn(1,N) + q;  % ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  % ��������

x_0  = 0;                          % ��ʼ״̬     
X_state_real(:,1) = x_0;           % ��ʵ״̬��ʼֵ
Z_state_meas(1,1) = V_noise(1,1);  % �������ݳ�ʼֵ

% ϵͳ�������������
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;

for k=2:1:N
    % ����ʵ��״̬
    diff_X_real = f(X_state_real(:,k-1)) + W_noise(1,k-1); 
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_state_real(1,k+1-i);
    end
    X_state_real(:,k) = diff_X_real - rema;

    % ʵ�ʹ۲�ֵ
    Z_state_meas(1,k) = h(X_state_real(:,k)) + V_noise(1,k); 
end

%*************************************************************************%
%-------------------------�޼��������˲������ܲ���------------------------%
%*************************************************************************%

X_state_esti = zeros(1,N);      % ״̬���Ź���ֵ
P_xesti      = cell(1,N);       % ����������

% ��ʼֵ���ã���ʼ������Ϊ�㣩
P_pred_0     = eye(1,1);        % ��ʼԤ�ⷽ����
P_xesti{1,1} = P_pred_0;        % ��ʼ���Ʒ�����

state_dim = 1;
L_sample  = 2 * state_dim +1;

SigmaPoints = zeros(1, L_sample);
SigmaWeight = zeros(1, L_sample);
%GammaPoints = zeros(1, L_sample);
ChiPoints   = zeros(1, L_sample);

tic

for k=2:1:N
    % �ԳƲ���
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_state_esti(:,k-1), P_xesti{1,k-1}); 
    
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
               rema = rema + bino_fir(1,i)*X_state_esti(1,k+1-i);
            end
        else
            for i = 2:1:k
                rema = rema + bino_fir(1,i)*X_state_esti(1,k+1-i);
            end
        end
    X_pre = X_pre - rema;

    % Predicted state error covariance 
    P_xpre = 0*I;
    for i = 1 : 1 : L_sample 
        P_xpre = P_xpre + ...
                 SigmaWeight(:,i)*(ChiPoints(:,i)-X_pre)*(ChiPoints(:,i)-X_pre)' + ...
                 SigmaWeight(:,i)*(SigmaPoints(:,i)-X_state_esti(1,k-1))*(ChiPoints(:,i)-X_pre)'*bino_fir(1,2)' ...
                 + bino_fir(1,2)*SigmaWeight(:,i)*(SigmaPoints(:,i)-X_state_esti(1,k-1))*(ChiPoints(:,i)-X_pre)';
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
    X_state_esti(:,k) = X_pre + Kk*( Z_state_meas(1,k) - Z_pre );

    % estimation error covariance
    P_xesti{1,k} = P_xpre - Kk*P_zpre*Kk';
    
end


 FUKF_SISO_TIME(1,kkk) = toc;
 FUKF_ERROR_norm1(1,kkk) = norm((X_state_real - X_state_esti),1);
 FUKF_ERROR_norm2(1,kkk) = norm((X_state_real - X_state_esti),2);
 FUKF_ERROR_RMSE(1,kkk) = sqrt(sum((X_state_real - X_state_esti).^2)/N);

 end

 FUKF_SISO_ERROR_TIME(1,:) = FUKF_ERROR_norm1(1,:);
 FUKF_SISO_ERROR_TIME(2,:) = FUKF_ERROR_norm2(1,:);
 FUKF_SISO_ERROR_TIME(3,:) = FUKF_SISO_TIME(1,:);
 %FUKF_SISO_ERROR_TIME(4,:) = FUKF_SISO_RMSE(1,:);
 FUKF_ERROR_norm1_average  = sum(FUKF_SISO_ERROR_TIME(1,:))/50
 FUKF_ERROR_norm2_average  = sum(FUKF_SISO_ERROR_TIME(2,:))/50
 FUKF_SISO_TIME_average    = sum(FUKF_SISO_ERROR_TIME(3,:))/50
% FUKF_SISO_RMSE_average    = sum(FUKF_SISO_ERROR_TIME(4,:))/50
%  save FUKF_SISO_ERROE_TIME1 FUKF_SISO_ERROR_TIME FUKF_ERROR_norm1_average ...
%       FUKF_ERROR_norm2_average FUKF_SISO_TIME_average



