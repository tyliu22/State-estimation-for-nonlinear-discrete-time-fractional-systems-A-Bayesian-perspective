%*************************************************************************%
%                  �������޼��������˲������渴��                         %
%_________________________________________________________________________%
%   ���� : 
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

clc
clear

% ���沽��
N = 50;

q = 1;                % ϵͳ������ֵ
r = 1;                % ����������ֵ
Q = 0.81;               % ϵͳ�����������
R = 0.25;                  % ���������������

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

% save��load��ʹ�÷�����
% save(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'),'data');
% load(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'));

% ָ����ǰ·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\SystemDataSet\';

% system parameter
SystemParameter = zeros(1,6);
SystemParameter(1) = N;     % ���沽��
SystemParameter(2) = q;     % ϵͳ������ֵ
SystemParameter(3) = r;     % ����������ֵ
SystemParameter(4) = Q;     % ϵͳ��������
SystemParameter(5) = R;     % ������������
SystemParameter(6) = alpha; % ϵͳ�״�

% real measurement
save(strcat(path,'RealMeasurement','.mat'),'Z_state_meas');
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');

% ���ܷ����ļ���·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\EstimatedState\';
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');
% real state
save(strcat(path,'RealState','.mat'),'X_state_real');



