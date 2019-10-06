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

q = [0 0 0]';                   % ϵͳ������ֵ
r = 0;                          % ����������ֵ
Q_temp = [0.0001 0.0001 0.0001];
Q = diag(Q_temp);              % ϵͳ�����������
R = 0.0001;                     % ���������������

 % GL �����¶̼���ԭ��ĳ���
 L = N+1;

alpha = [0.03 1.2 0.1];
 
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

I = eye(3,3);                % ���ɵ�λ��

% err_state_FEKF  = zeros(kk_N,N);

X_state_real = zeros(3,N);         % ��ʵ״̬
Z_state_meas = zeros(1,N);         % ʵ�ʹ۲�ֵ

% ����
W_noise = sqrt(Q)*randn(3,N) + q;  % ϵͳ����
V_noise = sqrt(R)*randn(1,N) + r;  % ��������

x_0  = [0 0 0]';                   % ��ʼ״̬     
X_state_real(:,1) = x_0;           % ��ʵ״̬��ʼֵ
Z_state_meas(1,1) = V_noise(1,1);  % �������ݳ�ʼֵ

% ϵͳ�������������
f=@(x)[0.2*x(1)*x(3) + 0.3*cos(x(2)) + 0.04; ...
       0.2*cos(x(1)) - 0.5*sin(x(3)) + 0.04; ...
       0.7*sin(x(2)) * cos(x(3)) + 0.04 ];
h=@(x)0.3*x(1)*x(2) - 0.2*sin(x(3));

for k=2:1:N
    % ����ʵ��״̬
    diff_X_real = f(X_state_real(:,k-1)) + W_noise(:,k-1); 
    rema = [0 0 0]';
    for i = 2:1:k
        rema = rema + gamma{1,i}*X_state_real(:,k+1-i);
    end
    X_state_real(:,k) = diff_X_real - rema;

    % ʵ�ʹ۲�ֵ
    Z_state_meas(:,k) = h(X_state_real(:,k)) + V_noise(:,k); 
end

% save��load��ʹ�÷�����
% save(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'),'data');
% load(strcat('D:\�洢�ļ�����\',filesep,'example','.mat'));

% ָ����ǰ·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_FPF\MIMOSystemDataSet\';

% system parameter
SystemParameter = cell(1,6);
SystemParameter{1} = N;     % ���沽��
SystemParameter{2} = q;     % ϵͳ������ֵ
SystemParameter{3} = r;     % ����������ֵ
SystemParameter{4} = Q;     % ϵͳ��������
SystemParameter{5} = R;     % ������������
SystemParameter{6} = alpha; % ϵͳ�״�

% real measurement
save(strcat(path,'RealMeasurement','.mat'),'Z_state_meas');
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');

% ���ܷ����ļ���·��
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_FPF\MIMOEstimatedState\';
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');
% real state
save(strcat(path,'RealState','.mat'),'X_state_real');



