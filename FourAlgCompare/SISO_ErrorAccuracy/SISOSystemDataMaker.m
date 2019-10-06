%*************************************************************************%
%                  分数阶无迹卡尔曼滤波器仿真复现                         %
%_________________________________________________________________________%
%   论文 : 
%   目的 : 分数阶无迹卡尔曼滤波器仿真复现
%   函数实验 :
%               D^{0.7} x_k = 3*sin(2*x_{k-1}) -x_{k-1} + w_k
%                       y_k = x_k + v_k
%
%   结果 : 
%
%   备注 : 
%
%*************************************************************************%

clc
clear

% 仿真步长
N = 50;

q = 1;                % 系统噪声均值
r = 1;                % 测量噪声均值
Q = 0.81;               % 系统噪声方差矩阵
R = 0.25;                  % 测量噪声方差矩阵

 % GL 定义下短记忆原理的长度
 L = N+1;

% 计算 alpha 阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);       % 微分阶次为0.7时 GL 定义下的系数
alpha = 0.7;
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha+1)/(i-1))*bino_fir(1,i-1);  
end

I = eye(1,1);                % 生成单位阵

% err_state_FEKF  = zeros(kk_N,N);

X_state_real = zeros(1,N);         % 真实状态
Z_state_meas = zeros(1,N);         % 实际观测值

% 噪声
W_noise = sqrt(Q)*randn(1,N) + q;  % 系统噪声
V_noise = sqrt(R)*randn(1,N) + r;  % 测量噪声

x_0  = 0;                          % 初始状态     
X_state_real(:,1) = x_0;           % 真实状态初始值
Z_state_meas(1,1) = V_noise(1,1);  % 测量数据初始值

% 系统函数与测量函数
f=@(x)3*sin(2*x)-exp(x);
h=@(x)x;

for k=2:1:N
    % 计算实际状态
    diff_X_real = f(X_state_real(:,k-1)) + W_noise(1,k-1); 
    rema = 0;
    for i = 2:1:k
        rema = rema + bino_fir(1,i)*X_state_real(1,k+1-i);
    end
    X_state_real(:,k) = diff_X_real - rema;

    % 实际观测值
    Z_state_meas(1,k) = h(X_state_real(:,k)) + V_noise(1,k); 
end

% save和load的使用方法：
% save(strcat('D:\存储文件夹名\',filesep,'example','.mat'),'data');
% load(strcat('D:\存储文件夹名\',filesep,'example','.mat'));

% 指定当前路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\SystemDataSet\';

% system parameter
SystemParameter = zeros(1,6);
SystemParameter(1) = N;     % 仿真步长
SystemParameter(2) = q;     % 系统噪声均值
SystemParameter(3) = r;     % 测量噪声均值
SystemParameter(4) = Q;     % 系统噪声方差
SystemParameter(5) = R;     % 测量噪声方差
SystemParameter(6) = alpha; % 系统阶次

% real measurement
save(strcat(path,'RealMeasurement','.mat'),'Z_state_meas');
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\SISO_ErrorAccuracy\EstimatedState\';
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');
% real state
save(strcat(path,'RealState','.mat'),'X_state_real');



