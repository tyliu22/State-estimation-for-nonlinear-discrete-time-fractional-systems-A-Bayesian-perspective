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

q = [0 0 0]';                   % 系统噪声均值
r = 0;                          % 测量噪声均值
Q_temp = [0.0001 0.0001 0.0001];
Q = diag(Q_temp);              % 系统噪声方差矩阵
R = 0.0001;                     % 测量噪声方差矩阵

 % GL 定义下短记忆原理的长度
 L = N+1;

alpha = [0.03 1.2 0.1];
 
%计算alpha阶次对应的GL定义系数 binomial coefficient 
bino_fir = zeros(1,N);          %微分阶次为0.03时GL定义下的系数
bino_fir(1,1) = 1;
for i = 2:1:N
    bino_fir(1,i) = (1-(alpha(1,1)+1)/(i-1))*bino_fir(1,i-1);  
end

bino_sec = zeros(1,N);          %微分阶次为1.2时GL定义下的系数
bino_sec(1,1) = 1;
for i = 2:1:N
    bino_sec(1,i) = (1-(alpha(1,2)+1)/(i-1))*bino_sec(1,i-1);  
end

bino_thi = zeros(1,N);          %微分阶次为0.1时GL定义下的系数
bino_thi(1,1) = 1;
for i = 2:1:N
    bino_thi(1,i) = (1-(alpha(1,3)+1)/(i-1))*bino_thi(1,i-1);  
end

%计算GL微分定义下系数矩阵
gamma = cell(1,N);
temp_matrx = zeros(3,3);
for i = 1:1:N 
    temp_matrx(1,1) = bino_fir(1,i);
    temp_matrx(2,2) = bino_sec(1,i);
    temp_matrx(3,3) = bino_thi(1,i);
    gamma{1,i} = temp_matrx;
end

I = eye(3,3);                % 生成单位阵

% err_state_FEKF  = zeros(kk_N,N);

X_state_real = zeros(3,N);         % 真实状态
Z_state_meas = zeros(1,N);         % 实际观测值

% 噪声
W_noise = sqrt(Q)*randn(3,N) + q;  % 系统噪声
V_noise = sqrt(R)*randn(1,N) + r;  % 测量噪声

x_0  = [0 0 0]';                   % 初始状态     
X_state_real(:,1) = x_0;           % 真实状态初始值
Z_state_meas(1,1) = V_noise(1,1);  % 测量数据初始值

% 系统函数与测量函数
f=@(x)[0.2*x(1)*x(3) + 0.3*cos(x(2)) + 0.04; ...
       0.2*cos(x(1)) - 0.5*sin(x(3)) + 0.04; ...
       0.7*sin(x(2)) * cos(x(3)) + 0.04 ];
h=@(x)0.3*x(1)*x(2) - 0.2*sin(x(3));

for k=2:1:N
    % 计算实际状态
    diff_X_real = f(X_state_real(:,k-1)) + W_noise(:,k-1); 
    rema = [0 0 0]';
    for i = 2:1:k
        rema = rema + gamma{1,i}*X_state_real(:,k+1-i);
    end
    X_state_real(:,k) = diff_X_real - rema;

    % 实际观测值
    Z_state_meas(:,k) = h(X_state_real(:,k)) + V_noise(:,k); 
end

% save和load的使用方法：
% save(strcat('D:\存储文件夹名\',filesep,'example','.mat'),'data');
% load(strcat('D:\存储文件夹名\',filesep,'example','.mat'));

% 指定当前路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_FPF\MIMOSystemDataSet\';

% system parameter
SystemParameter = cell(1,6);
SystemParameter{1} = N;     % 仿真步长
SystemParameter{2} = q;     % 系统噪声均值
SystemParameter{3} = r;     % 测量噪声均值
SystemParameter{4} = Q;     % 系统噪声方差
SystemParameter{5} = R;     % 测量噪声方差
SystemParameter{6} = alpha; % 系统阶次

% real measurement
save(strcat(path,'RealMeasurement','.mat'),'Z_state_meas');
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_FPF\MIMOEstimatedState\';
% system parameter
save(strcat(path,'SystemParameter','.mat'),'SystemParameter');
% real state
save(strcat(path,'RealState','.mat'),'X_state_real');



