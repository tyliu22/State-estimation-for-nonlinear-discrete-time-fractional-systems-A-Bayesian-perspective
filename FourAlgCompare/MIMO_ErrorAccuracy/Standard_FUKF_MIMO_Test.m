%*************************************************************************%
%                  分数阶无迹卡尔曼滤波器仿真复现                         %
%_________________________________________________________________________%
%   论文 : fractional order FUKF
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

% save和load的使用方法：
% save(strcat('D:\存储文件夹名\',filesep,'example','.mat'),'data');
% load(strcat('D:\存储文件夹名\',filesep,'example','.mat'));

% SystemParameter = zeros(1,6);
% SystemParameter(1) = N;     % 仿真步长
% SystemParameter(2) = q;     % 系统噪声均值
% SystemParameter(3) = r;     % 测量噪声均值
% SystemParameter(4) = Q;     % 系统噪声方差
% SystemParameter(5) = R;     % 测量噪声方差
% SystemParameter(6) = alpha; % 系统阶次

clc
clear

% 指定当前路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOSystemDataSet\';
SystemParameter = struct2cell(load(strcat(path,'SystemParameter','.mat')));
% SystemParameter = cell2mat(struct2cell(load(strcat(path,'SystemParameter','.mat'))));
Z_meas = cell2mat(struct2cell(load(strcat(path,'RealMeasurement','.mat'))));

%仿真步长
N = SystemParameter{1,1}{1};     % 仿真步长
q = SystemParameter{1,1}{2};     % 系统噪声均值
r = SystemParameter{1,1}{3};     % 测量噪声均值
Q = SystemParameter{1,1}{4};     % 系统噪声方差
R = SystemParameter{1,1}{5};     % 测量噪声方差
alpha = SystemParameter{1,1}{6}; % 系统阶次

%GL定义下短记忆原理的长度
L = N+1;

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

%系统矩阵设置
I = eye(3,3);                %生成单位阵

% 系统函数与测量函数
f=@(x)[0.2*x(1)*x(3) + 0.3*cos(x(2)) + 0.04; ...
       0.2*cos(x(1)) - 0.5*sin(x(3)) + 0.04; ...
       0.7*sin(x(2)) * cos(x(3)) + 0.04 ];
h=@(x)0.3*x(1)*x(2) - 0.2*sin(x(3));

%*************************************************************************%
%-------------------------无迹卡尔曼滤波器性能测试------------------------%
%*************************************************************************%
X_esti = zeros(3,N);          % 状态最优估计值
X_esti(:,1) = [0.1 0.1 0.1]'; % 初始状态估计
P_xesti = cell(1,N);          % 估计误差方差阵

% 初始值设置（初始矩阵不能为零）
P_pred_0 = [10,0,0;0,10,0;0,0,10];   % 初始预测方差阵
P_xesti{1,1} = P_pred_0;             % 初始估计方差阵

state_dim = 3;
L_sample  = 2 * state_dim +1;

SigmaPoints = zeros(state_dim, L_sample);
SigmaWeight = zeros(1, L_sample);
GammaPoints = zeros(state_dim, L_sample);
ChiPoints   = zeros(1, L_sample);

for k=2:1:N
    % 对称采样
    [SigmaWeight, SigmaPoints] = ukf_sample(state_dim, X_esti(:,k-1), P_xesti{1,k-1}); 

    for i = 1 : 1 : L_sample
        GammaPoints(:,i) = f(SigmaPoints(:,i)) + q;
    end

    % Predicted state
    X_pre = [0; 0; 0];
    for i = 1 : 1 : L_sample
        X_pre = X_pre +  SigmaWeight(:,i)*GammaPoints(:,i);
    end
        % 计算余项
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
    % 计算余项
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


% % 输入与测量输出图
% k = 1:1:N;
% LineWidth = 1.5;
% 
% % square error
% figure;
% plot(k,X_state_real(1,:),'r',k,X_state_esti(1,:),'b--','linewidth',LineWidth);
% % set(gcf,'Position',[200 200 400 300]); 
% % axis([xmin xmax ymin ymax])设置坐标轴在指定的区间
%  axis normal
%  axis([ 0 N -6 6 ])
% ylabel('x','FontSize',8)
% xlabel('time(sec)','FontSize',8)
% % 设置坐标轴刻度字体名称，大小
% set(gca,'FontName','Helvetica','FontSize',8)
% legend('real state','estimated state','Location','best');
% legend('Real state 1','Estimation state 1','Location','best');


FUKF_X_esti = X_esti;

% 性能分析文件夹路径
path = 'C:\software\matlab\matlab_example\KalmanFilter\FKF Framework\MIMO_ErrorAccuracy\MIMOEstimatedState\';
save(strcat(path,'FUKF_EstimatedState','.mat'),'FUKF_X_esti');


